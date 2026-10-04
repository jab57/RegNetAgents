#!/usr/bin/env python3
"""
Install the TCGA tumor-state networks used by RegNetAgents.

The TCGA ARACNe networks are the data of the Bioconductor package
``aracne.networks`` (Giorgi & Alvarez). Its authors publish the network files on
Zenodo (https://doi.org/10.5281/zenodo.22918956) under CC BY-NC-ND 4.0
(attribution, non-commercial use, no sharing of modified versions). RegNetAgents
does not ship these networks: this script downloads them from the official source
and builds the network caches locally.

What it does:
  1. Shows the license notice and requires --accept-license.
  2. Downloads the network files from the Zenodo record (default), or reads them
     from a Bioconductor aracne.networks 1.36.0/1.38.0 tarball (--source
     bioconductor or --tarball), which carries the package's Columbia license.
  3. Checks each network's data file against the SHA-256 recorded below
     (identical in the Zenodo record and aracne.networks 1.36.0 and 1.38.0).
  4. Converts Entrez IDs to gene symbols with a frozen mapping
     (scripts/data/tcga_entrez_to_symbol.json.gz), writes
     models/networks/tcga/<type>/network.csv and checks its SHA-256, so every
     install reproduces exactly the networks used by RegNetAgents and its paper.
  5. Builds models/networks/tcga/<type>/network_index.pkl.

Usage:
    pip install -e ".[tcga]"
    python scripts/setup_tcga_networks.py --accept-license
    python scripts/setup_tcga_networks.py --accept-license --cancer-type brca coad
    python scripts/setup_tcga_networks.py --accept-license --source bioconductor
    python scripts/setup_tcga_networks.py --accept-license --tarball aracne.networks_1.38.0.tar.gz
"""

import argparse
import csv
import io
import gzip
import hashlib
import json
import os
import sys
import tarfile
import tempfile
import urllib.request

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPTS_DIR)
sys.path.insert(0, SCRIPTS_DIR)
sys.path.insert(0, REPO_ROOT)

ZENODO_RECORD = "https://doi.org/10.5281/zenodo.22918956"
ZENODO_FILE_URL = "https://zenodo.org/records/22918956/files/{name}?download=1"
ZENODO_LICENSE_URL = "https://creativecommons.org/licenses/by-nc-nd/4.0/"
BIOC_LICENSE_URL = ("https://bioconductor.org/packages/release/data/experiment/"
                    "licenses/aracne.networks/LICENSE")
# Version-pinned Bioconductor URLs (the generic "release/" URL changes at every
# Bioconductor release). Network data are identical in 1.36.0 and 1.38.0; from 1.39.x
# the package no longer contains the .rda data files.
DOWNLOAD_URLS = [
    "https://bioconductor.org/packages/3.23/data/experiment/src/contrib/aracne.networks_1.38.0.tar.gz",
    "https://bioconductor.org/packages/3.22/data/experiment/src/contrib/aracne.networks_1.36.0.tar.gz",
]
MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "tcga_entrez_to_symbol.json.gz")

LICENSE_NOTICE = {
    "zenodo": f"""
The TCGA networks are the data files of the Bioconductor package aracne.networks
(Giorgi FM, Alvarez MJ), published by its authors on Zenodo ({ZENODO_RECORD})
under the Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 license
(CC BY-NC-ND 4.0). In summary (read the full text before continuing): credit the
authors and the record; non-commercial use only; do not share modified versions of
the networks (the caches this script builds are for your own use).

Full license: {ZENODO_LICENSE_URL}

RegNetAgents does not redistribute these networks. By passing --accept-license you
confirm that you have read the license and that your use complies with it.
""",
    "bioconductor": f"""
The TCGA networks come from the Bioconductor package aracne.networks, distributed
under a Columbia University software evaluation license. In summary (read the full
text before continuing): use is limited to non-commercial academic or educational
research; you may not redistribute the package or make it available to third
parties; commercial use requires a license from Columbia University. The same
network data are also available under CC BY-NC-ND 4.0 (--source zenodo, the default).

Full license: {BIOC_LICENSE_URL}

RegNetAgents does not redistribute these networks. By passing --accept-license you
confirm that you have read the license and that your use complies with it.
""",
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_symbol_maps() -> dict:
    with gzip.open(MAP_PATH, "rt", encoding="utf-8") as fh:
        data = json.load(fh)
    return data


def symbol_map_for(data: dict, cancer_type: str) -> dict:
    t = data["types"][cancer_type]
    unmapped = set(t["unmapped"])
    m = {e: s for e, s in data["base"].items() if e not in unmapped}
    m.update(t["override"])
    return m


def download(dest: str) -> str:
    for url in DOWNLOAD_URLS:
        try:
            print(f"Downloading {url} ...")
            urllib.request.urlretrieve(url, dest)
            return dest
        except Exception as exc:  # try the next mirror/version
            print(f"  failed: {exc}")
    sys.exit("ERROR: could not download aracne.networks. Download it manually from "
             "https://bioconductor.org/packages/aracne.networks and pass --tarball.")


def download_zenodo_rda(name: str) -> bytes:
    url = ZENODO_FILE_URL.format(name=name)
    print(f"Downloading {url} ...")
    with urllib.request.urlopen(url) as resp:
        return resp.read()


def parse_rda(raw: bytes, cancer_type: str):
    import rdata

    parsed = rdata.read_rda(io.BytesIO(raw))
    # The .rda exports one variable: regulon{ct} (e.g. regulonbrca)
    return parsed.get(f"regulon{cancer_type}", next(iter(parsed.values())))


def write_csv(edges: list, path: str) -> bytes:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["Regulator", "Target", "MoA", "Likelihood"],
                                lineterminator="\r\n")
        writer.writeheader()
        writer.writerows(edges)
    with open(path, "rb") as fh:
        return fh.read()


def main() -> None:
    from extract_tcga_networks import CANCER_TYPE_MAP, RDA_NAMES, regulon_to_edges

    ap = argparse.ArgumentParser(description="Install the TCGA networks (aracne.networks data).")
    ap.add_argument("--accept-license", action="store_true",
                    help="Confirm you have read and comply with the license of the chosen source.")
    ap.add_argument("--source", choices=["zenodo", "bioconductor"], default="zenodo",
                    help="Zenodo record (CC BY-NC-ND 4.0, default) or Bioconductor tarball "
                         "(Columbia license).")
    ap.add_argument("--tarball", help="Use a local aracne.networks_*.tar.gz (implies --source bioconductor).")
    ap.add_argument("--cancer-type", nargs="+", choices=list(CANCER_TYPE_MAP), metavar="TYPE",
                    help=f"Only these types (default: all 14): {', '.join(CANCER_TYPE_MAP)}")
    ap.add_argument("--output-dir", default="models/networks/tcga")
    ap.add_argument("--keep-tarball", action="store_true", help="Do not delete a downloaded tarball.")
    args = ap.parse_args()

    if args.tarball:
        args.tarball = os.path.abspath(args.tarball)
        args.source = "bioconductor"
    os.chdir(REPO_ROOT)  # the cache builder uses repo-relative paths

    print(LICENSE_NOTICE[args.source])
    if not args.accept_license:
        sys.exit("Not accepted: re-run with --accept-license once you have read the license.")

    try:
        import rdata  # noqa: F401
    except ImportError:
        sys.exit('ERROR: this script needs the "rdata" package: pip install -e ".[tcga]"')

    maps = load_symbol_maps()
    types = args.cancer_type or list(CANCER_TYPE_MAP)

    from build_tcga_cache import build_tcga_cache

    # Source data: the .rda bytes of each requested network.
    raw = {}
    tmpdir = None
    tarball = args.tarball
    try:
        if args.source == "zenodo":
            for ct in types:
                try:
                    raw[ct] = download_zenodo_rda(RDA_NAMES[ct])
                except Exception as exc:
                    sys.exit(f"ERROR: could not download {RDA_NAMES[ct]} from Zenodo ({exc}). "
                             "Retry, or use --source bioconductor.")
        else:
            if not tarball:
                tmpdir = tempfile.mkdtemp(prefix="aracne_networks_")
                tarball = download(os.path.join(tmpdir, "aracne.networks.tar.gz"))
            with tarfile.open(tarball, "r:gz") as tf:
                names = set(tf.getnames())
                missing = [ct for ct in types if f"aracne.networks/data/{RDA_NAMES[ct]}" not in names]
                if missing:
                    sys.exit("ERROR: this aracne.networks tarball has no network data files "
                             f"({', '.join(missing)}). Use version 1.36.0 or 1.38.0, e.g.\n  "
                             + DOWNLOAD_URLS[0] + "\nor use --source zenodo.")
                for ct in types:
                    member = tf.extractfile(f"aracne.networks/data/{RDA_NAMES[ct]}")
                    raw[ct] = member.read() if member else b""
    finally:
        if tmpdir and not args.keep_tarball:
            try:
                os.remove(tarball)
                os.rmdir(tmpdir)
            except OSError:
                pass

    failures = []
    for ct in types:
        expected = maps["types"][ct]
        print(f"\n[{ct}] checking source data ...")
        if sha256(raw[ct]) != expected["rda_sha256"]:
            failures.append(ct)
            print(f"  ERROR: {RDA_NAMES[ct]} differs from the version RegNetAgents was built "
                  "with; skipping (results would not match).")
            continue
        regulon = parse_rda(raw[ct], ct)
        edges, _ = regulon_to_edges(regulon, symbol_map_for(maps, ct))
        csv_path = os.path.join(args.output_dir, ct, "network.csv")
        written = write_csv(edges, csv_path)
        if sha256(written) != expected["csv_sha256"]:
            failures.append(ct)
            print(f"  ERROR: rebuilt {csv_path} does not match the expected checksum.")
            continue
        print(f"  {len(edges):,} edges -> {csv_path} (checksum verified)")
        build_tcga_cache(ct, output_dir=args.output_dir, skip_validation=True)

    if failures:
        sys.exit(f"\nFAILED for: {', '.join(failures)}")
    print(f"\nTCGA networks installed: {', '.join(types)}")


if __name__ == "__main__":
    main()
