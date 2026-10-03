#!/usr/bin/env python3
"""
Install the TCGA tumor-state networks used by RegNetAgents.

The TCGA ARACNe networks come from the Bioconductor package ``aracne.networks``
(author Federico M. Giorgi), which is distributed under a Columbia University
software evaluation license (non-commercial academic research use only; no
redistribution). RegNetAgents therefore does not ship these networks. This script
downloads the package from Bioconductor -- so you obtain it from the official
source under its own license -- and builds the network caches locally.

What it does:
  1. Shows the license notice and requires --accept-license.
  2. Downloads aracne.networks from Bioconductor (or uses --tarball).
  3. Checks each network's data file against the SHA-256 recorded below
     (network data are identical in aracne.networks 1.36.0 and 1.38.0).
  4. Converts Entrez IDs to gene symbols with a frozen mapping
     (scripts/data/tcga_entrez_to_symbol.json.gz), writes
     models/networks/tcga/<type>/network.csv and checks its SHA-256, so every
     install reproduces exactly the networks used by RegNetAgents and its paper.
  5. Builds models/networks/tcga/<type>/network_index.pkl.

Usage:
    pip install -e ".[tcga]"
    python scripts/setup_tcga_networks.py --accept-license
    python scripts/setup_tcga_networks.py --accept-license --cancer-type brca coad
    python scripts/setup_tcga_networks.py --accept-license --tarball aracne.networks_1.38.0.tar.gz
"""

import argparse
import csv
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

LICENSE_URL = ("https://bioconductor.org/packages/release/data/experiment/"
               "licenses/aracne.networks/LICENSE")
DOWNLOAD_URLS = [
    "https://bioconductor.org/packages/release/data/experiment/src/contrib/aracne.networks_1.38.0.tar.gz",
    "https://bioconductor.org/packages/3.23/data/experiment/src/contrib/aracne.networks_1.38.0.tar.gz",
    "https://bioconductor.org/packages/3.22/data/experiment/src/contrib/aracne.networks_1.36.0.tar.gz",
]
MAP_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "tcga_entrez_to_symbol.json.gz")

LICENSE_NOTICE = f"""
The TCGA networks come from the Bioconductor package aracne.networks, distributed
under a Columbia University software evaluation license. In summary (read the full
text before continuing): use is limited to non-commercial academic or educational
research; you may not redistribute the package or make it available to third
parties; commercial use requires a license from Columbia University.

Full license: {LICENSE_URL}

RegNetAgents does not redistribute these networks. By passing --accept-license you
confirm that you have read the license and that your use complies with it.
"""


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
    from extract_tcga_networks import CANCER_TYPE_MAP, RDA_NAMES, load_rda_from_tarball, regulon_to_edges

    ap = argparse.ArgumentParser(description="Install the TCGA networks from Bioconductor aracne.networks.")
    ap.add_argument("--accept-license", action="store_true",
                    help="Confirm you have read and comply with the aracne.networks license.")
    ap.add_argument("--tarball", help="Use a local aracne.networks_*.tar.gz instead of downloading.")
    ap.add_argument("--cancer-type", nargs="+", choices=list(CANCER_TYPE_MAP), metavar="TYPE",
                    help=f"Only these types (default: all 14): {', '.join(CANCER_TYPE_MAP)}")
    ap.add_argument("--output-dir", default="models/networks/tcga")
    ap.add_argument("--keep-tarball", action="store_true", help="Do not delete a downloaded tarball.")
    args = ap.parse_args()

    if args.tarball:
        args.tarball = os.path.abspath(args.tarball)
    os.chdir(REPO_ROOT)  # the cache builder uses repo-relative paths

    print(LICENSE_NOTICE)
    if not args.accept_license:
        sys.exit("Not accepted: re-run with --accept-license once you have read the license.")

    try:
        import rdata  # noqa: F401
    except ImportError:
        sys.exit('ERROR: this script needs the "rdata" package: pip install -e ".[tcga]"')

    maps = load_symbol_maps()
    types = args.cancer_type or list(CANCER_TYPE_MAP)

    tmpdir = None
    tarball = args.tarball
    if not tarball:
        tmpdir = tempfile.mkdtemp(prefix="aracne_networks_")
        tarball = download(os.path.join(tmpdir, "aracne.networks.tar.gz"))

    from build_tcga_cache import build_tcga_cache

    failures = []
    try:
        with tarfile.open(tarball, "r:gz") as tf:
            rda_sha = {}
            for ct in types:
                member = tf.extractfile(f"aracne.networks/data/{RDA_NAMES[ct]}")
                rda_sha[ct] = sha256(member.read()) if member else ""
        for ct in types:
            expected = maps["types"][ct]
            print(f"\n[{ct}] checking source data ...")
            if rda_sha[ct] != expected["rda_sha256"]:
                failures.append(ct)
                print(f"  ERROR: {RDA_NAMES[ct]} differs from the version RegNetAgents was built "
                      "with; skipping (results would not match).")
                continue
            regulon = load_rda_from_tarball(tarball, ct)
            edges, _ = regulon_to_edges(regulon, symbol_map_for(maps, ct))
            csv_path = os.path.join(args.output_dir, ct, "network.csv")
            written = write_csv(edges, csv_path)
            if sha256(written) != expected["csv_sha256"]:
                failures.append(ct)
                print(f"  ERROR: rebuilt {csv_path} does not match the expected checksum.")
                continue
            print(f"  {len(edges):,} edges -> {csv_path} (checksum verified)")
            build_tcga_cache(ct, output_dir=args.output_dir, skip_validation=True)
    finally:
        if tmpdir and not args.keep_tarball:
            try:
                os.remove(tarball)
                os.rmdir(tmpdir)
            except OSError:
                pass

    if failures:
        sys.exit(f"\nFAILED for: {', '.join(failures)}")
    print(f"\nTCGA networks installed: {', '.join(types)}")


if __name__ == "__main__":
    main()
