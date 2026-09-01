"""
driver_gene_client.py -- cancer-driver gene annotation from the IntOGen compendium.

Provides an in-memory lookup of known cancer-driver genes and their consensus
mode-of-action, sourced from IntOGen's Compendium of Mutational Cancer Driver
Genes (release 2024.09.20, CC0 1.0 Universal / public domain).

    Martinez-Jimenez et al., "A compendium of mutational cancer driver genes",
    Nature Reviews Cancer 2020 (doi:10.1038/s41568-020-0290-x)

The reference data ships with the repo as a trimmed, pre-collapsed snapshot at
``regnetagents/reference_data/intogen_drivers.tsv`` -- no network access or setup
step is required. The snapshot has two columns, ``symbol`` and ``role``, with
``#``-prefixed provenance comments at the top.

Role values (one per gene): ``oncogene``, ``tumor_suppressor``, ``mixed``
(genes with an exact Act/LoF split across cohorts), ``ambiguous`` (genes IntOGen
only ever called ambiguous). Every gene present in the snapshot is a driver
regardless of its role value -- treat ``role`` as advisory directional context
and presence/absence as the reliable signal. For genes IntOGen detected in only
one or two cohorts the directional call can be unreliable (it is a plain
majority vote of per-cohort calls); this is inherent to a mutational
positive-selection compendium.

Coverage caveat: IntOGen identifies drivers by positive selection in somatic
point mutations, so it under-covers genes driven mainly by fusion, copy-number
alteration, or epigenetic silencing (which curated resources such as OncoKB /
COSMIC CGC include). A gene absent from this annotation is not a
positive-selection driver in IntOGen's compendium -- not necessarily "not a
cancer driver".

Snapshot regeneration (routine maintenance, not automated): download the current
release from https://www.intogen.org/download, and from
``Compendium_Cancer_Genes.tsv`` collapse the per-cohort ``ROLE`` calls
(``Act`` / ``LoF`` / ``ambiguous``) to one value per ``SYMBOL`` by majority vote
-- ``Act`` majority -> ``oncogene``, ``LoF`` majority -> ``tumor_suppressor``,
exact tie -> ``mixed``, only-ambiguous -> ``ambiguous`` -- then write the
two-column snapshot with refreshed provenance comments.
"""

from __future__ import annotations

import csv
import logging
from pathlib import Path

logger = logging.getLogger(__name__)

SNAPSHOT_PATH = Path(__file__).parent / "reference_data" / "intogen_drivers.tsv"

CITATION = (
    "Martinez-Jimenez et al., A compendium of mutational cancer driver genes, "
    "Nature Reviews Cancer 2020 (doi:10.1038/s41568-020-0290-x)"
)

# Valid collapsed role values (see module docstring).
ROLE_ONCOGENE = "oncogene"
ROLE_TUMOR_SUPPRESSOR = "tumor_suppressor"
ROLE_MIXED = "mixed"
ROLE_AMBIGUOUS = "ambiguous"
VALID_ROLES = frozenset(
    {ROLE_ONCOGENE, ROLE_TUMOR_SUPPRESSOR, ROLE_MIXED, ROLE_AMBIGUOUS}
)

# Lazily populated by get_driver_roles(). None = not loaded yet; {} = load was
# attempted and failed (permanent for this process -- no per-call retries).
_ROLES_CACHE: dict[str, str] | None = None


def load_driver_roles(path: str | Path = SNAPSHOT_PATH) -> dict[str, str]:
    """Parse the trimmed IntOGen snapshot into ``{SYMBOL_UPPER: role}``.

    The snapshot is already collapsed to one role per gene at build time, so this
    is a straight read of the two-column ``symbol``/``role`` TSV, skipping
    ``#`` provenance comments. Any unrecognized role value is normalized to
    ``"ambiguous"`` so the gene still annotates as a driver rather than being
    dropped. Raises ``FileNotFoundError`` if the snapshot is missing.
    """
    roles: dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        data_lines = (ln for ln in fh if not ln.lstrip().startswith("#"))
        reader = csv.DictReader(data_lines, delimiter="\t")
        for row in reader:
            symbol = (row.get("symbol") or "").strip().upper()
            if not symbol:
                continue
            role = (row.get("role") or "").strip().lower()
            if role not in VALID_ROLES:
                role = ROLE_AMBIGUOUS
            roles[symbol] = role
    return roles


def _ensure_loaded() -> None:
    """Populate _ROLES_CACHE on first use; pin to {} permanently on failure."""
    global _ROLES_CACHE
    if _ROLES_CACHE is not None:
        return
    try:
        # Look up SNAPSHOT_PATH as a module global at call time (not as a
        # bound default argument) so tests can monkeypatch it.
        _ROLES_CACHE = load_driver_roles(SNAPSHOT_PATH)
        logger.info(
            "Loaded IntOGen driver annotation: %d genes from %s",
            len(_ROLES_CACHE),
            SNAPSHOT_PATH,
        )
    except Exception as exc:  # missing/corrupt file -- degrade, don't crash
        logger.warning(
            "IntOGen driver annotation unavailable (%s: %s); "
            "driver fields will be reported as unknown.",
            type(exc).__name__,
            exc,
        )
        _ROLES_CACHE = {}


def get_driver_roles() -> dict[str, str]:
    """Return the cached ``{SYMBOL_UPPER: role}`` map (empty dict if unavailable)."""
    _ensure_loaded()
    assert _ROLES_CACHE is not None  # set by _ensure_loaded
    return _ROLES_CACHE


def driver_data_available() -> bool:
    """True once the reference snapshot has loaded with at least one gene.

    Triggers the same lazy load as :func:`get_driver_roles` first, so the answer
    is correct regardless of call order. Lets callers tell "loaded fine, no
    drivers in this gene set" apart from "the reference file failed to load".
    """
    _ensure_loaded()
    return _ROLES_CACHE not in (None, {})


def annotate_genes(genes) -> dict[str, dict]:
    """Annotate an iterable of gene symbols.

    Returns ``{SYMBOL_UPPER: {"is_driver": bool, "role": str | None}}`` where
    ``role`` is one of the collapsed IntOGen roles for drivers, else ``None``.
    """
    roles = get_driver_roles()
    out: dict[str, dict] = {}
    for g in genes:
        sym = str(g).strip().upper()
        if not sym:
            continue
        role = roles.get(sym)
        out[sym] = {"is_driver": role is not None, "role": role}
    return out
