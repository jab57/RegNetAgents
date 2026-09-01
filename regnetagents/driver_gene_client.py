"""
driver_gene_client.py -- cancer-driver gene annotation from the IntOGen compendium.

Provides an in-memory lookup of known cancer-driver genes and their consensus
mode-of-action, sourced from IntOGen's Compendium of Mutational Cancer Driver
Genes (release 2024.09.20, CC0 1.0 Universal / public domain).

    Martinez-Jimenez et al., "A compendium of mutational cancer driver genes",
    Nature Reviews Cancer 2020 (doi:10.1038/s41568-020-0290-x)

The reference data ships with the repo as a trimmed, pre-collapsed snapshot at
``regnetagents/reference_data/intogen_drivers.tsv`` -- no network access or setup
step is required. The snapshot has four columns, with ``#``-prefixed provenance
comments at the top:

- ``symbol`` / ``role`` -- the pan-cancer driver call (see below). These two are
  positionally first and never change meaning; :func:`load_driver_roles` and
  :func:`get_driver_roles` read only these.
- ``n_cancer_types`` -- breadth: how many of IntOGen's distinct cancer types
  (max 86 in this release) called the gene a driver. A "how established /
  pan-cancer is this driver" signal; deliberately not ``len(cancer_types)``.
- ``cancer_types`` -- semicolon-separated RegNetAgents TCGA codes (subset of the
  14) the gene was called in, via :data:`INTOGEN_TO_TCGA_CANCER_TYPE`. Presence
  only -- no per-cancer-type role is stored (per-cohort role calls are too noisy
  to trust; e.g. they flag CCND1, a canonical oncogene, as LoF in HNSC). Empty
  for a driver with no RegNetAgents-mappable cancer type. The pan-cancer
  ``role`` stays the sole authoritative role regardless of tissue match.

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

from .tcga_registry import TCGA_CANCER_TYPES

logger = logging.getLogger(__name__)

SNAPSHOT_PATH = Path(__file__).parent / "reference_data" / "intogen_drivers.tsv"

# RegNetAgents TCGA cancer-type codes eligible for tissue-matching (the 14
# committed TCGA networks). Source of truth: regnetagents.tcga_registry.
TCGA_CANCER_TYPE_CODES = frozenset(TCGA_CANCER_TYPES)

# IntOGen CANCER_TYPE code -> RegNetAgents TCGA code(s). Used to BUILD the
# committed snapshot's `cancer_types` column, not at runtime (the snapshot
# already stores mapped codes) -- kept here as the explicit, documented mapping.
# Verified against the 2024.09.20 IntOGen release. Judgment calls:
#   - READ (rectal) -> coad; RegNetAgents has no separate rectal network.
#   - EGC (esophagogastric) -> stad.
#   - NSCLC / LUNG are pooled, non-subtype-specific cohorts -> counted toward
#     BOTH luad and lusc (the call can't be attributed to one subtype).
# Deliberately NOT mapped (distinct tumor biology despite similar organ/name):
#   NPC (vs hnsc), PRCC / CHRCC (papillary/chromophobe renal, vs kirc's clear
#   cell), LIHB (hepatoblastoma, vs lihc's hepatocellular), PANET
#   (neuroendocrine, vs paad's adenocarcinoma), UCS, MGCT.
INTOGEN_TO_TCGA_CANCER_TYPE: dict[str, tuple[str, ...]] = {
    "BLCA": ("blca",), "BLADDER": ("blca",),
    "BRCA": ("brca",),
    "CESC": ("cesc",), "CEAD": ("cesc",),
    "COAD": ("coad",), "COADREAD": ("coad",), "READ": ("coad",),
    "HNSC": ("hnsc",),
    "CCRCC": ("kirc",), "RCC": ("kirc",),
    "HCC": ("lihc",),
    "LUAD": ("luad",), "LUSC": ("lusc",),
    "NSCLC": ("luad", "lusc"), "LUNG": ("luad", "lusc"),
    "OVT": ("ov",),
    "PAAD": ("paad",), "PANCREAS": ("paad",),
    "PRAD": ("prad",), "PROSTATE": ("prad",),
    "STAD": ("stad",), "STOMACH": ("stad",), "EGC": ("stad",),
    "UCEC": ("ucec",),
}

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

# Lazily populated alongside _ROLES_CACHE. {SYMBOL: {tcga_code, ...}} from the
# `cancer_types` column; an empty set means "IntOGen driver, but not called in
# any RegNetAgents TCGA cancer type" (e.g. a blood-cancer-only driver).
_CANCER_TYPES_CACHE: dict[str, set[str]] | None = None


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


def load_driver_cancer_types(path: str | Path = SNAPSHOT_PATH) -> dict[str, set[str]]:
    """Parse the ``cancer_types`` column into ``{SYMBOL_UPPER: {tcga_code, ...}}``.

    Presence only -- no role. An empty cell yields an empty set (the gene is an
    IntOGen driver but was not called in any RegNetAgents TCGA cancer type). A
    snapshot without the column (an older file) yields all-empty sets rather
    than raising. Raises ``FileNotFoundError`` if the snapshot is missing.
    """
    out: dict[str, set[str]] = {}
    with open(path, encoding="utf-8") as fh:
        data_lines = (ln for ln in fh if not ln.lstrip().startswith("#"))
        reader = csv.DictReader(data_lines, delimiter="\t")
        for row in reader:
            symbol = (row.get("symbol") or "").strip().upper()
            if not symbol:
                continue
            raw = (row.get("cancer_types") or "").strip()
            # "".split(";") is [""], not [] -- guard both ways.
            out[symbol] = {c for c in raw.split(";") if c} if raw else set()
    return out


def _ensure_loaded() -> None:
    """Populate the caches on first use; pin to {} permanently on failure."""
    global _ROLES_CACHE, _CANCER_TYPES_CACHE
    if _ROLES_CACHE is not None:
        return
    try:
        # Look up SNAPSHOT_PATH as a module global at call time (not as a
        # bound default argument) so tests can monkeypatch it.
        _ROLES_CACHE = load_driver_roles(SNAPSHOT_PATH)
        _CANCER_TYPES_CACHE = load_driver_cancer_types(SNAPSHOT_PATH)
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
        _CANCER_TYPES_CACHE = {}


def get_driver_roles() -> dict[str, str]:
    """Return the cached ``{SYMBOL_UPPER: role}`` map (empty dict if unavailable)."""
    _ensure_loaded()
    assert _ROLES_CACHE is not None  # set by _ensure_loaded
    return _ROLES_CACHE


def get_driver_cancer_types() -> dict[str, set[str]]:
    """Return the cached ``{SYMBOL_UPPER: {tcga_code, ...}}`` map (empty if unavailable)."""
    _ensure_loaded()
    assert _CANCER_TYPES_CACHE is not None  # set by _ensure_loaded
    return _CANCER_TYPES_CACHE


def is_tissue_matched(gene: str, cancer_type: str) -> bool:
    """True if ``gene`` was called an IntOGen driver specifically in ``cancer_type``.

    ``cancer_type`` is a RegNetAgents TCGA code (e.g. ``"coad"``). This is a pure
    presence check and returns no role. A gene that is not a driver at all, and a
    gene that is a driver only in *other* cancer types, both return ``False`` --
    use :func:`get_driver_roles` for pan-cancer driver/role status. Unknown or
    unmapped ``cancer_type`` values simply return ``False``.
    """
    ct = (cancer_type or "").strip().lower()
    sym = str(gene).strip().upper()
    return ct in get_driver_cancer_types().get(sym, frozenset())


def driver_data_available() -> bool:
    """True once the reference snapshot has loaded with at least one gene.

    Triggers the same lazy load as :func:`get_driver_roles` first, so the answer
    is correct regardless of call order. Lets callers tell "loaded fine, no
    drivers in this gene set" apart from "the reference file failed to load".
    """
    _ensure_loaded()
    return _ROLES_CACHE not in (None, {})


def annotate_genes(genes, cancer_type: str | None = None) -> dict[str, dict]:
    """Annotate an iterable of gene symbols.

    Returns ``{SYMBOL_UPPER: {"is_driver": bool, "role": str | None}}`` where
    ``role`` is one of the collapsed IntOGen roles for drivers, else ``None``.

    When ``cancer_type`` (a RegNetAgents TCGA code) is given, each entry also
    carries ``"tissue_matched": bool`` -- whether IntOGen called the gene a
    driver specifically in that cancer type. ``role`` is unchanged by this: it
    is always the pan-cancer role. Caller is responsible for validating
    ``cancer_type`` (see :data:`TCGA_CANCER_TYPE_CODES`).
    """
    roles = get_driver_roles()
    ct = cancer_type.strip().lower() if cancer_type else None
    ct_map: dict[str, set[str]] = get_driver_cancer_types() if ct is not None else {}
    out: dict[str, dict] = {}
    for g in genes:
        sym = str(g).strip().upper()
        if not sym:
            continue
        role = roles.get(sym)
        entry: dict = {"is_driver": role is not None, "role": role}
        if ct is not None:
            entry["tissue_matched"] = ct in ct_map.get(sym, frozenset())
        out[sym] = entry
    return out
