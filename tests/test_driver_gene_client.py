"""
Tests for regnetagents.driver_gene_client (cancer-driver annotation).

The reference data (regnetagents/reference_data/intogen_drivers.tsv) is committed
to the repo, so these tests need no network access and no skip-guard for a
missing file in the normal case.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from regnetagents import driver_gene_client as dgc


# ---------------------------------------------------------------------------
# load_driver_roles / snapshot parsing
# ---------------------------------------------------------------------------

def test_load_driver_roles_parses_committed_snapshot():
    roles = dgc.load_driver_roles()
    assert isinstance(roles, dict)
    assert len(roles) > 500  # ~633 genes in the 2024.09.20 release
    assert all(v in dgc.VALID_ROLES for v in roles.values())


def test_known_oncogene_role():
    # KRAS: 84 Act calls, 0 LoF calls across cohorts -- unambiguous oncogene.
    roles = dgc.load_driver_roles()
    assert roles["KRAS"] == dgc.ROLE_ONCOGENE


def test_known_tumor_suppressor_role():
    # TP53: 161 LoF calls vs. 51 Act calls -- majority-vote tumor suppressor.
    roles = dgc.load_driver_roles()
    assert roles["TP53"] == dgc.ROLE_TUMOR_SUPPRESSOR


def test_symbol_churn_gene_present_under_current_symbol():
    # KMT2D (formerly MLL2) -- IntOGen uses the current HGNC symbol.
    roles = dgc.load_driver_roles()
    assert roles["KMT2D"] == dgc.ROLE_TUMOR_SUPPRESSOR
    assert "MLL2" not in roles  # confirms IntOGen is on current-symbol side


def test_non_driver_gene_absent():
    roles = dgc.load_driver_roles()
    assert "ACTB" not in roles
    assert "GAPDH" not in roles


def test_missing_file_raises():
    import pytest
    with pytest.raises(FileNotFoundError):
        dgc.load_driver_roles("does/not/exist.tsv")


# ---------------------------------------------------------------------------
# get_driver_roles / driver_data_available (cache semantics)
# ---------------------------------------------------------------------------

def test_driver_data_available_true_in_normal_case():
    assert dgc.driver_data_available() is True


def test_get_driver_roles_returns_cached_dict():
    roles = dgc.get_driver_roles()
    assert roles.get("KRAS") == dgc.ROLE_ONCOGENE


def test_driver_data_available_false_when_snapshot_missing(monkeypatch, tmp_path):
    # Reset the module cache and point at a nonexistent file; monkeypatch
    # restores both attributes after the test, so the real cache is unaffected.
    monkeypatch.setattr(dgc, "_ROLES_CACHE", None)
    monkeypatch.setattr(dgc, "SNAPSHOT_PATH", tmp_path / "missing.tsv")

    assert dgc.driver_data_available() is False
    assert dgc.get_driver_roles() == {}


def test_driver_data_available_triggers_lazy_load_regardless_of_call_order(monkeypatch):
    # Calling driver_data_available() before get_driver_roles() must not report
    # "unavailable" just because the cache hasn't been warmed yet.
    monkeypatch.setattr(dgc, "_ROLES_CACHE", None)
    assert dgc.driver_data_available() is True


# ---------------------------------------------------------------------------
# annotate_genes
# ---------------------------------------------------------------------------

def test_annotate_genes_mixed_set():
    result = dgc.annotate_genes(["KRAS", "TP53", "ACTB", "kras"])
    assert result["KRAS"] == {"is_driver": True, "role": dgc.ROLE_ONCOGENE}
    assert result["TP53"] == {"is_driver": True, "role": dgc.ROLE_TUMOR_SUPPRESSOR}
    assert result["ACTB"] == {"is_driver": False, "role": None}
    # case-insensitive, and repeated symbol collapses to one entry
    assert len(result) == 3


def test_annotate_genes_empty_input():
    assert dgc.annotate_genes([]) == {}


def test_annotate_genes_with_cancer_type_adds_tissue_matched():
    result = dgc.annotate_genes(["KRAS", "MYC", "ACTB"], cancer_type="coad")
    # KRAS is called by IntOGen in colorectal -> tissue-matched
    assert result["KRAS"] == {"is_driver": True, "role": dgc.ROLE_ONCOGENE, "tissue_matched": True}
    # MYC is an IntOGen driver, but only in blood cancers -> not tissue-matched for coad
    assert result["MYC"]["is_driver"] is True
    assert result["MYC"]["tissue_matched"] is False
    # ACTB isn't a driver at all
    assert result["ACTB"] == {"is_driver": False, "role": None, "tissue_matched": False}


def test_annotate_genes_without_cancer_type_omits_tissue_matched():
    result = dgc.annotate_genes(["KRAS"])
    assert "tissue_matched" not in result["KRAS"]


# ---------------------------------------------------------------------------
# tissue matching (cancer_types column / is_tissue_matched)
# ---------------------------------------------------------------------------

def test_load_driver_cancer_types_parses_column():
    m = dgc.load_driver_cancer_types()
    assert len(m) > 500
    assert all(isinstance(v, set) for v in m.values())
    # every code is a valid RegNetAgents TCGA code
    all_codes = set().union(*m.values()) if m else set()
    assert all_codes <= dgc.TCGA_CANCER_TYPE_CODES


def test_tissue_match_three_states_for_coad():
    # (a) tissue-matched, (b) driver elsewhere only, (c) not a driver
    assert dgc.is_tissue_matched("KRAS", "coad") is True           # (a)
    assert dgc.is_tissue_matched("MYC", "coad") is False           # (b) driver, but blood-cancer only
    assert dgc.get_driver_roles().get("MYC") is not None           #     ...still a known driver
    assert dgc.is_tissue_matched("ACTB", "coad") is False          # (c)
    assert dgc.get_driver_roles().get("ACTB") is None              #     ...and not a driver


def test_tissue_match_uses_intogen_to_tcga_mapping():
    # VHL is IntOGen-called under CCRCC/RCC, which map to kirc; exact-code
    # matching on "KIRC" would miss it.
    assert dgc.is_tissue_matched("VHL", "kirc") is True
    assert "CCRCC" in dgc.INTOGEN_TO_TCGA_CANCER_TYPE
    assert dgc.INTOGEN_TO_TCGA_CANCER_TYPE["CCRCC"] == ("kirc",)


def test_tissue_match_case_insensitive():
    assert dgc.is_tissue_matched("kras", "COAD") is True


def test_tissue_match_unknown_cancer_type_is_false_not_error():
    assert dgc.is_tissue_matched("KRAS", "melanoma") is False
    assert dgc.is_tissue_matched("KRAS", "") is False


def test_blood_cancer_only_driver_has_empty_cancer_types():
    # MYC: IntOGen driver (blood cancers), no RegNetAgents-mappable cancer type.
    m = dgc.load_driver_cancer_types()
    assert m["MYC"] == set()


def test_empty_cancer_types_cell_parses_to_empty_set_not_singleton(tmp_path):
    snapshot = tmp_path / "custom.tsv"
    snapshot.write_text(
        "# header comment\n"
        "symbol\trole\tn_cancer_types\tcancer_types\n"
        "FOO\toncogene\t3\t\n"
        "BAR\ttumor_suppressor\t5\tcoad;brca\n",
        encoding="utf-8",
    )
    m = dgc.load_driver_cancer_types(snapshot)
    assert m["FOO"] == set()          # not {""}
    assert m["BAR"] == {"coad", "brca"}


# ---------------------------------------------------------------------------
# blank/unrecognized role handling
# ---------------------------------------------------------------------------

def test_unrecognized_role_normalizes_to_ambiguous(tmp_path):
    snapshot = tmp_path / "custom.tsv"
    snapshot.write_text(
        "# test snapshot\nsymbol\trole\nFOO\tsome_unexpected_value\nBAR\t\n",
        encoding="utf-8",
    )
    roles = dgc.load_driver_roles(snapshot)
    assert roles["FOO"] == dgc.ROLE_AMBIGUOUS
    assert roles["BAR"] == dgc.ROLE_AMBIGUOUS
