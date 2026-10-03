"""Tests for scripts/setup_tcga_networks.py (offline; no download, no rdata needed)."""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

import setup_tcga_networks as S  # noqa: E402
from regnetagents.tcga_registry import TCGA_CANCER_TYPES  # noqa: E402


def test_requires_accept_license_and_never_downloads(monkeypatch):
    def no_download(*_args, **_kwargs):
        raise AssertionError("download attempted without --accept-license")

    monkeypatch.setattr(S, "download", no_download)
    monkeypatch.setattr(sys, "argv", ["setup_tcga_networks.py"])
    with pytest.raises(SystemExit) as exc:
        S.main()
    assert "accept-license" in str(exc.value)


def test_license_notice_points_to_license():
    assert "evaluation license" in S.LICENSE_NOTICE
    assert S.LICENSE_URL in S.LICENSE_NOTICE


def test_frozen_symbol_map_covers_all_cancer_types():
    data = S.load_symbol_maps()
    assert set(data["types"]) == set(TCGA_CANCER_TYPES)
    for ct, t in data["types"].items():
        assert len(t["rda_sha256"]) == 64 and len(t["csv_sha256"]) == 64, ct
        m = S.symbol_map_for(data, ct)
        assert len(m) > 19000, ct
        assert all(isinstance(k, str) and k.isdigit() for k in list(m)[:100]), ct


def test_frozen_symbol_map_known_genes():
    data = S.load_symbol_maps()
    brca = S.symbol_map_for(data, "brca")
    assert brca["7157"] == "TP53"
    assert brca["4609"] == "MYC"
    assert brca["9668"] == "ZNF432"
