"""
TCGA ARACNe network registry for RegNetAgents.

Maps the 8 supported TCGA cancer types (all epithelial-origin) to their
network CSV file paths and human-readable labels. GBM and LAML are
intentionally excluded — no reference network of the appropriate cell
lineage exists in RegNetAgents for these cancer types.

Download CSVs from Figshare (https://figshare.com/s/5d1ffd9f8b2e86e37ed6)
and place them at the paths listed below, then run:

    python scripts/build_tcga_cache.py --all
"""

TCGA_NETWORK_REGISTRY = {
    "brca": {
        "label": "Breast Invasive Carcinoma",
        "csv": "models/networks/tcga/brca/network.csv",
    },
    "coad": {
        "label": "Colon Adenocarcinoma",
        "csv": "models/networks/tcga/coad/network.csv",
    },
    "hnsc": {
        "label": "Head/Neck Squamous Cell Carcinoma",
        "csv": "models/networks/tcga/hnsc/network.csv",
    },
    "luad": {
        "label": "Lung Adenocarcinoma",
        "csv": "models/networks/tcga/luad/network.csv",
    },
    "lusc": {
        "label": "Lung Squamous Cell Carcinoma",
        "csv": "models/networks/tcga/lusc/network.csv",
    },
    "ov": {
        "label": "Ovarian Carcinoma",
        "csv": "models/networks/tcga/ov/network.csv",
    },
    "prad": {
        "label": "Prostate Adenocarcinoma",
        "csv": "models/networks/tcga/prad/network.csv",
    },
    "ucec": {
        "label": "Uterine Corpus Endometrial Carcinoma",
        "csv": "models/networks/tcga/ucec/network.csv",
    },
    # GBM excluded: no glial/CNS reference network in RegNetAgents
    # LAML excluded: monocyte networks are too sparse / biologically mismatched
}

TCGA_CANCER_TYPES = list(TCGA_NETWORK_REGISTRY.keys())
