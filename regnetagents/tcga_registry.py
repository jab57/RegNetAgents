"""
TCGA ARACNe network registry for RegNetAgents.

Maps the 14 supported TCGA cancer types (all epithelial-origin) to their
network CSV file paths and human-readable labels. GBM and LAML are
intentionally excluded — no reference network of the appropriate cell
lineage exists in RegNetAgents for these cancer types.

To rebuild CSVs from source, use the Bioconductor aracne.networks tarball:

    python scripts/extract_tcga_networks.py \
        --tarball /tmp/aracne.networks.tar.gz \
        --output-dir models/networks/tcga
    python scripts/build_tcga_cache.py --all

See docs/DATA_SOURCES.md for full instructions.
"""

TCGA_NETWORK_REGISTRY = {
    "blca": {
        "label": "Bladder Urothelial Carcinoma",
        "csv": "models/networks/tcga/blca/network.csv",
    },
    "brca": {
        "label": "Breast Invasive Carcinoma",
        "csv": "models/networks/tcga/brca/network.csv",
    },
    "cesc": {
        "label": "Cervical Squamous Cell Carcinoma",
        "csv": "models/networks/tcga/cesc/network.csv",
    },
    "coad": {
        "label": "Colon Adenocarcinoma",
        "csv": "models/networks/tcga/coad/network.csv",
    },
    "hnsc": {
        "label": "Head/Neck Squamous Cell Carcinoma",
        "csv": "models/networks/tcga/hnsc/network.csv",
    },
    "kirc": {
        "label": "Kidney Renal Clear Cell Carcinoma",
        "csv": "models/networks/tcga/kirc/network.csv",
    },
    "lihc": {
        "label": "Liver Hepatocellular Carcinoma",
        "csv": "models/networks/tcga/lihc/network.csv",
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
    "paad": {
        "label": "Pancreatic Adenocarcinoma",
        "csv": "models/networks/tcga/paad/network.csv",
    },
    "prad": {
        "label": "Prostate Adenocarcinoma",
        "csv": "models/networks/tcga/prad/network.csv",
    },
    "stad": {
        "label": "Stomach Adenocarcinoma",
        "csv": "models/networks/tcga/stad/network.csv",
    },
    "ucec": {
        "label": "Uterine Corpus Endometrial Carcinoma",
        "csv": "models/networks/tcga/ucec/network.csv",
    },
    # GBM excluded: no glial/CNS reference network in RegNetAgents
    # LAML excluded: monocyte networks are too sparse / biologically mismatched
}

TCGA_CANCER_TYPES = list(TCGA_NETWORK_REGISTRY.keys())
