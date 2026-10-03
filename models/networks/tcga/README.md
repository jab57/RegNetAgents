# TCGA tumor-state networks (not included)

The 14 TCGA ARACNe networks used by RegNetAgents come from the Bioconductor package
[`aracne.networks`](https://bioconductor.org/packages/aracne.networks) (author Federico M.
Giorgi). That package is distributed under a Columbia University software evaluation
license: non-commercial academic research use only, and no redistribution. RegNetAgents
therefore does not ship the networks.

Install them locally (downloads the package from Bioconductor; you accept its license):

```bash
pip install -e ".[tcga]"
python scripts/setup_tcga_networks.py --accept-license
```

This writes `<type>/network.csv` and `<type>/network_index.pkl` for blca, brca, cesc, coad,
hnsc, kirc, lihc, luad, lusc, ov, paad, prad, stad and ucec, and checks each network against
a recorded checksum, so every install matches the networks RegNetAgents was built and
evaluated with. Everything else in RegNetAgents (GREmLN cell-type networks, driver
annotation) works without them.
