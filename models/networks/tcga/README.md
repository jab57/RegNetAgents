# TCGA tumor-state networks (not included)

The 14 TCGA ARACNe networks used by RegNetAgents come from the Bioconductor package
[`aracne.networks`](https://bioconductor.org/packages/aracne.networks) (Giorgi FM,
Alvarez MJ). Its authors publish the network files on Zenodo
(https://doi.org/10.5281/zenodo.22918956) under CC BY-NC-ND 4.0: attribution, non-commercial use only, and no sharing of
modified versions. RegNetAgents therefore does not ship the networks.

Install them locally (downloads the files from Zenodo; you accept their license):

```bash
pip install -e ".[tcga]"
python scripts/setup_tcga_networks.py --accept-license
```

This writes `<type>/network.csv` and `<type>/network_index.pkl` for blca, brca, cesc, coad,
hnsc, kirc, lihc, luad, lusc, ov, paad, prad, stad and ucec, and checks each network against
a recorded checksum, so every install matches the networks RegNetAgents was built and
evaluated with. Everything else in RegNetAgents (GREmLN cell-type networks, driver
annotation) works without them.
