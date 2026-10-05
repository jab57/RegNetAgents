# Supplementary tables

| File | Content |
|---|---|
| `table_s1_brca_candidates.csv` | Source-labeled IntOGen-overlapping regulator candidates, BRCA focal genes |
| `table_s2_coad_candidates.csv` | Source-labeled IntOGen-overlapping regulator candidates, COAD focal genes |
| `table_s3_tier_specificity.csv` | Focal genes vs random genes, per network tier and background |

All three are written by `scripts/experiment_regulatory_candidates.py`.

## Data licensing

These tables are not covered by the MIT license of the RegNetAgents source code.

- **TCGA network values** (the `TCGA-only` / `Both` regulator–target relationships and the
  `moa` column in Tables S1–S2) come from the TCGA ARACNe networks of `aracne.networks`
  (Giorgi FM, Alvarez MJ; Zenodo doi:10.5281/zenodo.22918956), licensed
  [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/): attribution,
  non-commercial use only, no distribution of modified versions. Those values remain under that
  license.
- **GREmLN network values** (`GREmLN-only` / `Both` regulators) come from the GREmLN project's
  pre-computed ARACNe networks (Zhang et al., 2025), obtained from the CZI Virtual Cells
  Platform for scientific research use and derived from CELLxGENE Census data
  ([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); CZI Cell Science Program et al.,
  *Nucleic Acids Research* 2025, doi:10.1093/nar/gkae1142).
- **Driver annotations** (`intogen_role`) come from the IntOGen Compendium of Mutational Cancer
  Driver Genes, release 2024.09.20 ([CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/);
  Martínez-Jiménez et al., *Nature Reviews Cancer* 2020).
