import pandas as pd

# Step 1: Load the epithelial cell ARACNe network file
path = "C:/Dev/RegNetAgents/models/networks/epithelial_cell/network.tsv"
network = pd.read_csv(path, sep="\t")

# Step 2: Extract all unique Ensembl IDs
genes = set(network['regulator.values']).union(set(network['target.values']))
print(f"Total unique Ensembl IDs in epithelial network: {len(genes)}")

# Step 3: Define cervical cancer–associated genes (Ensembl IDs)
cervical_ensembl_ids = {
    "ENSG00000141510",  # TP53
    "ENSG00000171862",  # PTEN
    "ENSG00000121879",  # PIK3CA
    "ENSG00000157764",  # BRAF
    "ENSG00000133703",  # KRAS
    "ENSG00000155657",  # RB1
    "ENSG00000139618",  # BRCA1
    "ENSG00000139687",  # BRCA2
    "ENSG00000157764",  # BRAF
    "ENSG00000141579",  # MYC
    "ENSG00000146648",  # EGFR
    "ENSG00000141736",  # CCND1
    "ENSG00000147889",  # CDKN2A (p16)
    "ENSG00000121879",  # PIK3CA
    "ENSG00000171862"   # PTEN
}

# Step 4: Compare directly
present = [eid for eid in cervical_ensembl_ids if eid in genes]
missing = [eid for eid in cervical_ensembl_ids if eid not in genes]

# Step 5: Print results with gene symbols
lookup = {
    "ENSG00000141510":"TP53","ENSG00000171862":"PTEN","ENSG00000121879":"PIK3CA","ENSG00000157764":"BRAF",
    "ENSG00000133703":"KRAS","ENSG00000155657":"RB1","ENSG00000139618":"BRCA1","ENSG00000139687":"BRCA2",
    "ENSG00000141579":"MYC","ENSG00000146648":"EGFR","ENSG00000141736":"CCND1","ENSG00000147889":"CDKN2A"
}

print("\n=== Results ===")
print("Cervical cancer genes found in epithelial network:")
for eid in present:
    print(f"  {eid} ({lookup.get(eid, 'Unknown')})")

print("\nCervical cancer genes missing:")
for eid in missing:
    print(f"  {eid} ({lookup.get(eid, 'Unknown')})")

print("\nScript finished successfully.")
