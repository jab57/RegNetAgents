import pandas as pd

# Step 1: Load the epithelial cell ARACNe network file
path = "C:/Dev/RegNetAgents/models/networks/epithelial_cell/network.tsv"
network = pd.read_csv(path, sep="\t")

# Step 2: Extract all unique Ensembl IDs
genes = set(network['regulator.values']).union(set(network['target.values']))
print(f"Total unique Ensembl IDs in epithelial network: {len(genes)}")

# Step 3: Define CRC-associated genes (Ensembl IDs)
crc_ensembl_ids = {
    "ENSG00000134982",  # APC
    "ENSG00000141510",  # TP53
    "ENSG00000133703",  # KRAS
    "ENSG00000213281",  # NRAS
    "ENSG00000157764",  # BRAF
    "ENSG00000141646",  # SMAD4
    "ENSG00000121879",  # PIK3CA
    "ENSG00000076242",  # MLH1
    "ENSG00000116141",  # MSH2
    "ENSG00000116062",  # MSH6
    "ENSG00000148737",  # TCF7L2
    "ENSG00000168646",  # AXIN2
    "ENSG00000109685",  # FBXW7
    "ENSG00000171862",  # PTEN
    "ENSG00000039068"   # CDH1
}

# Step 4: Compare directly
present = [eid for eid in crc_ensembl_ids if eid in genes]
missing = [eid for eid in crc_ensembl_ids if eid not in genes]

# Step 5: Print results with gene symbols
lookup = {
    "ENSG00000134982":"APC","ENSG00000141510":"TP53","ENSG00000133703":"KRAS","ENSG00000213281":"NRAS",
    "ENSG00000157764":"BRAF","ENSG00000141646":"SMAD4","ENSG00000121879":"PIK3CA","ENSG00000076242":"MLH1",
    "ENSG00000116141":"MSH2","ENSG00000116062":"MSH6","ENSG00000148737":"TCF7L2","ENSG00000168646":"AXIN2",
    "ENSG00000109685":"FBXW7","ENSG00000171862":"PTEN","ENSG00000039068":"CDH1"
}

print("\n=== Results ===")
print("CRC genes found in epithelial network:")
for eid in present:
    print(f"  {eid} ({lookup[eid]})")

print("\nCRC genes missing:")
for eid in missing:
    print(f"  {eid} ({lookup[eid]})")

print("\nScript finished successfully.")
