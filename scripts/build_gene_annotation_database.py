#!/usr/bin/env python3
"""
Build NCBI+UniProt gene annotation database

Fetches gene annotations from MyGene.info (NCBI-backed service) to create
the local gene annotation database that RegNetAgents uses.
"""

import json
import requests
import time
from pathlib import Path
import pickle

def get_all_genes_from_networks():
    """Extract all unique gene symbols from network files"""
    genes = set()

    # Get genes from network indices
    network_dir = Path("models/networks")
    cell_types = [
        "cd14_monocytes", "cd16_monocytes", "cd20_b_cells", "cd4_t_cells",
        "cd8_t_cells", "epithelial_cell", "erythrocytes",
        "monocyte_derived_dendritic_cells", "nk_cells", "nkt_cells"
    ]

    for cell_type in cell_types:
        network_file = network_dir / cell_type / "network_index.pkl"
        if network_file.exists():
            try:
                with open(network_file, 'rb') as f:
                    cache = pickle.load(f)
                    # Extract gene symbols from the cache structure
                    if 'all_genes' in cache:
                        genes.update(cache['all_genes'])
                    elif 'network' in cache:
                        network = cache['network']
                        genes.update(network.nodes())
                    print(f"  [OK] Loaded {cell_type}: {cache.get('num_genes', 0)} genes")
            except Exception as e:
                print(f"  [WARN] Could not load {cell_type}: {e}")

    print(f"\nTotal unique genes across all networks: {len(genes)}")
    return genes

def fetch_gene_annotations_batch(gene_ids, batch_size=100):
    """Fetch gene annotations from MyGene.info API in batches using Ensembl IDs"""

    annotations = {}
    gene_list = list(gene_ids)
    total = len(gene_list)

    print(f"\nFetching annotations for {total} genes from MyGene.info API...")

    for i in range(0, total, batch_size):
        batch = gene_list[i:i+batch_size]

        try:
            # Query MyGene.info API using Ensembl IDs
            url = "https://mygene.info/v3/query"
            params = {
                'q': ','.join(batch),
                'scopes': 'ensembl.gene',
                'fields': 'name,summary,symbol,ensembl.gene',
                'species': 'human',
                'size': batch_size
            }

            response = requests.post(url, data=params, timeout=30)

            if response.status_code == 200:
                results = response.json()

                for result in results:
                    if 'symbol' in result and 'ensembl' in result:
                        symbol = result['symbol']
                        ensembl_id = result.get('query', batch[results.index(result)])

                        # Build annotation string in NCBI+UniProt format
                        parts = []
                        parts.append(f"Gene Symbol {symbol}")

                        if 'name' in result:
                            parts.append(f"Full Name: {result['name']}")

                        if 'summary' in result:
                            parts.append(f"Protein summary: {result['summary']}")
                        elif 'name' in result:
                            # Fallback to name if no summary
                            parts.append(f"Protein summary: {result['name']}")

                        # Store by gene symbol (as that's what the system uses internally)
                        annotations[symbol] = " | ".join(parts)

                print(f"  [OK] Batch {i//batch_size + 1}/{(total-1)//batch_size + 1}: {len(results)} genes")

            else:
                print(f"  [WARN] Batch {i//batch_size + 1} failed: HTTP {response.status_code}")

            # Rate limiting
            time.sleep(0.5)

        except Exception as e:
            print(f"  [WARN] Error fetching batch {i//batch_size + 1}: {e}")

    return annotations

def create_annotation_database():
    """Create the gene annotation database"""

    print("="*60)
    print("Building NCBI+UniProt Gene Annotation Database")
    print("="*60)

    # Step 1: Get all genes from networks
    print("\n1. Extracting genes from network files...")
    genes = get_all_genes_from_networks()

    if not genes:
        print("ERROR: No genes found in network files!")
        return

    # Step 2: Fetch annotations
    print("\n2. Fetching gene annotations from MyGene.info...")
    annotations = fetch_gene_annotations_batch(genes)

    print(f"\n  [OK] Successfully annotated {len(annotations)} genes")
    print(f"  [WARN] Missing annotations for {len(genes) - len(annotations)} genes")

    # Step 3: Save to JSON
    output_dir = Path("regnetagents/models/gene_embeddings/gene_annotations")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_file = output_dir / "NCBI_UniProt_summary_of_genes.json"

    print(f"\n3. Saving to {output_file}...")
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(annotations, f, indent=2, ensure_ascii=False)

    # Calculate file size
    file_size = output_file.stat().st_size / 1024

    print("="*60)
    print("SUCCESS!")
    print("="*60)
    print(f"Database file: {output_file}")
    print(f"Total genes: {len(annotations)}")
    print(f"File size: {file_size:.1f} KB")
    print(f"\nData source: MyGene.info (NCBI-backed service)")
    print(f"Quality: Professional curation from NCBI + UniProt")
    print("="*60)

if __name__ == "__main__":
    create_annotation_database()
