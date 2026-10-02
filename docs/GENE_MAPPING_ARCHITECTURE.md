# Gene ID Mapping Architecture

## Overview
The GREmLN cell-type networks use **Ensembl Gene IDs** (ENSG format) as node IDs, but users can query with familiar **gene symbols** (like "APC", "TP53"). `GeneIDMapper` (`regnetagents/gene_id_mapper.py`) converts between the two using a local cache only — no network call is made for symbol↔ID lookups.

TCGA cancer-type networks are symbol-native (nodes keyed by uppercase gene symbol), so TCGA queries do not go through the mapper.

## Storage

### Local Cache
- **File**: `cache/gene_id_cache.pkl` (tracked in the repository)
- **Format**: Python pickle (binary)
- **Content**: Bidirectional mapping dictionaries (`symbol_to_ensembl`, `ensembl_to_symbol`)
- **Speed**: Instant lookup (no network calls)
- **Built by**: `scripts/build_network_cache.py --enrich-gene-cache` bulk-resolves all
  GREmLN ENSG IDs via MyGene.info at cache build time (done automatically with `--all`)

### Placeholder IDs
`GeneIDMapper` also pre-populates the cache with symbols from the local UniProt gene
database. A symbol with no known Ensembl ID gets a synthetic placeholder ID of the form
`ENSG_CACHED_<symbol>`. Placeholders are never network node IDs; they only let the
symbol be recognized as a real gene.

### Load-time repair
On load, `GeneIDMapper` points `symbol_to_ensembl` at the real Ensembl ID wherever
`ensembl_to_symbol` has one, replacing any placeholder (in memory; existing real
mappings are kept, and a symbol with several real IDs takes the lowest, for
determinism). This runs before placeholders are added, so **every GREmLN network gene
resolves by symbol**. Cache rebuilds (`update_gene_id_cache()` in
`scripts/build_network_cache.py`) likewise replace placeholders with real IDs.

> Versions before this repair stored placeholders for ~9,400 of the 14,621 GREmLN
> genes (~64%), so symbol queries for those genes wrongly reported
> "Gene not found in network". Upgrade if you see that for a gene you expect to be
> present.

## Cache Structure

```python
{
    "symbol_to_ensembl": {
        "APC": "ENSG00000134982",
        "TP53": "ENSG00000141510",
        "BRCA1": "ENSG00000012048",
        "MYC": "ENSG00000136997",
        "GAPDH": "ENSG00000111640"
    },
    "ensembl_to_symbol": {
        "ENSG00000134982": "APC",
        "ENSG00000141510": "TP53",
        "ENSG00000012048": "BRCA1",
        "ENSG00000136997": "MYC",
        "ENSG00000111640": "GAPDH"
    }
}
```

## Usage Flow

1. **User Input**: `"APC"` (gene symbol)
2. **Check Local Cache**: Found in `gene_id_cache.pkl`
3. **Return Ensembl ID**: `"ENSG00000134982"`
4. **Query RegNetAgents Network**: Using Ensembl ID
5. **Return Results**: Include both symbol and Ensembl ID

## Alias Fallback (Unrecognized Symbol)

If a symbol is not found, the workflow calls `GeneIDMapper.resolve_aliases()`, which
queries MyGene.info for the canonical HGNC symbol and known aliases (5 s timeout,
in-memory cache only) and retries the lookup with those. This handles older HGNC
symbols and synonyms; results are not written to `gene_id_cache.pkl`.

## Cache Management

- **Bulk enrichment**: Run `python scripts/build_network_cache.py --enrich-gene-cache` to
  pre-populate all GREmLN ENSG IDs via MyGene.info (done automatically with `--all`)
- **Cache location**: `cache/gene_id_cache.pkl` (project root)
- **Offline**: All symbol↔ID lookups work offline; only the alias fallback needs network access

## Files Involved

- `regnetagents/gene_id_mapper.py` - Main mapping class
- `cache/gene_id_cache.pkl` - Local cache storage
- `scripts/build_network_cache.py` - Builds and enriches the cache
- `regnetagents_langgraph_mcp_server.py` - Integration with MCP server
- `regnetagents_langgraph_workflow.py` - Workflow implementation (alias fallback)
