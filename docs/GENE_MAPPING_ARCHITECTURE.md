# Gene ID Mapping Architecture

## Overview
The RegNetAgents system uses **Ensembl Gene IDs** (ENSG format) internally, but users can query with familiar **gene symbols** (like "APC", "TP53"). The mapping system automatically converts between formats.

## Storage Layers

### 1. Local Cache (Primary)
- **File**: `cache/gene_id_cache.pkl`
- **Format**: Python pickle (binary)
- **Content**: Bidirectional mapping dictionaries (`symbol_to_ensembl`, `ensembl_to_symbol`)
- **Persistence**: Saved between sessions
- **Speed**: Instant lookup (no network calls)
- **Pre-populated**: `scripts/build_network_cache.py --enrich-gene-cache` bulk-resolves all
  GREmLN ENSG IDs via MyGene.info at cache build time (~100% coverage, 29,000+ entries)

### 2. Ensembl REST API (Fallback)
- **URL**: `https://rest.ensembl.org/lookup/`
- **Used when**: Gene not in local cache
- **Speed**: ~1-2 seconds per gene
- **Auto-caching**: Results saved to local cache
- **Reliability**: Official Ensembl database

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

## Alternative Flow (New Gene)

1. **User Input**: `"BRCA2"` (not in cache)
2. **API Query**: `https://rest.ensembl.org/lookup/symbol/homo_sapiens/BRCA2`
3. **Cache Result**: Save to `gene_id_cache.pkl`
4. **Query RegNetAgents**: Using retrieved Ensembl ID
5. **Return Results**: Include both formats

## Cache Management

- **Auto-saving**: Every successful API lookup
- **Bulk enrichment**: Run `python scripts/build_network_cache.py --enrich-gene-cache` to
  pre-populate all GREmLN ENSG IDs via MyGene.info (done automatically with `--all`)
- **Error handling**: Graceful fallback if API fails
- **Cache location**: `cache/gene_id_cache.pkl` (project root)
- **Growth**: Expands as new genes are queried at runtime

## Benefits

1. **Fast lookups**: Local cache = instant results
2. **Network resilient**: Works offline for cached genes
3. **User-friendly**: Accept familiar gene symbols
4. **Automatic**: No manual mapping needed
5. **Persistent**: Cache survives server restarts

## Files Involved

- `gene_id_mapper.py` - Main mapping class
- `gene_id_cache.pkl` - Local cache storage
- `regnetagents_langgraph_mcp_server.py` - Integration with MCP server
- `regnetagents_langgraph_workflow.py` - Workflow implementation

## API Endpoints Used

- Symbol → Ensembl: `GET /lookup/symbol/homo_sapiens/{symbol}`
- Ensembl → Symbol: `GET /lookup/id/{ensembl_id}`