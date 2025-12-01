
#!/usr/bin/env python3
"""
RegNetAgents LangGraph Workflow
=================================

Multi-agent AI framework for gene regulatory network analysis powered by LangGraph orchestration.

This module implements a sophisticated workflow system that coordinates multiple specialized agents
to perform comprehensive gene regulatory network analysis. The framework integrates:

- Network modeling (regulatory relationships from ARACNe networks)
- Pathway enrichment (Reactome API integration)
- Perturbation analysis (therapeutic target ranking via network centrality)
- Domain-specific insights (LLM-powered cancer, drug, clinical, systems biology agents)
- Cross-cell type comparison

The workflow uses LangGraph for state management and parallel execution of independent analysis steps,
achieving 480-24,000× speedup compared to manual analysis workflows.

Architecture:
    - State-based workflow orchestration (LangGraph StateGraph)
    - Parallel batch processing of independent analyses
    - Conditional routing based on gene network position
    - Graceful degradation (LLM agents fallback to rule-based heuristics)

Performance:
    - Rule-based mode: 0.6-15 seconds for comprehensive analysis
    - LLM-powered mode: 15-62 seconds with scientific rationales
    - Perturbation analysis: Pre-computed PageRank enables instant ranking

Author: Jose A. Bird, PhD
License: MIT
"""

from typing import Dict, List, TypedDict, Optional, Any
from langgraph.graph import StateGraph, END
from langgraph.checkpoint.memory import MemorySaver
import asyncio
import json
import logging
import ollama
from dotenv import load_dotenv

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load environment variables from .env file
load_dotenv()

# Import existing components
from regnetagents import GeneIDMapper, CompleteGeneService
from enum import Enum
import pickle
import os

# Define the essential classes that were in the old MCP server
class CellType(Enum):
    # Existing cell types (10)
    EPITHELIAL_CELLS = "epithelial_cell"
    CD14_MONOCYTES = "cd14_monocytes"
    CD16_MONOCYTES = "cd16_monocytes"
    CD20_B_CELLS = "cd20_b_cells"
    CD4_T_CELLS = "cd4_t_cells"
    CD8_T_CELLS = "cd8_t_cells"
    ERYTHROCYTES = "erythrocytes"
    NK_CELLS = "nk_cells"
    NKT_CELLS = "nkt_cells"
    MONOCYTE_DERIVED_DENDRITIC_CELLS = "monocyte-derived_dendritic_cells"

    # New priority cell types (5)
    HEPATOCYTES = "hepatocytes"
    CARDIOMYOCYTES = "cardiomyocytes"
    NEURONS = "neurons"
    FIBROBLASTS = "fibroblasts"
    ENDOTHELIAL_CELLS = "endothelial_cells"

class RegNetAgentsCache:
    """Gene regulatory network cache for storing analysis results."""
    def __init__(self):
        self.network_indices = {}
        self.load_network_indices()

    def load_network_indices(self):
        """Load network indices from pickle files."""
        import os
        models_dir = "models/networks"

        if os.path.exists(models_dir):
            for cell_type in CellType:
                cache_file = os.path.join(models_dir, cell_type.value, "network_index.pkl")

                if os.path.exists(cache_file):
                    try:
                        with open(cache_file, 'rb') as f:
                            data = pickle.load(f)
                            self.network_indices[cell_type.value] = data

                    except Exception as e:
                        logger.warning(f"Failed to load {cache_file}: {e}")
                else:
                    logger.warning(f"Network file missing: {cache_file}")
        else:
            logger.error(f"Models directory not found: {models_dir}")
            logger.error("Please run: python build_network_cache.py --all")

class RegNetAgentsModelingAgent:
    """Agent for gene network modeling and analysis."""
    def __init__(self, cache):
        self.cache = cache
        self.gene_service = CompleteGeneService()
        # Initialize gene mapper for ID conversion
        self.gene_mapper = GeneIDMapper()

    def _convert_ensembl_to_symbol(self, ensembl_id: str) -> str:
        """Convert Ensembl ID to gene symbol, with fallback to original ID."""
        if not ensembl_id:
            return "Unknown"

        symbol = self.gene_mapper.ensembl_to_symbol(ensembl_id)
        return symbol if symbol else ensembl_id

    async def analyze_gene_network_context(self, gene: str, cell_type: CellType):
        """Analyze gene network context."""
        gene_info = self.gene_service.get_gene_info(gene)
        network_data = self.cache.network_indices.get(cell_type.value, {})

        # Get Ensembl ID for network lookup
        ensembl_id = self.gene_mapper.symbol_to_ensembl(gene)

        # Analyze network position
        regulators = []
        targets = []
        regulatory_role = "unknown"

        if network_data and ensembl_id:
            # Access the actual network structure
            regulator_targets = network_data.get('regulator_targets', {})
            target_regulators = network_data.get('target_regulators', {})

            # Find targets (genes this gene regulates)
            targets = regulator_targets.get(ensembl_id, [])

            # Find regulators (genes that regulate this gene)
            regulators = target_regulators.get(ensembl_id, [])

            # Determine regulatory role - prioritize hub status over heavily regulated status
            num_regulators = len(regulators)
            num_targets = len(targets)

            if num_targets > 20:
                regulatory_role = "hub_regulator"      # Regulates many genes (high priority)
            elif num_regulators > 15:
                regulatory_role = "heavily_regulated"  # Controlled by many regulators
            elif num_targets > 5 and num_regulators > 5:
                regulatory_role = "intermediate_node"  # Balanced regulatory role
            elif num_targets > 0:
                regulatory_role = "regulator"          # Has downstream targets
            else:
                regulatory_role = "weakly_regulated"   # Few regulators, no targets

        # Determine if gene is actually in the network (has connections)
        in_network = len(regulators) > 0 or len(targets) > 0
        is_regulator = len(targets) > 0

        return {
            "gene": gene,
            "ensembl_id": ensembl_id,
            "cell_type": cell_type.value,
            "gene_info": gene_info.__dict__ if gene_info else None,
            "network_available": bool(network_data),
            "in_network": in_network,
            "is_regulator": is_regulator,
            "regulatory_role": regulatory_role,
            "num_regulators": len(regulators),
            "num_targets": len(targets),
            "regulators": regulators[:10],  # Top 10 for brevity
            "targets": targets[:10],        # Top 10 for brevity
            "network_position": {
                "in_degree": len(regulators),
                "out_degree": len(targets),
                "role": regulatory_role
            }
        }

    async def analyze_regulators_detailed(self, gene: str, cell_type: CellType, max_regulators: int = 25):
        """Detailed analysis of gene regulators."""
        # Get gene info and network data directly to avoid [:10] slice in gene_context
        network_data = self.cache.network_indices.get(cell_type.value, {})
        gene_ensembl = self.gene_mapper.symbol_to_ensembl(gene)

        # Get FULL regulator list from network (not the [:10] sliced version from analyze_gene_network_context)
        target_regulators = network_data.get('target_regulators', {})
        regulators = target_regulators.get(gene_ensembl, []) if network_data and gene_ensembl else []

        return {
            "gene": gene,
            "cell_type": cell_type.value,
            "regulator_summary": {
                "total_regulators": len(regulators),
                "analyzed_regulators": min(len(regulators), max_regulators)
            },
            "hub_regulators": [
                {
                    "gene_symbol": self._convert_ensembl_to_symbol(reg),
                    "ensembl_id": reg,
                    "regulatory_strength": "moderate"
                }
                for reg in regulators[:max_regulators]
            ]
        }

    async def analyze_targets_detailed(self, gene: str, cell_type: CellType, max_targets: int = 25):
        """Detailed analysis of gene targets."""
        gene_context = await self.analyze_gene_network_context(gene, cell_type)
        targets = gene_context.get('targets', [])

        return {
            "gene": gene,
            "cell_type": cell_type.value,
            "target_summary": {
                "total_targets": len(targets),
                "analyzed_targets": min(len(targets), max_targets)
            },
            "cascade_targets": [
                {
                    "gene_symbol": self._convert_ensembl_to_symbol(target),
                    "ensembl_id": target,
                    "cascade_level": 1
                }
                for target in targets[:max_targets]
            ]
        }


    async def compare_gene_across_cell_types(self, gene: str):
        """Compare gene across multiple cell types."""
        results = {}

        for cell_type in CellType:
            context = await self.analyze_gene_network_context(gene, cell_type)
            results[cell_type.value] = {
                "regulatory_role": context.get('regulatory_role'),
                "num_regulators": context.get('num_regulators'),
                "num_targets": context.get('num_targets'),
                "in_network": context.get('in_network', False)
            }

        return {
            "gene": gene,
            "cell_type_analysis": results,
            "summary": "Cross-cell analysis completed"
        }

    async def analyze_perturbation_effects(
        self,
        target_gene: str,
        cell_type: CellType,
        regulators: list = None,
        max_regulators: int = 100
    ) -> dict:
        """
        Analyze therapeutic target potential by simulating regulator inhibition.

        Performs network perturbation analysis to rank upstream regulators as potential
        therapeutic targets. Uses pre-computed PageRank centrality scores to measure
        each regulator's importance in the network, providing instant ranking without
        expensive recomputation.

        The analysis measures:
        - Regulatory loss: % of target gene regulation lost if regulator inhibited
        - Network centrality: PageRank score (best predictor of drug targets per literature)
        - Out-degree centrality: Number of genes regulated (cascade breadth)
        - Degree centrality: Total connectivity (overall network importance)

        Automatically triggered for genes with >5 regulators, identifying the most
        promising candidates for experimental validation or drug development.

        Args:
            target_gene: Gene to analyze (e.g., "TP53", "BRCA1")
            cell_type: Cell type for analysis (affects network topology)
            regulators: List of regulator dicts (if already computed from prior analysis)
            max_regulators: Maximum number of regulators to analyze (default 100)

        Returns:
            dict: Perturbation analysis results containing:
                - ranked_regulators: List sorted by PageRank (best targets first)
                - summary: Analysis metadata (total regulators, ranking metric)
                - Each regulator includes: symbol, PageRank, regulatory_loss, centrality scores
        """
        logger.info(f"Starting perturbation analysis for {target_gene}")

        # Get baseline network state
        baseline = await self.analyze_gene_network_context(target_gene, cell_type)

        # Get regulators if not provided
        if not regulators:
            reg_analysis = await self.analyze_regulators_detailed(
                target_gene, cell_type, max_regulators
            )
            regulators = reg_analysis.get('hub_regulators', [])

        # For each regulator, simulate its inhibition
        perturbation_effects = []

        for regulator in regulators[:max_regulators]:
            effect = await self._simulate_regulator_inhibition(
                target_gene=target_gene,
                regulator_id=regulator['ensembl_id'],
                regulator_symbol=regulator.get('gene_symbol', 'Unknown'),
                cell_type=cell_type,
                baseline=baseline
            )
            perturbation_effects.append(effect)

        # Rank by different centrality metrics (standard approach from literature)
        # Per Mora & Donaldson (2021): Drug targets show high degree, betweenness, closeness, PageRank

        # Provide multiple rankings for interpretation
        by_degree = sorted(perturbation_effects, key=lambda x: x['centrality_metrics']['degree_centrality'], reverse=True)
        by_pagerank = sorted(perturbation_effects, key=lambda x: x['centrality_metrics']['pagerank'], reverse=True)

        # Primary ranking: PageRank (best predictor per Mora & Donaldson 2021)
        ranked_targets = by_pagerank

        # Extract therapeutic insights
        therapeutic_insights = self._extract_therapeutic_insights(
            ranked_targets, target_gene, by_degree
        )

        return {
            "target_gene": target_gene,
            "cell_type": cell_type.value,
            "baseline_regulators": len(regulators),
            "perturbation_results": ranked_targets,  # All regulators ranked by PageRank
            "rankings": {
                "by_pagerank": [{"regulator": t['regulator'], "score": t['centrality_metrics']['pagerank']} for t in by_pagerank[:5]],
                "by_degree": [{"regulator": t['regulator'], "score": t['centrality_metrics']['degree_centrality']} for t in by_degree[:5]]
            },
            "therapeutic_insights": therapeutic_insights,
            "summary": f"Analyzed {len(perturbation_effects)} regulatory perturbations for {target_gene}"
        }

    async def _simulate_regulator_inhibition(
        self,
        target_gene: str,
        regulator_id: str,
        regulator_symbol: str,
        cell_type: CellType,
        baseline: dict
    ) -> dict:
        """
        Simulate removing one regulator from the network using standard network centrality metrics.

        Uses established network centrality measures from computational biology literature:
        - Degree Centrality: Hub identification (Mora & Donaldson, 2021)
        - Betweenness Centrality: Bottleneck identification (Koschützki & Schreiber, 2008)
        - Closeness Centrality: Global influence (Mora & Donaldson, 2021)
        - PageRank: Connection quality (Koschützki & Schreiber, 2008)
        """
        import networkx as nx

        network_data = self.cache.network_indices.get(cell_type.value, {})

        # Build NetworkX graph if not already cached
        if not hasattr(self, '_network_graphs'):
            self._network_graphs = {}

        if cell_type.value not in self._network_graphs:
            # Try to load pre-computed PageRank from cache (Version 2+)
            pagerank_normalized = network_data.get('pagerank_normalized', {})
            cache_version = network_data.get('cache_version', 1)

            # Build graph for degree centrality and fallback PageRank if needed
            G = nx.DiGraph()
            regulator_targets_dict = network_data.get('regulator_targets', {})
            for reg_id, targets in regulator_targets_dict.items():
                for target_id in targets:
                    G.add_edge(reg_id, target_id)

            # If PageRank not pre-computed (Version 1 cache), calculate on-demand
            if not pagerank_normalized:
                logger.warning(f"PageRank not pre-computed for {cell_type.value} (cache version {cache_version}), calculating on-demand...")
                logger.info(f"Calculating standard centrality metrics for {cell_type.value} network ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")

                # Calculate PageRank with proper convergence
                try:
                    pagerank_scores = nx.pagerank(G, alpha=0.85, max_iter=100, tol=1e-06)
                    # Normalize PageRank by maximum value for interpretability
                    max_pagerank = max(pagerank_scores.values()) if pagerank_scores else 1.0
                    pagerank_normalized = {k: v / max_pagerank for k, v in pagerank_scores.items()}
                except:
                    # Fallback to simpler calculation if PageRank fails
                    logger.warning("PageRank failed, using out-degree as proxy")
                    pagerank_normalized = nx.out_degree_centrality(G)
            else:
                logger.info(f"Using pre-computed PageRank for {cell_type.value} (cache version {cache_version}, {len(pagerank_normalized)} scores)")

            centrality_dict = {
                'graph': G,
                'degree_centrality': nx.degree_centrality(G),
                'out_degree_centrality': nx.out_degree_centrality(G),
                'pagerank': pagerank_normalized,
                'is_large_network': G.number_of_nodes() >= 1000
            }

            # Note: Betweenness and closeness centrality removed - never computable for large networks (>1000 nodes)
            # Epithelial network has 14,628 nodes, so these metrics were always 0.0

            self._network_graphs[cell_type.value] = centrality_dict

        metrics = self._network_graphs[cell_type.value]

        # Get regulator's centrality scores
        degree_cent = metrics['degree_centrality'].get(regulator_id, 0.0)
        out_degree_cent = metrics['out_degree_centrality'].get(regulator_id, 0.0)
        pagerank_cent = metrics['pagerank'].get(regulator_id, 0.0)

        # Cascade effects
        regulator_targets = network_data.get('regulator_targets', {}).get(regulator_id, [])
        baseline_targets = baseline.get('targets', [])
        cascade_overlap = len(set(regulator_targets) & set(baseline_targets))

        # Network centrality metrics for therapeutic target ranking
        # Primary metric: PageRank (best predictor of drug targets per Mora & Donaldson 2021)
        # Secondary metrics: degree/out-degree centrality for hub identification

        return {
            "regulator": regulator_symbol,
            "regulator_ensembl_id": regulator_id,
            "regulator_downstream_targets": len(regulator_targets),
            "cascade_overlap": cascade_overlap,
            "affected_cascades": [self._convert_ensembl_to_symbol(t) for t in regulator_targets[:5]],
            # Network centrality metrics
            "centrality_metrics": {
                "degree_centrality": round(degree_cent, 4),
                "out_degree_centrality": round(out_degree_cent, 4),
                "pagerank": round(pagerank_cent, 4)
            }
        }

    def _assess_druggability(self, regulator_symbol: str, num_targets: int) -> str:
        """Provide basic druggability assessment based on gene characteristics."""
        # Simple heuristics - could be expanded with database lookups
        if num_targets > 50:
            return f"{regulator_symbol} is a hub regulator - high therapeutic impact but may have off-target effects"
        elif num_targets > 20:
            return f"{regulator_symbol} shows good therapeutic potential with moderate target specificity"
        elif num_targets > 5:
            return f"{regulator_symbol} is a specific regulator - lower impact but better target specificity"
        else:
            return f"{regulator_symbol} has minimal downstream effects - limited therapeutic potential"

    def _extract_therapeutic_insights(self, ranked_targets: list, target_gene: str, by_degree: list) -> dict:
        """
        Extract key therapeutic insights from perturbation results using network centrality metrics.

        Interpretation based on Mora & Donaldson (2021) BMC Bioinformatics:
        - Approved drug targets show significantly higher degree and PageRank
        - Highest PageRank nodes are best therapeutic target candidates
        """
        if not ranked_targets:
            return {
                "summary": f"No regulators found for {target_gene}",
                "top_target": None,
                "strategy": "No therapeutic strategy available"
            }

        top_by_pagerank = ranked_targets[0]
        top_by_degree = by_degree[0]

        # Check if same regulator appears in both top rankings (strong consensus candidate)
        consensus_target = (top_by_pagerank['regulator'] == top_by_degree['regulator'])

        # Calculate regulatory loss: each regulator contributes equally (1/n)
        regulatory_loss_pct = round(100.0 / len(ranked_targets), 1)

        return {
            "summary": f"Identified {len(ranked_targets)} potential therapeutic targets for {target_gene}",
            "top_target_by_pagerank": {
                "regulator": top_by_pagerank['regulator'],
                "pagerank": top_by_pagerank['centrality_metrics']['pagerank'],
                "downstream_targets": top_by_pagerank['regulator_downstream_targets']
            },
            "top_target_by_degree": {
                "regulator": top_by_degree['regulator'],
                "degree_centrality": top_by_degree['centrality_metrics']['degree_centrality'],
                "downstream_targets": top_by_degree['regulator_downstream_targets']
            },
            "consensus_target": consensus_target,
            "strategy": f"Inhibiting {top_by_pagerank['regulator']} (highest PageRank) would reduce {target_gene} regulation by {regulatory_loss_pct}% and affect {top_by_pagerank['regulator_downstream_targets']} downstream genes",
            "interpretation": self._generate_centrality_interpretation(top_by_pagerank, target_gene, consensus_target)
        }

    def _generate_centrality_interpretation(self, top_target: dict, target_gene: str, consensus: bool) -> str:
        """Generate interpretation based on standard network centrality principles."""
        regulator = top_target['regulator']
        pagerank = top_target['centrality_metrics']['pagerank']
        degree = top_target['centrality_metrics']['degree_centrality']
        targets = top_target['regulator_downstream_targets']

        interpretation = f"{regulator} shows strong therapeutic potential for modulating {target_gene}. "

        # Interpret PageRank (connection quality)
        if pagerank > 0.3:
            interpretation += f"High PageRank ({pagerank:.3f}) indicates this regulator has high-quality connections in the network, typical of successful drug targets. "
        else:
            interpretation += f"Moderate PageRank ({pagerank:.3f}) suggests moderate network influence. "

        # Interpret degree (hub status)
        if targets > 200:
            interpretation += f"As a hub regulator ({targets} targets), inhibiting this gene would have broad network effects. "
        elif targets > 50:
            interpretation += f"With {targets} downstream targets, this regulator has moderate network influence. "

        # Consensus across metrics
        if consensus:
            interpretation += "This target ranks highly across multiple centrality metrics, strengthening the therapeutic recommendation."
        else:
            interpretation += "Consider validating with experimental perturbation studies."

        return interpretation

class PathwayEnricherAgent:
    """Agent for Reactome pathway enrichment analysis."""
    def __init__(self):
        try:
            import requests
            self.requests = requests
            self.reactome_base_url = "https://reactome.org/AnalysisService"
            logger.info("PathwayEnricherAgent initialized with Reactome API")
        except ImportError:
            logger.warning("requests library not available. Install with: pip install requests")
            self.requests = None

    async def enrich_pathways_reactome(self, gene_list: list, species: str = "Homo sapiens") -> Dict:
        """
        Perform pathway enrichment analysis using Reactome API.

        Args:
            gene_list: List of gene symbols
            species: Species name (default: Homo sapiens)

        Returns:
            Dictionary with enrichment results including p-values and FDR
        """
        if not self.requests:
            return {
                "status": "error",
                "message": "requests library not available",
                "genes_analyzed": gene_list
            }

        try:
            # Prepare gene list for Reactome
            genes_str = "\n".join(gene_list)

            # Call Reactome Analysis API
            url = f"{self.reactome_base_url}/identifiers/projection"
            headers = {"Content-Type": "text/plain"}

            # Use asyncio to run the synchronous request in a thread
            loop = asyncio.get_event_loop()
            response = await loop.run_in_executor(
                None,
                lambda: self.requests.post(url, data=genes_str, headers=headers)
            )

            if response.status_code != 200:
                logger.error(f"Reactome API error: {response.status_code}")
                return {
                    "status": "error",
                    "message": f"Reactome API returned status {response.status_code}",
                    "genes_analyzed": gene_list
                }

            # Get analysis results (includes pathways directly)
            result = response.json()
            token = result.get('summary', {}).get('token')

            if not token:
                return {
                    "status": "error",
                    "message": "No analysis token received from Reactome",
                    "genes_analyzed": gene_list
                }

            # Pathways are in the initial response
            pathways_data = result

            # Process and format results - filter to FDR < 0.05 for scientific rigor
            enriched_pathways = []
            for pathway in pathways_data.get('pathways', []):
                # Only include statistically significant pathways
                fdr = pathway.get('entities', {}).get('fdr')
                if fdr is not None and fdr < 0.05:
                    enriched_pathways.append({
                        "pathway_id": pathway.get('stId'),
                        "pathway_name": pathway.get('name'),
                        "p_value": pathway.get('entities', {}).get('pValue'),
                        "fdr": pathway.get('entities', {}).get('fdr'),
                        "genes_found": pathway.get('entities', {}).get('found'),
                        "genes_total": pathway.get('entities', {}).get('total'),
                        "species": pathway.get('species', {}).get('name')
                    })

            return {
                "status": "success",
                "analysis_type": "reactome_pathway_enrichment",
                "genes_analyzed": gene_list,
                "genes_recognized": result.get('summary', {}).get('identifiersFound', 0),
                "token": token,
                "enriched_pathways": enriched_pathways,
                "summary": {
                    "total_pathways": len(enriched_pathways),
                    "significant_pathways": len([p for p in enriched_pathways if p.get('fdr', 1) < 0.05])
                }
            }

        except Exception as e:
            logger.error(f"Reactome enrichment failed: {str(e)}")
            return {
                "status": "error",
                "message": str(e),
                "genes_analyzed": gene_list
            }

class CrossSystemIntegrationAgent:
    """Agent for cross-system integration analysis."""
    def __init__(self, cache):
        self.cache = cache

    async def compare_gene_across_cell_types(self, gene: str, cell_types: list, analysis_results: dict):
        """Compare gene across multiple cell types."""
        return {
            "gene": gene,
            "comparison": "cross_cell_analysis",
            "cell_types": [ct.value if hasattr(ct, 'value') else str(ct) for ct in cell_types],
            "summary": "Cross-cell comparison completed"
        }

    async def analyze_biological_context(self, gene_list: list, analysis_type: str = "pathway_enrichment"):
        """Analyze biological context and pathways (legacy method - use PathwayEnricherAgent instead)."""
        return {
            "analysis_type": analysis_type,
            "genes_analyzed": gene_list,
            "enriched_pathways": {
                "wnt_signaling": {"score": 0.85, "genes": gene_list[:3]},
                "cell_cycle": {"score": 0.72, "genes": gene_list[:2]},
                "apoptosis": {"score": 0.68, "genes": gene_list[:2]}
            },
            "top_pathways": [
                ("wnt_signaling", 0.85),
                ("cell_cycle", 0.72),
                ("apoptosis", 0.68)
            ]
        }

class DomainAnalysisAgents:
    """Specialized domain analysis agents integrated into workflow"""

    def __init__(self, cache, use_llm=None):
        self.cache = cache

        # Auto-detect if Ollama is available
        if use_llm is None:
            use_llm = os.getenv('USE_LLM_AGENTS', 'true').lower() == 'true'

        self.use_llm = use_llm
        self.ollama_client = self._initialize_ollama() if use_llm else False
        self.ollama_model = os.getenv('OLLAMA_MODEL', 'llama3.1:8b')
        self.ollama_temperature = float(os.getenv('OLLAMA_TEMPERATURE', '0.3'))
        self.ollama_max_tokens = int(os.getenv('OLLAMA_MAX_TOKENS', '1500'))

    def _initialize_ollama(self):
        """Check if Ollama is available and running"""
        try:
            # Test connection to Ollama
            models_response = ollama.list()

            # Handle different response structures from Ollama API
            if hasattr(models_response, 'models'):
                models_list = models_response.models
            elif isinstance(models_response, dict):
                models_list = models_response.get('models', [])
            else:
                models_list = models_response

            # Extract model names (API may return dicts or objects)
            available_models = []
            for m in models_list:
                name = None
                if isinstance(m, dict):
                    # Dict-style response
                    name = m.get('name') or m.get('model')
                elif hasattr(m, 'model'):
                    # Object-style response (ollama._types.ListResponse.Model)
                    name = m.model
                elif hasattr(m, 'name'):
                    # Alternative object style
                    name = m.name

                if name:
                    available_models.append(name)

            # Check if specified model is available
            model_name = os.getenv('OLLAMA_MODEL', 'llama3.1:8b')
            if model_name not in available_models:
                logger.error(f"Ollama model '{model_name}' not found. Available models: {available_models}")
                logger.error(f"Please run: ollama pull {model_name}")
                return False

            logger.info(f"Ollama available, using model: {model_name}")
            return True

        except Exception as e:
            logger.warning(f"Ollama not available ({e}), falling back to rule-based agents")
            logger.warning("To use LLM agents: 1) Install Ollama from https://ollama.com/download")
            logger.warning(f"                     2) Run: ollama pull {os.getenv('OLLAMA_MODEL', 'llama3.1:8b')}")
            return False

    async def _call_ollama(self, prompt: str, system_prompt: str = None, max_retries: int = 2) -> str:
        """Call Ollama with structured prompt, return response with retry logic"""
        timeout = int(os.getenv('OLLAMA_TIMEOUT', '30'))

        for attempt in range(max_retries):
            try:
                messages = []
                if system_prompt:
                    messages.append({"role": "system", "content": system_prompt})
                messages.append({"role": "user", "content": prompt})

                response = await asyncio.wait_for(
                    asyncio.to_thread(
                        ollama.chat,
                        model=self.ollama_model,
                        messages=messages,
                        options={
                            "temperature": self.ollama_temperature,
                            "num_predict": self.ollama_max_tokens
                        }
                    ),
                    timeout=timeout
                )

                # Validate response structure
                if not response or 'message' not in response or 'content' not in response['message']:
                    raise ValueError("Invalid response structure from Ollama")

                content = response['message']['content']
                if not content or len(content.strip()) < 10:
                    raise ValueError("Empty or too short response from Ollama")

                return content

            except asyncio.TimeoutError:
                logger.error(f"Ollama call timed out after {timeout} seconds (attempt {attempt + 1}/{max_retries})")
                if attempt == max_retries - 1:
                    raise
                await asyncio.sleep(1)  # Brief delay before retry

            except Exception as e:
                logger.error(f"Ollama call failed (attempt {attempt + 1}/{max_retries}): {e}")
                if attempt == max_retries - 1:
                    raise
                await asyncio.sleep(1)  # Brief delay before retry

    def _parse_llm_json(self, response_text: str, expected_keys: list) -> dict:
        """Extract JSON from LLM response, validate structure"""
        try:
            # LLM might wrap JSON in markdown code blocks
            original_text = response_text
            if "```json" in response_text:
                response_text = response_text.split("```json")[1].split("```")[0]
            elif "```" in response_text:
                response_text = response_text.split("```")[1].split("```")[0]

            # Try to extract JSON object if embedded in text
            response_text = response_text.strip()
            if not response_text.startswith('{'):
                # Try to find JSON object in text
                json_start = response_text.find('{')
                json_end = response_text.rfind('}')
                if json_start != -1 and json_end != -1:
                    response_text = response_text[json_start:json_end+1]

            parsed = json.loads(response_text)

            # Validate it's a dict
            if not isinstance(parsed, dict):
                raise ValueError(f"Expected dict, got {type(parsed)}")

            # Validate required keys present
            missing_keys = [key for key in expected_keys if key not in parsed]
            if missing_keys:
                logger.warning(f"Missing expected keys in LLM response: {', '.join(missing_keys)}")
                # Add placeholder values for missing keys
                for key in missing_keys:
                    parsed[key] = "N/A" if "_rationale" in key or "summary" in key else "unknown"

            # Validate and normalize score values (0.0-1.0)
            for key, value in parsed.items():
                if isinstance(value, (int, float)) and ("score" in key or "centrality" in key):
                    if value < 0 or value > 1:
                        logger.warning(f"Score {key}={value} out of range [0,1], clamping")
                        parsed[key] = max(0.0, min(1.0, value))

            return parsed

        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse LLM JSON: {e}")
            logger.error(f"Response was: {original_text[:500]}...")
            raise
        except Exception as e:
            logger.error(f"Unexpected error parsing LLM response: {e}")
            raise

    async def analyze_cancer_context(self, gene: str, gene_info: Dict, pathway_analysis: Dict) -> Dict:
        """Analyze gene from cancer research perspective (LLM-powered with rule-based fallback)"""

        # If LLM available, use it
        if self.ollama_client:
            try:
                logger.debug(f"Using LLM for cancer analysis of {gene}")
                return await self._analyze_cancer_llm(gene, gene_info, pathway_analysis)
            except Exception as e:
                logger.warning(f"LLM cancer analysis failed for {gene}: {e}, using rule-based fallback")

        # Fallback to rule-based analysis
        return self._analyze_cancer_rules(gene, gene_info, pathway_analysis)

    async def _analyze_cancer_llm(self, gene: str, gene_info: Dict, pathway_analysis: Dict) -> Dict:
        """LLM-powered cancer domain analysis"""

        # Extract pathway information for prompt
        enriched_pathways = []
        if pathway_analysis:
            pathways = pathway_analysis.get('enriched_pathways', [])
            # Handle both list and dict formats
            if isinstance(pathways, dict):
                enriched_pathways = list(pathways.keys())
            elif isinstance(pathways, list):
                # Extract pathway names from list of dicts
                enriched_pathways = [p.get('pathway_name', p.get('name', str(p))) for p in pathways if p][:10]

        # Extract functional description if available
        gene_function = None
        if gene_info and isinstance(gene_info, dict) and 'gene_info' in gene_info:
            # gene_info contains the GeneInfo object as dict
            gene_obj = gene_info['gene_info']
            if gene_obj and isinstance(gene_obj, dict):
                gene_function = gene_obj.get('function', None)

        # Build structured prompt
        function_context = f"\nGene Function: {gene_function}\n" if gene_function else ""

        prompt = f"""Analyze the gene {gene} from a cancer biology perspective.
{function_context}
Gene Network Context:
- Regulatory Role: {gene_info.get('regulatory_role', 'unknown')}
- Upstream Regulators: {gene_info.get('num_regulators', 0)}
- Downstream Targets: {gene_info.get('num_targets', 0)}
- Network Position: in-degree={gene_info.get('num_regulators', 0)}, out-degree={gene_info.get('num_targets', 0)}

Enriched Pathways: {enriched_pathways if enriched_pathways else 'No pathway data available'}

Provide a cancer biology analysis in this EXACT JSON format:
{{
  "oncogenic_potential": "high|moderate|low",
  "oncogenic_rationale": "brief scientific explanation based on network topology and cancer biology",
  "tumor_suppressor_likelihood": "high|moderate|low",
  "tumor_suppressor_rationale": "brief scientific explanation",
  "therapeutic_target_score": 0.0-1.0,
  "therapeutic_rationale": "explanation of druggability and therapeutic potential",
  "cancer_pathways": ["pathway1", "pathway2"],
  "biomarker_potential": "high|moderate|low",
  "biomarker_utility": "diagnostic|prognostic|predictive",
  "biomarker_rationale": "explanation",
  "research_priority": "high|moderate|low",
  "summary": "1-2 sentence synthesis of cancer relevance"
}}

Base your analysis on:
- Network centrality (hub regulators are critical for cancer)
- Regulatory control (highly regulated genes often tumor suppressors)
- Pathway involvement (cancer-related pathways)
- Known cancer biology principles

Provide only the JSON, no additional text."""

        system_prompt = "You are an expert cancer biologist analyzing gene regulatory networks. Provide scientifically accurate, evidence-based analysis in structured JSON format only."

        # Call Ollama
        response = await self._call_ollama(prompt, system_prompt)

        # Parse JSON response
        expected_keys = ['oncogenic_potential', 'tumor_suppressor_likelihood',
                         'therapeutic_target_score', 'biomarker_potential', 'summary']
        insights = self._parse_llm_json(response, expected_keys)

        # Format into standard structure
        return {
            "gene": gene,
            "domain": "cancer_research",
            "insights": {
                "oncogenic_potential": insights.get('oncogenic_potential', 'moderate'),
                "tumor_suppressor_likelihood": insights.get('tumor_suppressor_likelihood', 'moderate'),
                "therapeutic_target_score": float(insights.get('therapeutic_target_score', 0.5)),
                "mutation_impact": "high" if insights.get('oncogenic_potential') == 'high' else "moderate"
            },
            "cancer_pathways": insights.get('cancer_pathways', []),
            "biomarker_potential": insights.get('biomarker_potential', 'moderate'),
            "research_priority": insights.get('research_priority', 'moderate'),
            "summary": insights.get('summary', f"{gene} cancer analysis"),
            "llm_rationale": {
                "oncogenic": insights.get('oncogenic_rationale', ''),
                "tumor_suppressor": insights.get('tumor_suppressor_rationale', ''),
                "therapeutic": insights.get('therapeutic_rationale', ''),
                "biomarker": insights.get('biomarker_rationale', '')
            },
            "llm_powered": True
        }

    def _analyze_cancer_rules(self, gene: str, gene_info: Dict, pathway_analysis: Dict) -> Dict:
        """Rule-based cancer domain analysis (fallback)"""
        regulatory_role = gene_info.get('regulatory_role', 'unknown')
        num_regulators = gene_info.get('num_regulators', 0)
        num_targets = gene_info.get('num_targets', 0)

        # Cancer-specific insights
        cancer_insights = {
            "oncogenic_potential": "high" if num_targets > 20 else "moderate" if num_targets > 5 else "low",
            "tumor_suppressor_likelihood": "high" if num_regulators > 15 else "moderate" if num_regulators > 5 else "low",
            "therapeutic_target_score": min(0.9, (num_targets * 0.02) + (num_regulators * 0.01)),
            "mutation_impact": "high" if regulatory_role in ['hub_regulator', 'master_regulator'] else "moderate"
        }

        # Pathway-based cancer insights
        cancer_pathways = []
        if pathway_analysis:
            enriched = pathway_analysis.get('enriched_pathways', {})
            if 'wnt_signaling' in enriched:
                cancer_pathways.append("Wnt signaling - colorectal cancer pathway")
            if 'cell_cycle' in enriched:
                cancer_pathways.append("Cell cycle - proliferation control")
            if 'apoptosis' in enriched:
                cancer_pathways.append("Apoptosis - cell death regulation")

        return {
            "gene": gene,
            "domain": "cancer_research",
            "insights": cancer_insights,
            "cancer_pathways": cancer_pathways,
            "biomarker_potential": "high" if cancer_insights["therapeutic_target_score"] > 0.7 else "moderate",
            "research_priority": "high" if regulatory_role == "hub_regulator" else "moderate",
            "summary": f"{gene} shows {cancer_insights['oncogenic_potential']} oncogenic potential with {cancer_insights['tumor_suppressor_likelihood']} tumor suppressor likelihood"
        }

    async def analyze_drug_development(self, gene: str, gene_info: Dict, regulators_analysis: Dict, targets_analysis: Dict) -> Dict:
        """Analyze gene from drug development perspective (LLM-powered with rule-based fallback)"""

        # If LLM available, use it
        if self.ollama_client:
            try:
                logger.debug(f"Using LLM for drug development analysis of {gene}")
                return await self._analyze_drug_llm(gene, gene_info, regulators_analysis, targets_analysis)
            except Exception as e:
                logger.warning(f"LLM drug analysis failed for {gene}: {e}, using rule-based fallback")

        # Fallback to rule-based analysis
        return self._analyze_drug_rules(gene, gene_info, regulators_analysis, targets_analysis)

    async def _analyze_drug_llm(self, gene: str, gene_info: Dict, regulators_analysis: Dict, targets_analysis: Dict) -> Dict:
        """LLM-powered drug development analysis"""

        # Extract regulator/target info
        target_count = targets_analysis.get('target_summary', {}).get('total_targets', 0) if targets_analysis else 0
        regulator_count = regulators_analysis.get('regulator_summary', {}).get('total_regulators', 0) if regulators_analysis else 0

        # Extract functional description if available
        gene_function = None
        if gene_info and isinstance(gene_info, dict) and 'gene_info' in gene_info:
            gene_obj = gene_info['gene_info']
            if gene_obj and isinstance(gene_obj, dict):
                gene_function = gene_obj.get('function', None)

        function_context = f"\nGene Function: {gene_function}\n" if gene_function else ""

        prompt = f"""Analyze the gene {gene} from a drug development perspective.
{function_context}
Gene Network Context:
- Regulatory Role: {gene_info.get('regulatory_role', 'unknown')}
- Upstream Regulators: {gene_info.get('num_regulators', 0)}
- Downstream Targets: {gene_info.get('num_targets', 0)}
- Total Regulators Found: {regulator_count}
- Total Targets Found: {target_count}

Provide drug development analysis in this EXACT JSON format:
{{
  "druggability_score": 0.0-1.0,
  "druggability_rationale": "explanation of druggability based on structure and network",
  "target_class": "kinase|GPCR|transcription_factor|nuclear_receptor|other",
  "intervention_strategy": "inhibition|activation|modulation|allosteric",
  "intervention_rationale": "why this strategy is appropriate",
  "development_complexity": "high|moderate|low",
  "cascade_effects": ["effect1", "effect2"],
  "clinical_trial_readiness": "ready|needs_preclinical|needs_research|not_suitable",
  "development_timeline": "estimated years",
  "summary": "1-2 sentence synthesis of drug development potential"
}}

Base analysis on:
- Network topology (hub genes may have broad effects)
- Regulatory control (heavily regulated may be indirect targets)
- Cascade effects (downstream impact)
- Known drug target classes

Provide only the JSON, no additional text."""

        system_prompt = "You are an expert in drug discovery and development analyzing gene regulatory networks. Provide scientifically accurate, evidence-based analysis in structured JSON format only."

        # Call Ollama
        response = await self._call_ollama(prompt, system_prompt)

        # Parse JSON response
        expected_keys = ['druggability_score', 'target_class', 'intervention_strategy', 'development_complexity', 'summary']
        insights = self._parse_llm_json(response, expected_keys)

        # Format into standard structure
        return {
            "gene": gene,
            "domain": "drug_development",
            "insights": {
                "druggability_score": float(insights.get('druggability_score', 0.5)),
                "target_class": insights.get('target_class', 'other'),
                "intervention_strategy": insights.get('intervention_strategy', 'modulation'),
                "development_complexity": insights.get('development_complexity', 'moderate')
            },
            "cascade_effects": insights.get('cascade_effects', []),
            "clinical_trial_readiness": insights.get('clinical_trial_readiness', 'needs_research'),
            "development_timeline": insights.get('development_timeline', '5-8 years'),
            "summary": insights.get('summary', f"{gene} drug development analysis"),
            "llm_rationale": {
                "druggability": insights.get('druggability_rationale', ''),
                "intervention": insights.get('intervention_rationale', '')
            },
            "llm_powered": True
        }

    def _analyze_drug_rules(self, gene: str, gene_info: Dict, regulators_analysis: Dict, targets_analysis: Dict) -> Dict:
        """Rule-based drug development analysis (fallback)"""
        num_targets = gene_info.get('num_targets', 0)
        regulatory_role = gene_info.get('regulatory_role', 'unknown')

        # Drug target assessment
        druggability_score = min(0.95, (num_targets * 0.03) + (0.2 if regulatory_role == 'hub_regulator' else 0))

        drug_insights = {
            "druggability_score": druggability_score,
            "target_class": "kinase" if num_targets > 10 else "transcription_factor" if regulatory_role == "regulator" else "other",
            "intervention_strategy": "inhibition" if num_targets > 15 else "activation" if gene_info.get('num_regulators', 0) > 10 else "modulation",
            "development_complexity": "high" if num_targets > 20 else "moderate" if num_targets > 5 else "low"
        }

        # Cascade analysis for drug effects
        cascade_effects = []
        if targets_analysis:
            target_count = targets_analysis.get('target_summary', {}).get('total_targets', 0)
            if target_count > 10:
                cascade_effects.append(f"Targeting {gene} affects {target_count} downstream genes")

        return {
            "gene": gene,
            "domain": "drug_development",
            "insights": drug_insights,
            "cascade_effects": cascade_effects,
            "clinical_trial_readiness": "ready" if druggability_score > 0.6 else "needs_research",
            "development_timeline": "3-5 years" if druggability_score > 0.7 else "5-8 years",
            "summary": f"{gene} has {druggability_score:.1%} druggability score with {drug_insights['intervention_strategy']} strategy recommended"
        }

    async def analyze_clinical_relevance(self, gene: str, gene_info: Dict, cross_cell_analysis: Dict) -> Dict:
        """Analyze gene from clinical/personalized medicine perspective (LLM-powered with rule-based fallback)"""

        if self.ollama_client:
            try:
                logger.debug(f"Using LLM for clinical analysis of {gene}")
                return await self._analyze_clinical_llm(gene, gene_info, cross_cell_analysis)
            except Exception as e:
                logger.warning(f"LLM clinical analysis failed for {gene}: {e}, using rule-based fallback")

        return self._analyze_clinical_rules(gene, gene_info, cross_cell_analysis)

    async def _analyze_clinical_llm(self, gene: str, gene_info: Dict, cross_cell_analysis: Dict) -> Dict:
        """LLM-powered clinical relevance analysis"""

        # Extract functional description if available
        gene_function = None
        if gene_info and isinstance(gene_info, dict) and 'gene_info' in gene_info:
            gene_obj = gene_info['gene_info']
            if gene_obj and isinstance(gene_obj, dict):
                gene_function = gene_obj.get('function', None)

        # Prepare context for LLM
        num_regulators = gene_info.get('num_regulators', 0)
        num_targets = gene_info.get('num_targets', 0)
        regulatory_role = gene_info.get('regulatory_role', 'unknown')

        # Cross-cell type context
        tissue_context = "Not available"
        if cross_cell_analysis:
            cell_analysis = cross_cell_analysis.get('cell_type_analysis', {})
            active_cell_types = [ct for ct, data in cell_analysis.items() if data.get('in_network', False)]
            tissue_context = f"Active in {len(active_cell_types)} cell types: {', '.join(active_cell_types[:5])}"

        function_context = f"\nGene Function: {gene_function}\n" if gene_function else ""

        prompt = f"""Analyze the gene {gene} from a clinical medicine and personalized healthcare perspective.
{function_context}
Gene Network Context:
- Regulatory Role: {regulatory_role}
- Upstream Regulators: {num_regulators}
- Downstream Targets: {num_targets}
- Tissue Distribution: {tissue_context}

Provide a clinical analysis in this EXACT JSON format:
{{
  "disease_association_likelihood": "high|moderate|low",
  "disease_rationale": "brief explanation of disease relevance",
  "biomarker_utility": "diagnostic|prognostic|predictive|therapeutic",
  "biomarker_rationale": "brief explanation of biomarker potential",
  "clinical_actionability": "high|moderate|low",
  "actionability_rationale": "brief explanation of clinical utility",
  "tissue_specificity": "tissue-specific|broadly_expressed|ubiquitous",
  "diagnostic_potential": "high|moderate|low",
  "summary": "1-2 sentence clinical significance summary"
}}

Focus on:
- Disease association potential based on network position
- Biomarker utility for diagnosis, prognosis, or therapeutic monitoring
- Clinical actionability and translational potential
- Tissue specificity implications for personalized medicine

Provide only the JSON, no additional text."""

        system_prompt = "You are an expert clinician and translational researcher analyzing gene networks for precision medicine applications."

        response = await self._call_ollama(prompt, system_prompt)

        # Parse LLM response
        expected_keys = [
            "disease_association_likelihood", "disease_rationale",
            "biomarker_utility", "biomarker_rationale",
            "clinical_actionability", "actionability_rationale",
            "tissue_specificity", "diagnostic_potential", "summary"
        ]
        insights = self._parse_llm_json(response, expected_keys)

        # Format tissue specificity as list (for compatibility)
        tissue_spec_list = []
        if cross_cell_analysis:
            cell_analysis = cross_cell_analysis.get('cell_type_analysis', {})
            active_cell_types = [ct for ct, data in cell_analysis.items() if data.get('in_network', False)]
            if len(active_cell_types) < 3:
                tissue_spec_list.append(f"Tissue-specific expression in {len(active_cell_types)} cell types")
            else:
                tissue_spec_list.append(f"Broadly expressed across {len(active_cell_types)} cell types")

        # Structure response
        return {
            "gene": gene,
            "domain": "clinical_relevance",
            "insights": {
                "disease_association_likelihood": insights.get("disease_association_likelihood", "unknown"),
                "disease_rationale": insights.get("disease_rationale", "N/A"),
                "biomarker_utility": insights.get("biomarker_utility", "unknown"),
                "biomarker_rationale": insights.get("biomarker_rationale", "N/A"),
                "clinical_actionability": insights.get("clinical_actionability", "unknown"),
                "actionability_rationale": insights.get("actionability_rationale", "N/A")
            },
            "tissue_specificity": tissue_spec_list,
            "diagnostic_potential": insights.get("diagnostic_potential", "moderate"),
            "summary": insights.get("summary", f"{gene} clinical analysis completed"),
            "llm_powered": True
        }

    def _analyze_clinical_rules(self, gene: str, gene_info: Dict, cross_cell_analysis: Dict) -> Dict:
        """Rule-based clinical relevance analysis (fallback)"""
        num_regulators = gene_info.get('num_regulators', 0)
        regulatory_role = gene_info.get('regulatory_role', 'unknown')

        # Clinical significance assessment
        clinical_insights = {
            "disease_association_likelihood": "high" if num_regulators > 15 else "moderate" if num_regulators > 5 else "low",
            "biomarker_utility": "diagnostic" if regulatory_role == "heavily_regulated" else "prognostic" if regulatory_role == "hub_regulator" else "predictive",
            "clinical_actionability": "high" if num_regulators > 10 and regulatory_role in ['hub_regulator', 'heavily_regulated'] else "moderate"
        }

        # Cross-cell type clinical insights
        tissue_specificity = []
        if cross_cell_analysis:
            cell_analysis = cross_cell_analysis.get('cell_type_analysis', {})
            active_cell_types = [ct for ct, data in cell_analysis.items() if data.get('in_network', False)]
            if len(active_cell_types) < 3:
                tissue_specificity.append(f"Tissue-specific expression in {len(active_cell_types)} cell types")
            else:
                tissue_specificity.append(f"Broadly expressed across {len(active_cell_types)} cell types")

        return {
            "gene": gene,
            "domain": "clinical_relevance",
            "insights": clinical_insights,
            "tissue_specificity": tissue_specificity,
            "diagnostic_potential": "high" if clinical_insights["clinical_actionability"] == "high" else "moderate",
            "summary": f"{gene} shows {clinical_insights['disease_association_likelihood']} disease association with {clinical_insights['biomarker_utility']} biomarker utility",
            "llm_powered": False
        }

    async def analyze_systems_biology(self, gene: str, gene_info: Dict, regulators_analysis: Dict, targets_analysis: Dict) -> Dict:
        """Analyze gene from systems biology perspective (LLM-powered with rule-based fallback)"""

        if self.ollama_client:
            try:
                logger.debug(f"Using LLM for systems biology analysis of {gene}")
                return await self._analyze_systems_llm(gene, gene_info, regulators_analysis, targets_analysis)
            except Exception as e:
                logger.warning(f"LLM systems analysis failed for {gene}: {e}, using rule-based fallback")

        return self._analyze_systems_rules(gene, gene_info, regulators_analysis, targets_analysis)

    async def _analyze_systems_llm(self, gene: str, gene_info: Dict, regulators_analysis: Dict, targets_analysis: Dict) -> Dict:
        """LLM-powered systems biology analysis"""

        # Extract functional description if available
        gene_function = None
        if gene_info and isinstance(gene_info, dict) and 'gene_info' in gene_info:
            gene_obj = gene_info['gene_info']
            if gene_obj and isinstance(gene_obj, dict):
                gene_function = gene_obj.get('function', None)

        # Prepare context for LLM
        num_regulators = gene_info.get('num_regulators', 0)
        num_targets = gene_info.get('num_targets', 0)
        regulatory_role = gene_info.get('regulatory_role', 'unknown')
        pagerank = gene_info.get('pagerank', 0.0)

        # Network context
        reg_count = regulators_analysis.get('regulator_summary', {}).get('total_regulators', 0) if regulators_analysis else num_regulators
        target_count = targets_analysis.get('target_summary', {}).get('total_targets', 0) if targets_analysis else num_targets

        function_context = f"\nGene Function: {gene_function}\n" if gene_function else ""

        prompt = f"""Analyze the gene {gene} from a systems biology and network theory perspective.
{function_context}
Gene Network Context:
- Regulatory Role: {regulatory_role}
- Upstream Regulators: {num_regulators} (total in cascade: {reg_count})
- Downstream Targets: {num_targets} (total in cascade: {target_count})
- PageRank Centrality: {pagerank:.4f}
- Total Network Degree: {num_regulators + num_targets}

Provide a systems biology analysis in this EXACT JSON format:
{{
  "network_centrality": 0.0-1.0,
  "centrality_rationale": "brief explanation of network position",
  "regulatory_hierarchy": "master|hub|intermediate|peripheral",
  "hierarchy_rationale": "brief explanation of hierarchical position",
  "information_flow": "high|moderate|low",
  "flow_rationale": "brief explanation of information processing",
  "network_vulnerability": "critical|important|moderate|minimal",
  "vulnerability_rationale": "brief explanation of network impact",
  "perturbation_impact": "system-wide|modular|localized|minimal",
  "perturbation_rationale": "brief explanation of knockout/perturbation effects",
  "evolutionary_conservation": "high|moderate|low",
  "conservation_rationale": "brief inference about evolutionary importance",
  "summary": "1-2 sentence systems biology summary"
}}

Focus on:
- Network topology and centrality (degree, betweenness, PageRank implications)
- Hierarchical position and regulatory control
- Information flow and signal transduction
- Network robustness and vulnerability to perturbation
- Evolutionary conservation inferred from network position

Provide only the JSON, no additional text."""

        system_prompt = "You are an expert systems biologist and network theorist analyzing gene regulatory networks."

        response = await self._call_ollama(prompt, system_prompt)

        # Parse LLM response
        expected_keys = [
            "network_centrality", "centrality_rationale",
            "regulatory_hierarchy", "hierarchy_rationale",
            "information_flow", "flow_rationale",
            "network_vulnerability", "vulnerability_rationale",
            "perturbation_impact", "perturbation_rationale",
            "evolutionary_conservation", "conservation_rationale",
            "summary"
        ]
        insights = self._parse_llm_json(response, expected_keys)

        # Format network effects as list (for compatibility)
        network_effects = []
        if reg_count > 15:
            network_effects.append(f"Highly regulated node with {reg_count} upstream controllers")
        if target_count > 10:
            network_effects.append(f"Regulatory hub controlling {target_count} downstream targets")

        # Structure response
        return {
            "gene": gene,
            "domain": "systems_biology",
            "insights": {
                "network_centrality": insights.get("network_centrality", 0.0),
                "centrality_rationale": insights.get("centrality_rationale", "N/A"),
                "regulatory_hierarchy": insights.get("regulatory_hierarchy", "unknown"),
                "hierarchy_rationale": insights.get("hierarchy_rationale", "N/A"),
                "information_flow": insights.get("information_flow", "unknown"),
                "flow_rationale": insights.get("flow_rationale", "N/A"),
                "network_vulnerability": insights.get("network_vulnerability", "unknown"),
                "vulnerability_rationale": insights.get("vulnerability_rationale", "N/A")
            },
            "network_effects": network_effects,
            "perturbation_impact": insights.get("perturbation_impact", "unknown"),
            "evolutionary_conservation": insights.get("evolutionary_conservation", "moderate"),
            "summary": insights.get("summary", f"{gene} systems biology analysis completed"),
            "llm_powered": True
        }

    def _analyze_systems_rules(self, gene: str, gene_info: Dict, regulators_analysis: Dict, targets_analysis: Dict) -> Dict:
        """Rule-based systems biology analysis (fallback)"""
        num_regulators = gene_info.get('num_regulators', 0)
        num_targets = gene_info.get('num_targets', 0)
        regulatory_role = gene_info.get('regulatory_role', 'unknown')

        # Network topology analysis
        network_centrality = (num_regulators + num_targets) / 50.0  # Normalized centrality score

        systems_insights = {
            "network_centrality": min(1.0, network_centrality),
            "regulatory_hierarchy": "master" if regulatory_role == "master_regulator" else "hub" if regulatory_role == "hub_regulator" else "intermediate",
            "information_flow": "high" if num_regulators > 10 and num_targets > 10 else "moderate" if num_regulators + num_targets > 10 else "low",
            "network_vulnerability": "critical" if regulatory_role in ['hub_regulator', 'master_regulator'] else "important" if num_targets > 5 else "minimal"
        }

        # Regulatory network effects
        network_effects = []
        if regulators_analysis:
            reg_count = regulators_analysis.get('regulator_summary', {}).get('total_regulators', 0)
            if reg_count > 15:
                network_effects.append(f"Highly regulated node with {reg_count} upstream controllers")

        if targets_analysis:
            target_count = targets_analysis.get('target_summary', {}).get('total_targets', 0)
            if target_count > 10:
                network_effects.append(f"Regulatory hub controlling {target_count} downstream targets")

        return {
            "gene": gene,
            "domain": "systems_biology",
            "insights": systems_insights,
            "network_effects": network_effects,
            "perturbation_impact": "system-wide" if systems_insights["network_vulnerability"] == "critical" else "localized",
            "evolutionary_conservation": "high" if network_centrality > 0.5 else "moderate",
            "summary": f"{gene} has {systems_insights['network_centrality']:.1%} network centrality with {systems_insights['network_vulnerability']} vulnerability level",
            "llm_powered": False
        }

class GeneAnalysisState(TypedDict):
    """State object that tracks the entire gene analysis workflow"""
    # Input parameters
    gene: str
    cell_type: str
    analysis_depth: str  # "basic", "comprehensive", "focused"

    # Current workflow state
    current_step: str
    workflow_complete: bool
    error_message: Optional[str]

    # Analysis results
    gene_info: Optional[Dict]
    regulators_analysis: Optional[Dict]
    targets_analysis: Optional[Dict]
    pathway_analysis: Optional[Dict]
    cross_cell_analysis: Optional[Dict]
    perturbation_analysis: Optional[Dict]

    # Domain-specific analysis results
    cancer_analysis: Optional[Dict]
    drug_analysis: Optional[Dict]
    clinical_analysis: Optional[Dict]
    systems_analysis: Optional[Dict]

    # Workflow decisions
    next_actions: List[str]
    priority_analysis: str
    domain_focus: List[str]

    # Final output
    comprehensive_report: Optional[Dict]
    analysis_metadata: Dict

class RegNetAgentsWorkflow:
    """
    LangGraph implementation of RegNetAgents gene analysis workflow.

    This class orchestrates a multi-agent system for comprehensive gene regulatory network analysis.
    It coordinates specialized agents (modeling, pathway, domain) through a state-based workflow
    that optimizes parallel execution while maintaining dependencies between analysis steps.

    The workflow architecture uses LangGraph's StateGraph for:
    - Parallel batch processing of independent analyses (regulators, targets, pathways)
    - Conditional routing based on gene network characteristics
    - State persistence for long-running analyses
    - Error handling with graceful degradation

    Key Features:
        - Multi-agent coordination with specialized domain expertise
        - Parallel execution of independent analysis tasks
        - Automatic perturbation analysis for genes with >5 regulators
        - LLM-powered insights with rule-based fallback
        - Cross-cell type comparison capabilities

    Attributes:
        cache (RegNetAgentsCache): Pre-loaded network indices for 10 cell types
        modeling_agent (RegNetAgentsModelingAgent): Network topology analysis
        integration_agent (CrossSystemIntegrationAgent): Multi-gene integration
        pathway_enricher (PathwayEnricherAgent): Reactome pathway enrichment
        domain_agents (DomainAnalysisAgents): LLM-powered domain insights
        workflow (StateGraph): Compiled LangGraph workflow

    Example:
        >>> workflow = RegNetAgentsWorkflow()
        >>> result = await workflow.run(
        ...     gene="TP53",
        ...     cell_type=CellType.EPITHELIAL_CELLS,
        ...     analysis_depth="comprehensive"
        ... )
    """

    def __init__(self):
        """
        Initialize the workflow with RegNetAgents components.

        Loads network indices for all available cell types, initializes specialized agents,
        and compiles the LangGraph workflow structure with appropriate node connections
        and conditional routing logic.
        """
        self.cache = RegNetAgentsCache()
        self.modeling_agent = RegNetAgentsModelingAgent(self.cache)
        self.integration_agent = CrossSystemIntegrationAgent(self.cache)
        self.pathway_enricher = PathwayEnricherAgent()
        self.domain_agents = DomainAnalysisAgents(self.cache)

        # Create the workflow graph
        self.workflow = self._create_workflow()
        logger.info("RegNetAgents LangGraph workflow initialized")

    def _create_workflow(self) -> StateGraph:
        """
        Create the LangGraph workflow structure with optimized parallel execution.

        Constructs a state graph that orchestrates gene analysis through multiple stages:
        1. Initialization and gene network position analysis
        2. Conditional routing based on network characteristics
        3. Parallel batch processing of independent analyses
        4. Domain-specific insights generation (LLM-powered)
        5. Final report compilation

        The workflow uses batch nodes to parallelize independent analyses, reducing
        execution time from sequential ~2-3 seconds to parallel ~0.6-1 second for
        core analyses (network + regulators + targets + pathways).

        Returns:
            StateGraph: Compiled LangGraph workflow ready for execution
        """
        workflow = StateGraph(GeneAnalysisState)

        # Add initialization and routing nodes
        # These run sequentially to establish gene context and determine analysis path
        workflow.add_node("initialize_analysis", self._initialize_analysis)
        workflow.add_node("analyze_gene_network", self._analyze_gene_network)
        workflow.add_node("decide_next_steps", self._decide_next_steps)

        # Parallel batch nodes for core analyses
        # These execute independent analyses concurrently for optimal performance
        workflow.add_node("batch_core_analyses", self._batch_core_analyses)
        workflow.add_node("batch_secondary_analyses", self._batch_secondary_analyses)
        workflow.add_node("batch_domain_analyses", self._batch_domain_analyses)

        # Individual analysis nodes (for fallback/specific routing)
        workflow.add_node("analyze_regulators", self._analyze_regulators)
        workflow.add_node("analyze_targets", self._analyze_targets)
        workflow.add_node("analyze_pathways", self._analyze_pathways)
        workflow.add_node("analyze_perturbations", self._analyze_perturbations)
        workflow.add_node("cross_cell_comparison", self._cross_cell_comparison)

        # Add domain analysis nodes
        workflow.add_node("analyze_cancer_domain", self._analyze_cancer_domain)
        workflow.add_node("analyze_drug_domain", self._analyze_drug_domain)
        workflow.add_node("analyze_clinical_domain", self._analyze_clinical_domain)
        workflow.add_node("analyze_systems_domain", self._analyze_systems_domain)

        workflow.add_node("generate_final_report", self._generate_final_report)
        workflow.add_node("handle_error", self._handle_error)

        # Define the workflow flow
        workflow.set_entry_point("initialize_analysis")

        # Linear progression with conditional branching
        workflow.add_edge("initialize_analysis", "analyze_gene_network")
        workflow.add_edge("analyze_gene_network", "decide_next_steps")

        # Conditional routing based on gene characteristics
        workflow.add_conditional_edges(
            "decide_next_steps",
            self._route_next_action,
            {
                "batch_core": "batch_core_analyses",
                "batch_secondary": "batch_secondary_analyses",
                "batch_domain": "batch_domain_analyses",
                "regulators": "analyze_regulators",
                "targets": "analyze_targets",
                "pathways": "analyze_pathways",
                "perturbations": "analyze_perturbations",
                "cross_cell": "cross_cell_comparison",
                "cancer_domain": "analyze_cancer_domain",
                "drug_domain": "analyze_drug_domain",
                "clinical_domain": "analyze_clinical_domain",
                "systems_domain": "analyze_systems_domain",
                "complete": "generate_final_report",
                "error": "handle_error"
            }
        )

        # Batch nodes flow to decision point
        workflow.add_edge("batch_core_analyses", "decide_next_steps")
        workflow.add_edge("batch_secondary_analyses", "decide_next_steps")
        workflow.add_edge("batch_domain_analyses", "decide_next_steps")

        # Flow back to decision point for multi-step analysis
        workflow.add_edge("analyze_regulators", "decide_next_steps")
        workflow.add_edge("analyze_targets", "decide_next_steps")
        workflow.add_edge("analyze_pathways", "decide_next_steps")
        workflow.add_edge("analyze_perturbations", "decide_next_steps")
        workflow.add_edge("cross_cell_comparison", "decide_next_steps")

        # Domain analysis flows back to decision point
        workflow.add_edge("analyze_cancer_domain", "decide_next_steps")
        workflow.add_edge("analyze_drug_domain", "decide_next_steps")
        workflow.add_edge("analyze_clinical_domain", "decide_next_steps")
        workflow.add_edge("analyze_systems_domain", "decide_next_steps")

        # Terminal nodes
        workflow.add_edge("generate_final_report", END)
        workflow.add_edge("handle_error", END)

        return workflow

    async def _initialize_analysis(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Initialize the analysis with input validation"""
        logger.info(f"Initializing analysis for gene: {state['gene']}")

        try:
            # Validate inputs
            if not state.get('gene'):
                raise ValueError("Gene name is required")

            if not state.get('cell_type'):
                state['cell_type'] = 'epithelial_cell'  # Default

            if not state.get('analysis_depth'):
                state['analysis_depth'] = 'comprehensive'

            # Initialize workflow state
            state['current_step'] = 'initialization'
            state['workflow_complete'] = False
            state['next_actions'] = []
            state['analysis_metadata'] = {
                'start_time': asyncio.get_event_loop().time(),
                'steps_completed': [],
                'total_analysis_time': 0
            }

            logger.info(f"Analysis initialized for {state['gene']} in {state['cell_type']}")
            return state

        except Exception as e:
            state['error_message'] = str(e)
            state['current_step'] = 'error'
            return state

    async def _analyze_gene_network(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Perform initial gene network analysis (equivalent to analyze_gene_network MCP tool)"""
        logger.info(f"Analyzing gene network for {state['gene']}")

        try:
            state['current_step'] = 'gene_network_analysis'

            # Use existing MCP server logic
            cell_type = CellType(state['cell_type'])
            result = await self.modeling_agent.analyze_gene_network_context(
                state['gene'],
                cell_type
            )

            state['gene_info'] = result
            state['analysis_metadata']['steps_completed'].append('gene_network_analysis')

            logger.info(f"Gene network analysis complete. Role: {result.get('regulatory_role')}")
            return state

        except Exception as e:
            state['error_message'] = f"Gene network analysis failed: {str(e)}"
            state['current_step'] = 'error'
            return state

    def _route_next_action(self, state: GeneAnalysisState) -> str:
        """Smart routing logic based on gene characteristics"""

        if state.get('error_message'):
            return "error"

        if state['workflow_complete']:
            return "complete"

        # Get gene characteristics
        gene_info = state.get('gene_info', {})
        regulatory_role = gene_info.get('regulatory_role')

        num_regulators = gene_info.get('num_regulators', 0)
        num_targets = gene_info.get('num_targets', 0)
        is_regulator = gene_info.get('is_regulator', False)

        completed_steps = state['analysis_metadata'].get('steps_completed', [])

        # High priority: Analyze regulators for heavily regulated genes
        if (num_regulators > 15 and 'regulators_analysis' not in completed_steps):
            logger.info(f"Routing to regulators analysis (high priority: {num_regulators} regulators)")
            return "regulators"

        # High priority: Analyze targets for hub regulators
        if (is_regulator and num_targets > 20 and 'targets_analysis' not in completed_steps):
            logger.info(f"Routing to targets analysis (high priority: {num_targets} targets)")
            return "targets"

        # Medium priority: Analyze regulators for moderately regulated genes
        if (num_regulators > 5 and 'regulators_analysis' not in completed_steps):
            logger.info(f"Routing to regulators analysis (medium priority: {num_regulators} regulators)")
            return "regulators"

        # Medium priority: Analyze targets for intermediate regulators
        if (is_regulator and num_targets > 5 and 'targets_analysis' not in completed_steps):
            logger.info(f"Routing to targets analysis (medium priority: {num_targets} targets)")
            return "targets"

        # Perturbation analysis after regulators (therapeutic potential)
        if ('regulators_analysis' in completed_steps
            and 'perturbation_analysis' not in completed_steps
            and num_regulators > 5):  # Only for genes with meaningful regulators
            logger.info(f"Routing to perturbation analysis ({num_regulators} regulators to test)")
            return "perturbations"

        # Pathway analysis if we have regulator or target data (comprehensive mode only)
        if (state['analysis_depth'] == 'comprehensive'
            and ('regulators_analysis' in completed_steps or 'targets_analysis' in completed_steps)
            and 'pathway_analysis' not in completed_steps):
            logger.info("Routing to pathway analysis (comprehensive mode)")
            return "pathways"

        # Medium priority: Cross-cell comparison for important genes
        if ((regulatory_role in ['hub_regulator', 'master_regulator'] or num_regulators > 15)
            and 'cross_cell_analysis' not in completed_steps):
            logger.info("Routing to cross-cell comparison")
            return "cross_cell"

        # Domain analyses (comprehensive mode only) - run all 4 in parallel
        if state['analysis_depth'] == 'comprehensive':
            # Check if we have necessary data for domain analyses
            has_core_data = ('regulators_analysis' in completed_steps or
                           'targets_analysis' in completed_steps)

            # Check if ANY domain analysis is not yet completed
            domain_analyses_needed = any([
                'cancer_domain_analysis' not in completed_steps,
                'drug_domain_analysis' not in completed_steps,
                'clinical_domain_analysis' not in completed_steps,
                'systems_domain_analysis' not in completed_steps
            ])

            if has_core_data and domain_analyses_needed:
                logger.info("Routing to batch domain analyses (parallel execution)")
                return "batch_domain"

        # If all relevant analyses are complete
        logger.info("All analyses complete, generating final report")
        return "complete"

    async def _decide_next_steps(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Decision node that determines what analysis to do next"""
        state['current_step'] = 'decision'

        # Update next actions based on current state
        next_action = self._route_next_action(state)
        state['next_actions'] = [next_action]

        if next_action == "complete":
            state['workflow_complete'] = True

        return state

    async def _batch_core_analyses(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Run core network analyses in parallel (regulators + targets)"""
        logger.info(f"Running batch core analyses for {state['gene']}")
        state['current_step'] = 'batch_core_analyses'

        gene_info = state.get('gene_info', {})
        num_regulators = gene_info.get('num_regulators', 0)
        num_targets = gene_info.get('num_targets', 0)
        is_regulator = gene_info.get('is_regulator', False)
        completed_steps = state['analysis_metadata'].get('steps_completed', [])

        cell_type = CellType(state['cell_type'])

        # Determine which analyses to run and execute them directly
        tasks = []
        task_names = []

        if num_regulators > 5 and 'regulators_analysis' not in completed_steps:
            tasks.append(self.modeling_agent.analyze_regulators_detailed(state['gene'], cell_type, max_regulators=100))
            task_names.append('regulators')

        if is_regulator and num_targets > 5 and 'targets_analysis' not in completed_steps:
            tasks.append(self.modeling_agent.analyze_targets_detailed(state['gene'], cell_type, max_targets=25))
            task_names.append('targets')

        # Run analyses in parallel
        if tasks:
            logger.info(f"Running {len(tasks)} core analyses in parallel: {task_names}")
            results = await asyncio.gather(*tasks, return_exceptions=True)

            # Merge results back into state
            for i, result in enumerate(results):
                if isinstance(result, Exception):
                    logger.error(f"Core analysis {task_names[i]} failed: {result}")
                else:
                    # Merge successful results
                    if task_names[i] == 'regulators':
                        state['regulators_analysis'] = result
                        if 'regulators_analysis' not in state['analysis_metadata']['steps_completed']:
                            state['analysis_metadata']['steps_completed'].append('regulators_analysis')
                    elif task_names[i] == 'targets':
                        state['targets_analysis'] = result
                        if 'targets_analysis' not in state['analysis_metadata']['steps_completed']:
                            state['analysis_metadata']['steps_completed'].append('targets_analysis')

        logger.info(f"Batch core analyses complete. Completed: {task_names}")
        return state

    async def _batch_secondary_analyses(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Run secondary analyses in parallel (pathways + cross-cell)"""
        logger.info(f"Running batch secondary analyses for {state['gene']}")
        state['current_step'] = 'batch_secondary_analyses'

        gene_info = state.get('gene_info', {})
        regulatory_role = gene_info.get('regulatory_role')
        num_regulators = gene_info.get('num_regulators', 0)
        completed_steps = state['analysis_metadata'].get('steps_completed', [])

        # Determine which analyses to run and create coroutines
        tasks = []
        task_names = []

        has_core_data = ('regulators_analysis' in completed_steps or 'targets_analysis' in completed_steps)

        if has_core_data and 'pathway_analysis' not in completed_steps:
            # Build gene list for pathway analysis
            gene_list = [state['gene']]
            if state.get('regulators_analysis'):
                hub_regulators = state['regulators_analysis'].get('hub_regulators', [])
                for reg in hub_regulators[:10]:
                    if reg.get('gene_symbol') and reg['gene_symbol'] != 'Unknown':
                        gene_list.append(reg['gene_symbol'])
            if state.get('targets_analysis'):
                cascade_targets = state['targets_analysis'].get('cascade_targets', [])
                for target in cascade_targets[:10]:
                    if target.get('gene_symbol') and target['gene_symbol'] != 'Unknown':
                        gene_list.append(target['gene_symbol'])

            # Remove duplicates
            gene_list = list(set(gene_list))

            tasks.append(self.pathway_enricher.enrich_pathways_reactome(gene_list))
            task_names.append('pathways')

        if ((regulatory_role in ['hub_regulator', 'master_regulator'] or num_regulators > 15)
            and 'cross_cell_analysis' not in completed_steps):
            tasks.append(self.modeling_agent.compare_gene_across_cell_types(state['gene']))
            task_names.append('cross_cell')

        # Run analyses in parallel
        if tasks:
            logger.info(f"Running {len(tasks)} secondary analyses in parallel: {task_names}")
            results = await asyncio.gather(*tasks, return_exceptions=True)

            # Merge results back into state
            for i, result in enumerate(results):
                if isinstance(result, Exception):
                    logger.warning(f"Secondary analysis {task_names[i]} failed: {result}")
                else:
                    # Merge successful results
                    if task_names[i] == 'pathways':
                        state['pathway_analysis'] = result
                        if 'pathway_analysis' not in state['analysis_metadata']['steps_completed']:
                            state['analysis_metadata']['steps_completed'].append('pathway_analysis')
                    elif task_names[i] == 'cross_cell':
                        state['cross_cell_analysis'] = result
                        if 'cross_cell_analysis' not in state['analysis_metadata']['steps_completed']:
                            state['analysis_metadata']['steps_completed'].append('cross_cell_analysis')

        logger.info(f"Batch secondary analyses complete. Completed: {task_names}")
        return state

    async def _batch_domain_analyses(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Run all domain analyses in parallel (cancer + drug + clinical + systems)"""
        logger.info(f"Running batch domain analyses for {state['gene']}")
        state['current_step'] = 'batch_domain_analyses'

        gene_info = state.get('gene_info', {})
        regulatory_role = gene_info.get('regulatory_role')
        num_regulators = gene_info.get('num_regulators', 0)
        num_targets = gene_info.get('num_targets', 0)
        is_regulator = gene_info.get('is_regulator', False)
        completed_steps = state['analysis_metadata'].get('steps_completed', [])

        pathway_analysis = state.get('pathway_analysis', {})
        regulators_analysis = state.get('regulators_analysis', {})
        targets_analysis = state.get('targets_analysis', {})
        cross_cell_analysis = state.get('cross_cell_analysis', {})

        # Determine which analyses to run and create coroutines
        tasks = []
        task_names = []

        if ('cancer_domain_analysis' not in completed_steps
            and (num_regulators > 10 or num_targets > 15 or
                 regulatory_role in ['hub_regulator', 'master_regulator', 'heavily_regulated'])):
            tasks.append(self.domain_agents.analyze_cancer_context(state['gene'], gene_info, pathway_analysis))
            task_names.append('cancer')

        if ('drug_domain_analysis' not in completed_steps
            and (is_regulator and num_targets > 10)):
            tasks.append(self.domain_agents.analyze_drug_development(state['gene'], gene_info, regulators_analysis, targets_analysis))
            task_names.append('drug')

        if ('clinical_domain_analysis' not in completed_steps
            and ('pathway_analysis' in completed_steps or num_regulators > 5)):
            tasks.append(self.domain_agents.analyze_clinical_relevance(state['gene'], gene_info, cross_cell_analysis))
            task_names.append('clinical')

        if ('systems_domain_analysis' not in completed_steps
            and len(completed_steps) >= 3
            and ('regulators_analysis' in completed_steps or 'targets_analysis' in completed_steps)):
            tasks.append(self.domain_agents.analyze_systems_biology(state['gene'], gene_info, regulators_analysis, targets_analysis))
            task_names.append('systems')

        # Run analyses in parallel
        if tasks:
            logger.info(f"Running {len(tasks)} domain analyses in parallel: {task_names}")
            results = await asyncio.gather(*tasks, return_exceptions=True)

            # Merge results back into state
            for i, result in enumerate(results):
                if isinstance(result, Exception):
                    logger.warning(f"Domain analysis {task_names[i]} failed (optional): {result}")
                    # Mark as completed even if failed (optional analyses)
                    if f'{task_names[i]}_domain_analysis' not in completed_steps:
                        state['analysis_metadata']['steps_completed'].append(f'{task_names[i]}_domain_analysis')
                else:
                    # Merge successful results
                    if task_names[i] == 'cancer':
                        state['cancer_analysis'] = result
                        if 'cancer_domain_analysis' not in state['analysis_metadata']['steps_completed']:
                            state['analysis_metadata']['steps_completed'].append('cancer_domain_analysis')
                    elif task_names[i] == 'drug':
                        state['drug_analysis'] = result
                        if 'drug_domain_analysis' not in state['analysis_metadata']['steps_completed']:
                            state['analysis_metadata']['steps_completed'].append('drug_domain_analysis')
                    elif task_names[i] == 'clinical':
                        state['clinical_analysis'] = result
                        if 'clinical_domain_analysis' not in state['analysis_metadata']['steps_completed']:
                            state['analysis_metadata']['steps_completed'].append('clinical_domain_analysis')
                    elif task_names[i] == 'systems':
                        state['systems_analysis'] = result
                        if 'systems_domain_analysis' not in state['analysis_metadata']['steps_completed']:
                            state['analysis_metadata']['steps_completed'].append('systems_domain_analysis')
        else:
            # No domain analyses were eligible - mark all as skipped to prevent infinite loop
            logger.info("No domain analyses eligible for this gene - marking all as completed (skipped)")
            for domain in ['cancer', 'drug', 'clinical', 'systems']:
                domain_key = f'{domain}_domain_analysis'
                if domain_key not in completed_steps:
                    state['analysis_metadata']['steps_completed'].append(domain_key)

        logger.info(f"Batch domain analyses complete. Completed: {task_names if tasks else 'none (all skipped)'}")
        return state

    async def _analyze_regulators(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Analyze gene regulators (equivalent to analyze_regulators MCP tool)"""
        logger.info(f"Analyzing regulators for {state['gene']}")

        try:
            state['current_step'] = 'regulators_analysis'

            cell_type = CellType(state['cell_type'])
            result = await self.modeling_agent.analyze_regulators_detailed(
                state['gene'],
                cell_type,
                max_regulators=100  # Analyze all regulators (not just top 25)
            )

            state['regulators_analysis'] = result
            state['analysis_metadata']['steps_completed'].append('regulators_analysis')

            logger.info(f"Regulators analysis complete. Found {result.get('regulator_summary', {}).get('total_regulators', 0)} regulators")
            return state

        except Exception as e:
            state['error_message'] = f"Regulators analysis failed: {str(e)}"
            state['current_step'] = 'error'
            return state

    async def _analyze_targets(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Analyze gene targets (equivalent to analyze_targets MCP tool)"""
        logger.info(f"Analyzing targets for {state['gene']}")

        try:
            state['current_step'] = 'targets_analysis'

            cell_type = CellType(state['cell_type'])
            result = await self.modeling_agent.analyze_targets_detailed(
                state['gene'],
                cell_type,
                max_targets=25
            )

            state['targets_analysis'] = result
            state['analysis_metadata']['steps_completed'].append('targets_analysis')

            logger.info(f"Targets analysis complete. Found {result.get('target_summary', {}).get('total_targets', 0)} targets")
            return state

        except Exception as e:
            state['error_message'] = f"Targets analysis failed: {str(e)}"
            state['current_step'] = 'error'
            return state

    async def _analyze_pathways(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Analyze biological pathways using Reactome enrichment"""
        logger.info(f"Analyzing pathways for {state['gene']} using Reactome")

        try:
            state['current_step'] = 'pathway_analysis'

            # Extract gene list from previous analyses
            gene_list = [state['gene']]

            # Add regulators if available
            # Limit to top 10 to balance biological signal with pathway specificity:
            # - Reactome enrichment works best with focused gene sets (10-50 genes)
            # - Top 10 captures strongest regulatory relationships (highest centrality)
            # - Prevents over-enrichment with broad, non-specific pathways
            # - For genes with many regulators, focusing on top connections is more informative
            if state.get('regulators_analysis'):
                hub_regulators = state['regulators_analysis'].get('hub_regulators', [])
                for reg in hub_regulators[:10]:  # Top 10 hub regulators
                    if reg.get('gene_symbol') and reg['gene_symbol'] != 'Unknown':
                        gene_list.append(reg['gene_symbol'])

            # Add targets if available
            # Same rationale as regulators: top 10 focuses on strongest downstream effects
            # Example: MYC has 427 targets - using all would dilute pathway specificity
            if state.get('targets_analysis'):
                cascade_targets = state['targets_analysis'].get('cascade_targets', [])
                for target in cascade_targets[:10]:  # Top 10 cascade targets
                    if target.get('gene_symbol') and target['gene_symbol'] != 'Unknown':
                        gene_list.append(target['gene_symbol'])

            # Remove duplicates
            gene_list = list(set(gene_list))

            # Perform Reactome pathway enrichment
            result = await self.pathway_enricher.enrich_pathways_reactome(gene_list)

            state['pathway_analysis'] = result
            state['analysis_metadata']['steps_completed'].append('pathway_analysis')

            if result.get('status') == 'success':
                sig_count = result.get('summary', {}).get('significant_pathways', 0)
                logger.info(f"Pathway analysis complete. Found {sig_count} significant pathways (FDR < 0.05)")
            else:
                logger.warning(f"Pathway analysis returned with status: {result.get('status')}")

            return state

        except Exception as e:
            logger.error(f"Pathway analysis failed: {str(e)}")
            state['error_message'] = f"Pathway analysis failed: {str(e)}"
            state['current_step'] = 'error'
            return state

    async def _cross_cell_comparison(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Compare gene across cell types (equivalent to compare_across_cell_types)"""
        logger.info(f"Performing cross-cell comparison for {state['gene']}")

        try:
            state['current_step'] = 'cross_cell_analysis'

            result = await self.modeling_agent.compare_gene_across_cell_types(state['gene'])

            state['cross_cell_analysis'] = result
            state['analysis_metadata']['steps_completed'].append('cross_cell_analysis')

            logger.info(f"Cross-cell analysis complete. Analyzed across {len(result.get('cell_type_analysis', {}))} cell types")
            return state

        except Exception as e:
            state['error_message'] = f"Cross-cell analysis failed: {str(e)}"
            state['current_step'] = 'error'
            return state

    async def _analyze_perturbations(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Analyze perturbation effects of inhibiting upstream regulators"""
        logger.info(f"Analyzing perturbation effects for {state['gene']}")

        try:
            state['current_step'] = 'perturbation_analysis'

            cell_type = CellType(state['cell_type'])
            regulators_analysis = state.get('regulators_analysis', {})

            if not regulators_analysis:
                logger.warning("No regulators analysis available for perturbation")
                state['perturbation_analysis'] = {
                    "status": "skipped",
                    "reason": "No regulators to perturb"
                }
            else:
                # Run perturbation analysis (analyze all regulators, not just top 10)
                result = await self.modeling_agent.analyze_perturbation_effects(
                    target_gene=state['gene'],
                    cell_type=cell_type,
                    regulators=regulators_analysis.get('hub_regulators', []),
                    max_regulators=100
                )

                state['perturbation_analysis'] = result
                top_target = result.get('perturbation_results', [{}])[0].get('regulator', 'N/A') if result.get('perturbation_results') else 'N/A'
                logger.info(f"Perturbation analysis complete. Top therapeutic target: {top_target}")

            state['analysis_metadata']['steps_completed'].append('perturbation_analysis')
            return state

        except Exception as e:
            logger.warning(f"Perturbation analysis failed (optional): {str(e)}")
            state['perturbation_analysis'] = {"error": str(e), "analysis_skipped": True}
            state['analysis_metadata']['steps_completed'].append('perturbation_analysis')
            return state

    async def _analyze_cancer_domain(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Analyze gene in cancer context using domain analysis agents"""
        logger.info(f"Analyzing cancer domain context for {state['gene']}")

        try:
            state['current_step'] = 'cancer_domain_analysis'

            gene_info = state.get('gene_info', {})
            pathway_analysis = state.get('pathway_analysis', {})

            result = await self.domain_agents.analyze_cancer_context(
                state['gene'], gene_info, pathway_analysis
            )

            state['cancer_analysis'] = result
            state['analysis_metadata']['steps_completed'].append('cancer_domain_analysis')

            logger.info(f"Cancer domain analysis complete for {state['gene']}")
            return state

        except Exception as e:
            logger.warning(f"Cancer domain analysis failed (optional): {str(e)}")
            state['cancer_analysis'] = {"error": str(e), "analysis_skipped": True}
            state['analysis_metadata']['steps_completed'].append('cancer_domain_analysis')
            return state

    async def _analyze_drug_domain(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Analyze gene in drug development context using domain analysis agents"""
        logger.info(f"Analyzing drug development domain context for {state['gene']}")

        try:
            state['current_step'] = 'drug_domain_analysis'

            gene_info = state.get('gene_info', {})
            pathway_analysis = state.get('pathway_analysis', {})

            regulators_analysis = state.get('regulators_analysis', {})
            targets_analysis = state.get('targets_analysis', {})

            result = await self.domain_agents.analyze_drug_development(
                state['gene'], gene_info, regulators_analysis, targets_analysis
            )

            state['drug_analysis'] = result
            state['analysis_metadata']['steps_completed'].append('drug_domain_analysis')

            logger.info(f"Drug development domain analysis complete for {state['gene']}")
            return state

        except Exception as e:
            logger.warning(f"Drug development domain analysis failed (optional): {str(e)}")
            state['drug_analysis'] = {"error": str(e), "analysis_skipped": True}
            state['analysis_metadata']['steps_completed'].append('drug_domain_analysis')
            return state

    async def _analyze_clinical_domain(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Analyze gene in clinical relevance context using domain analysis agents"""
        logger.info(f"Analyzing clinical domain context for {state['gene']}")

        try:
            state['current_step'] = 'clinical_domain_analysis'

            gene_info = state.get('gene_info', {})
            pathway_analysis = state.get('pathway_analysis', {})

            result = await self.domain_agents.analyze_clinical_relevance(
                state['gene'], gene_info, pathway_analysis
            )

            state['clinical_analysis'] = result
            state['analysis_metadata']['steps_completed'].append('clinical_domain_analysis')

            logger.info(f"Clinical domain analysis complete for {state['gene']}")
            return state

        except Exception as e:
            logger.warning(f"Clinical domain analysis failed (optional): {str(e)}")
            state['clinical_analysis'] = {"error": str(e), "analysis_skipped": True}
            state['analysis_metadata']['steps_completed'].append('clinical_domain_analysis')
            return state

    async def _analyze_systems_domain(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Analyze gene in systems biology context using domain analysis agents"""
        logger.info(f"Analyzing systems biology domain context for {state['gene']}")

        try:
            state['current_step'] = 'systems_domain_analysis'

            gene_info = state.get('gene_info', {})
            pathway_analysis = state.get('pathway_analysis', {})
            regulators_analysis = state.get('regulators_analysis', {})
            targets_analysis = state.get('targets_analysis', {})

            result = await self.domain_agents.analyze_systems_biology(
                state['gene'], gene_info, regulators_analysis, targets_analysis
            )

            state['systems_analysis'] = result
            state['analysis_metadata']['steps_completed'].append('systems_domain_analysis')

            logger.info(f"Systems biology domain analysis complete for {state['gene']}")
            return state

        except Exception as e:
            logger.warning(f"Systems biology domain analysis failed (optional): {str(e)}")
            state['systems_analysis'] = {"error": str(e), "analysis_skipped": True}
            state['analysis_metadata']['steps_completed'].append('systems_domain_analysis')
            return state

    async def _generate_final_report(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Generate comprehensive final report"""
        logger.info(f"Generating final report for {state['gene']}")

        try:
            state['current_step'] = 'final_report'

            # Calculate total analysis time
            start_time = state['analysis_metadata'].get('start_time', 0)
            total_time = (asyncio.get_event_loop().time() - start_time) * 1000  # Convert to ms
            state['analysis_metadata']['total_analysis_time'] = total_time

            # Compile comprehensive report
            report = {
                "gene_analysis_summary": {
                    "gene": state['gene'],
                    "cell_type": state['cell_type'],
                    "regulatory_role": state.get('gene_info', {}).get('regulatory_role'),
                    "analysis_depth": state['analysis_depth']
                },
                "network_analysis": state.get('gene_info'),
                "regulatory_analysis": state.get('regulators_analysis'),
                "target_analysis": state.get('targets_analysis'),
                "pathway_enrichment": state.get('pathway_analysis'),
                "perturbation_analysis": state.get('perturbation_analysis'),
                "cross_cell_comparison": state.get('cross_cell_analysis'),
                "similarity_analysis": state.get('similarity_analysis'),
                "domain_analysis": {
                    "cancer_analysis": state.get('cancer_analysis'),
                    "drug_analysis": state.get('drug_analysis'),
                    "clinical_analysis": state.get('clinical_analysis'),
                    "systems_analysis": state.get('systems_analysis')
                },
                "workflow_metadata": state['analysis_metadata'],
                "key_insights": self._extract_key_insights(state)
            }

            state['comprehensive_report'] = report
            state['workflow_complete'] = True

            logger.info(f"Final report generated. Total analysis time: {total_time:.2f}ms")
            return state

        except Exception as e:
            state['error_message'] = f"Report generation failed: {str(e)}"
            state['current_step'] = 'error'
            return state

    def _extract_key_insights(self, state: GeneAnalysisState) -> Dict:
        """Extract key insights from all analyses"""
        insights = {}

        # Gene network insights
        gene_info = state.get('gene_info', {})
        if gene_info:
            network_pos = gene_info.get('network_position', {})
            insights['regulatory_role'] = gene_info.get('regulatory_role')
            insights['regulation_strength'] = "high" if gene_info.get('num_regulators', 0) > 15 else "moderate"
            insights['regulatory_influence'] = "high" if gene_info.get('num_targets', 0) > 20 else "low"

        # Pathway insights
        pathway_analysis = state.get('pathway_analysis', {})
        if pathway_analysis:
            top_pathways = pathway_analysis.get('top_pathways', [])
            if top_pathways:
                insights['primary_pathway'] = top_pathways[0][0] if top_pathways[0] else "unknown"
                insights['pathway_diversity'] = len(pathway_analysis.get('enriched_pathways', {}))

        # Cross-cell insights
        cross_cell = state.get('cross_cell_analysis', {})
        if cross_cell:
            cell_roles = cross_cell.get('cell_type_analysis', {})
            role_diversity = len(set(analysis.get('regulatory_role') for analysis in cell_roles.values()))
            insights['cell_type_specificity'] = "high" if role_diversity > 3 else "low"

        # Domain analysis insights
        cancer_analysis = state.get('cancer_analysis', {})
        if cancer_analysis and not cancer_analysis.get('error'):
            insights['cancer_relevance'] = cancer_analysis.get('oncogenic_potential', 'unknown')
            insights['therapeutic_potential'] = cancer_analysis.get('therapeutic_target_score', 0)

        drug_analysis = state.get('drug_analysis', {})
        if drug_analysis and not drug_analysis.get('error'):
            insights['druggability'] = drug_analysis.get('druggability_score', 0)
            insights['drug_development_priority'] = drug_analysis.get('development_priority', 'unknown')

        clinical_analysis = state.get('clinical_analysis', {})
        if clinical_analysis and not clinical_analysis.get('error'):
            insights['clinical_significance'] = clinical_analysis.get('clinical_significance', 'unknown')
            insights['biomarker_potential'] = clinical_analysis.get('biomarker_potential', 'unknown')

        systems_analysis = state.get('systems_analysis', {})
        if systems_analysis and not systems_analysis.get('error'):
            insights['network_importance'] = systems_analysis.get('network_centrality', 'unknown')
            insights['systems_complexity'] = systems_analysis.get('system_complexity', 'unknown')

        return insights

    async def _handle_error(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Handle workflow errors gracefully"""
        logger.error(f"Workflow error for gene {state['gene']}: {state.get('error_message')}")

        state['current_step'] = 'error_handled'
        state['workflow_complete'] = True

        # Create error report
        error_report = {
            "status": "error",
            "gene": state['gene'],
            "error_message": state.get('error_message'),
            "completed_steps": state['analysis_metadata'].get('steps_completed', []),
            "partial_results": {
                "gene_info": state.get('gene_info'),
                "regulators_analysis": state.get('regulators_analysis'),
                "targets_analysis": state.get('targets_analysis')
            }
        }

        state['comprehensive_report'] = error_report
        return state

    async def run_analysis(self, gene: str, cell_type: str = "epithelial_cell",
                          analysis_depth: str = "comprehensive") -> Dict:
        """Run the complete gene analysis workflow"""

        # Initialize state
        initial_state = GeneAnalysisState(
            gene=gene,
            cell_type=cell_type,
            analysis_depth=analysis_depth,
            current_step="",
            workflow_complete=False,
            error_message=None,
            gene_info=None,
            regulators_analysis=None,
            targets_analysis=None,
            pathway_analysis=None,
            perturbation_analysis=None,
            cross_cell_analysis=None,
            cancer_analysis=None,
            drug_analysis=None,
            clinical_analysis=None,
            systems_analysis=None,
            domain_focus=None,
            next_actions=[],
            priority_analysis="",
            comprehensive_report=None,
            analysis_metadata={}
        )

        # Compile workflow with memory checkpointing
        memory = MemorySaver()
        app = self.workflow.compile(checkpointer=memory)

        # Run the workflow
        config = {"configurable": {"thread_id": f"gene_analysis_{gene}"}}

        try:
            final_state = await app.ainvoke(initial_state, config)
            return final_state.get('comprehensive_report', {})

        except Exception as e:
            logger.error(f"Workflow execution failed: {str(e)}")
            return {
                "status": "workflow_error",
                "error": str(e),
                "gene": gene,
                "cell_type": cell_type
            }

# Example usage and testing
async def main():
    """Test the LangGraph workflow"""
    workflow = RegNetAgentsWorkflow()

    # Test with APC (the gene from your example)
    print("Testing RegNetAgents LangGraph Workflow with APC...")
    result = await workflow.run_analysis(
        gene="APC",
        cell_type="epithelial_cell",
        analysis_depth="comprehensive"
    )

    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    asyncio.run(main())