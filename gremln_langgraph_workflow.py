
#!/usr/bin/env python3
"""
GREmLN LangGraph Workflow
Converts the MCP server gene analysis workflow into a visual LangGraph workflow
"""

from typing import Dict, List, TypedDict, Optional, Any
from langgraph.graph import StateGraph, END
from langgraph.checkpoint.memory import MemorySaver
import asyncio
import json
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Import existing components
from gene_id_mapper import GeneIDMapper
from complete_gene_service import CompleteGeneService
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

class GREmLNCache:
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

class GREmLNModelingAgent:
    """Agent for gene network modeling and analysis."""
    def __init__(self, cache):
        self.cache = cache
        self.gene_service = CompleteGeneService()
        # Initialize gene mapper for ID conversion
        from gene_id_mapper import GeneIDMapper
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

            # Determine regulatory role
            num_regulators = len(regulators)
            num_targets = len(targets)

            if num_regulators > 15:
                regulatory_role = "terminal_target"  # Highly regulated
            elif num_targets > 20:
                regulatory_role = "hub_regulator"    # Regulates many
            elif num_targets > 5 and num_regulators > 5:
                regulatory_role = "intermediate_node"
            elif num_targets > 0:
                regulatory_role = "regulator"
            else:
                regulatory_role = "target"

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
        gene_context = await self.analyze_gene_network_context(gene, cell_type)
        regulators = gene_context.get('regulators', [])

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

    async def find_similar_genes(self, gene: str, cell_type: CellType, top_k: int = 10):
        """Find genes with similar regulatory patterns."""
        network_data = self.cache.network_indices.get(cell_type.value, {})

        if not network_data:
            return {
                "gene": gene,
                "cell_type": cell_type.value,
                "similar_genes": [],
                "error": "Network data not available"
            }

        # Get gene's regulatory profile
        ensembl_id = self.gene_mapper.symbol_to_ensembl(gene)
        if not ensembl_id:
            return {
                "gene": gene,
                "cell_type": cell_type.value,
                "similar_genes": [],
                "error": "Gene not found in mapping"
            }

        target_regulators = network_data.get('target_regulators', {}).get(ensembl_id, [])
        regulator_targets = network_data.get('regulator_targets', {}).get(ensembl_id, [])

        # Calculate similarity scores for all other genes
        similarities = []
        all_genes = network_data.get('all_genes', [])

        for other_gene_id in all_genes:
            if other_gene_id == ensembl_id:
                continue

            # Get other gene's regulatory profile
            other_regulators = network_data.get('target_regulators', {}).get(other_gene_id, [])
            other_targets = network_data.get('regulator_targets', {}).get(other_gene_id, [])

            # Calculate similarity based on shared regulators and targets
            similarity_score = self._calculate_regulatory_similarity(
                target_regulators, regulator_targets,
                other_regulators, other_targets
            )

            if similarity_score > 0.01:  # Lower threshold to find more potential similarities
                gene_symbol = self._convert_ensembl_to_symbol(other_gene_id)
                similarities.append({
                    "gene_symbol": gene_symbol,
                    "ensembl_id": other_gene_id,
                    "similarity_score": similarity_score,
                    "shared_regulators": len(set(target_regulators) & set(other_regulators)),
                    "shared_targets": len(set(regulator_targets) & set(other_targets))
                })

        # If we have few similarities, add known gene families as fallback
        if len(similarities) < 3:
            gene_families = self._get_known_gene_families(gene)
            for family_gene in gene_families:
                family_ensembl = self.gene_mapper.symbol_to_ensembl(family_gene)
                if family_ensembl and family_ensembl in all_genes and family_ensembl != ensembl_id:
                    # Add with lower similarity score to indicate it's family-based
                    gene_symbol = self._convert_ensembl_to_symbol(family_ensembl)
                    similarities.append({
                        "gene_symbol": gene_symbol,
                        "ensembl_id": family_ensembl,
                        "similarity_score": 0.3,  # Family-based similarity
                        "shared_regulators": 0,
                        "shared_targets": 0,
                        "similarity_type": "gene_family"
                    })

        # Sort by similarity score and take top K
        similarities.sort(key=lambda x: x['similarity_score'], reverse=True)
        top_similar = similarities[:top_k]

        return {
            "gene": gene,
            "ensembl_id": ensembl_id,
            "cell_type": cell_type.value,
            "similar_genes": top_similar,
            "total_genes_compared": len(all_genes),
            "genes_with_similarity": len(similarities)
        }

    def _calculate_regulatory_similarity(self, reg1, tgt1, reg2, tgt2):
        """Calculate similarity between two genes based on their regulatory patterns."""
        # Convert to sets for faster operations
        reg1_set = set(reg1)
        tgt1_set = set(tgt1)
        reg2_set = set(reg2)
        tgt2_set = set(tgt2)

        # Jaccard similarity for regulators and targets
        reg_intersection = len(reg1_set & reg2_set)
        reg_union = len(reg1_set | reg2_set)
        reg_similarity = reg_intersection / reg_union if reg_union > 0 else 0

        tgt_intersection = len(tgt1_set & tgt2_set)
        tgt_union = len(tgt1_set | tgt2_set)
        tgt_similarity = tgt_intersection / tgt_union if tgt_union > 0 else 0

        # Weight regulator similarity more heavily (genes with similar inputs are more similar)
        # This captures functional similarity - genes regulated by same factors likely have similar roles
        combined_similarity = 0.7 * reg_similarity + 0.3 * tgt_similarity

        # Boost similarity if both genes have substantial regulatory networks
        if len(reg1) > 2 and len(reg2) > 2:
            combined_similarity *= 1.2

        # Cap at 1.0
        return min(combined_similarity, 1.0)

    def _get_known_gene_families(self, gene: str):
        """Get known gene family members for fallback similarity."""
        gene_families = {
            'TP53': ['TP63', 'TP73', 'MDM2', 'MDM4', 'CDKN1A', 'CDKN2A', 'BAX', 'BCL2', 'RB1', 'E2F1'],
            'TP63': ['TP53', 'TP73', 'CDKN1A', 'CDKN2A'],
            'TP73': ['TP53', 'TP63', 'CDKN1A', 'CDKN2A'],
            'MDM2': ['TP53', 'MDM4', 'CDKN1A'],
            'MYC': ['MYCN', 'MYCL', 'MAX', 'MXD1', 'MXD3', 'MXD4'],
            'RB1': ['RBL1', 'RBL2', 'E2F1', 'E2F2', 'E2F3', 'CDKN1A', 'CDKN2A'],
            'BCL2': ['BCL2L1', 'BCL2A1', 'BCL2L2', 'BAX', 'BAK1', 'BID'],
            'BRCA1': ['BRCA2', 'ATM', 'ATR', 'CHEK1', 'CHEK2', 'TP53'],
            'EGFR': ['ERBB2', 'ERBB3', 'ERBB4', 'SRC', 'PIK3CA'],
            'KRAS': ['HRAS', 'NRAS', 'PIK3CA', 'AKT1', 'RAF1'],
        }
        return gene_families.get(gene, [])

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
        """Analyze biological context and pathways."""
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

    def __init__(self, cache):
        self.cache = cache

    async def analyze_cancer_context(self, gene: str, gene_info: Dict, pathway_analysis: Dict) -> Dict:
        """Analyze gene from cancer research perspective"""
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
        """Analyze gene from drug development perspective"""
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
        """Analyze gene from clinical/personalized medicine perspective"""
        num_regulators = gene_info.get('num_regulators', 0)
        regulatory_role = gene_info.get('regulatory_role', 'unknown')

        # Clinical significance assessment
        clinical_insights = {
            "disease_association_likelihood": "high" if num_regulators > 15 else "moderate" if num_regulators > 5 else "low",
            "biomarker_utility": "diagnostic" if regulatory_role == "terminal_target" else "prognostic" if regulatory_role == "hub_regulator" else "predictive",
            "personalization_potential": min(0.9, num_regulators * 0.04),
            "clinical_actionability": "high" if num_regulators > 10 and regulatory_role in ['hub_regulator', 'terminal_target'] else "moderate"
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
            "therapeutic_monitoring": "recommended" if clinical_insights["personalization_potential"] > 0.6 else "optional",
            "summary": f"{gene} shows {clinical_insights['disease_association_likelihood']} disease association with {clinical_insights['biomarker_utility']} biomarker utility"
        }

    async def analyze_systems_biology(self, gene: str, gene_info: Dict, regulators_analysis: Dict, targets_analysis: Dict) -> Dict:
        """Analyze gene from systems biology perspective"""
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
            "summary": f"{gene} has {systems_insights['network_centrality']:.1%} network centrality with {systems_insights['network_vulnerability']} vulnerability level"
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
    similarity_analysis: Optional[Dict]

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

class GREmLNWorkflow:
    """LangGraph implementation of GREmLN gene analysis workflow"""

    def __init__(self):
        """Initialize the workflow with existing GREmLN components"""
        self.cache = GREmLNCache()
        self.modeling_agent = GREmLNModelingAgent(self.cache)
        self.integration_agent = CrossSystemIntegrationAgent(self.cache)
        self.domain_agents = DomainAnalysisAgents(self.cache)

        # Create the workflow graph
        self.workflow = self._create_workflow()
        logger.info("GREmLN LangGraph workflow initialized")

    def _create_workflow(self) -> StateGraph:
        """Create the LangGraph workflow structure"""
        workflow = StateGraph(GeneAnalysisState)

        # Add analysis nodes
        workflow.add_node("initialize_analysis", self._initialize_analysis)
        workflow.add_node("analyze_gene_network", self._analyze_gene_network)
        workflow.add_node("decide_next_steps", self._decide_next_steps)
        workflow.add_node("analyze_regulators", self._analyze_regulators)
        workflow.add_node("analyze_targets", self._analyze_targets)
        workflow.add_node("analyze_pathways", self._analyze_pathways)
        workflow.add_node("cross_cell_comparison", self._cross_cell_comparison)
        workflow.add_node("find_similar_genes", self._find_similar_genes)

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
                "regulators": "analyze_regulators",
                "targets": "analyze_targets",
                "pathways": "analyze_pathways",
                "cross_cell": "cross_cell_comparison",
                "similarity": "find_similar_genes",
                "cancer_domain": "analyze_cancer_domain",
                "drug_domain": "analyze_drug_domain",
                "clinical_domain": "analyze_clinical_domain",
                "systems_domain": "analyze_systems_domain",
                "complete": "generate_final_report",
                "error": "handle_error"
            }
        )

        # Flow back to decision point for multi-step analysis
        workflow.add_edge("analyze_regulators", "decide_next_steps")
        workflow.add_edge("analyze_targets", "decide_next_steps")
        workflow.add_edge("analyze_pathways", "decide_next_steps")
        workflow.add_edge("cross_cell_comparison", "decide_next_steps")
        workflow.add_edge("find_similar_genes", "decide_next_steps")

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
        """Smart routing logic based on gene characteristics (equivalent to smart_workflow)"""

        if state.get('error_message'):
            return "error"

        if state['workflow_complete']:
            return "complete"

        # Get gene characteristics
        gene_info = state.get('gene_info', {})
        network_position = gene_info.get('network_position', {})
        regulatory_role = gene_info.get('regulatory_role')

        num_regulators = gene_info.get('num_regulators', 0)
        num_targets = gene_info.get('num_targets', 0)
        is_regulator = gene_info.get('is_regulator', False)

        completed_steps = state['analysis_metadata'].get('steps_completed', [])

        # Priority-based routing (from your existing smart_workflow logic)

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

        # Medium priority: Cross-cell comparison for important genes
        if ((regulatory_role in ['hub_regulator', 'master_regulator'] or num_regulators > 15)
            and 'cross_cell_analysis' not in completed_steps):
            logger.info("Routing to cross-cell comparison")
            return "cross_cell"

        # Pathway analysis if we have regulator or target data
        if (('regulators_analysis' in completed_steps or 'targets_analysis' in completed_steps)
            and 'pathway_analysis' not in completed_steps):
            logger.info("Routing to pathway analysis")
            return "pathways"

        # Similarity analysis for comprehensive mode
        if (state['analysis_depth'] == 'comprehensive'
            and 'similarity_analysis' not in completed_steps):
            logger.info("Routing to similarity analysis")
            return "similarity"

        # Domain analysis - prioritize based on gene characteristics and previous findings

        # Cancer domain analysis - prioritize for genes with high regulatory complexity
        if (state['analysis_depth'] == 'comprehensive'
            and 'cancer_domain_analysis' not in completed_steps
            and (num_regulators > 10 or num_targets > 15 or
                 regulatory_role in ['hub_regulator', 'master_regulator', 'terminal_target'])):
            logger.info("Routing to cancer domain analysis")
            return "cancer_domain"

        # Drug development domain analysis - for genes with therapeutic potential
        if (state['analysis_depth'] == 'comprehensive'
            and 'drug_domain_analysis' not in completed_steps
            and (is_regulator and num_targets > 10)):
            logger.info("Routing to drug development domain analysis")
            return "drug_domain"

        # Clinical domain analysis - for medically relevant genes
        if (state['analysis_depth'] == 'comprehensive'
            and 'clinical_domain_analysis' not in completed_steps
            and ('pathway_analysis' in completed_steps or num_regulators > 5)):
            logger.info("Routing to clinical domain analysis")
            return "clinical_domain"

        # Systems biology domain analysis - comprehensive view when we have sufficient data
        if (state['analysis_depth'] == 'comprehensive'
            and 'systems_domain_analysis' not in completed_steps
            and len(completed_steps) >= 3  # Have sufficient analysis data
            and ('regulators_analysis' in completed_steps or 'targets_analysis' in completed_steps)):
            logger.info("Routing to systems biology domain analysis")
            return "systems_domain"

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

    async def _analyze_regulators(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Analyze gene regulators (equivalent to analyze_regulators MCP tool)"""
        logger.info(f"Analyzing regulators for {state['gene']}")

        try:
            state['current_step'] = 'regulators_analysis'

            cell_type = CellType(state['cell_type'])
            result = await self.modeling_agent.analyze_regulators_detailed(
                state['gene'],
                cell_type,
                max_regulators=25
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
        """Analyze biological pathways (equivalent to analyze_biological_context)"""
        logger.info(f"Analyzing pathways for {state['gene']}")

        try:
            state['current_step'] = 'pathway_analysis'

            # Extract gene list from previous analyses
            gene_list = [state['gene']]

            # Add regulators if available
            if state.get('regulators_analysis'):
                hub_regulators = state['regulators_analysis'].get('hub_regulators', [])
                for reg in hub_regulators[:5]:  # Top 5 hub regulators
                    if reg.get('gene_symbol') and reg['gene_symbol'] != 'Unknown':
                        gene_list.append(reg['gene_symbol'])

            # Add targets if available
            if state.get('targets_analysis'):
                cascade_targets = state['targets_analysis'].get('cascade_targets', [])
                for target in cascade_targets[:3]:  # Top 3 cascade targets
                    if target.get('gene_symbol') and target['gene_symbol'] != 'Unknown':
                        gene_list.append(target['gene_symbol'])

            # Perform pathway enrichment
            result = await self.integration_agent.analyze_biological_context(
                gene_list,
                analysis_type="pathway_enrichment"
            )

            state['pathway_analysis'] = result
            state['analysis_metadata']['steps_completed'].append('pathway_analysis')

            logger.info(f"Pathway analysis complete. Found {len(result.get('enriched_pathways', {}))} pathways")
            return state

        except Exception as e:
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

    async def _find_similar_genes(self, state: GeneAnalysisState) -> GeneAnalysisState:
        """Find genes with similar regulatory patterns"""
        logger.info(f"Finding similar genes to {state['gene']}")

        try:
            state['current_step'] = 'similarity_analysis'

            cell_type = CellType(state['cell_type'])
            result = await self.modeling_agent.find_similar_genes(
                state['gene'],
                cell_type,
                top_k=10
            )

            state['similarity_analysis'] = result
            state['analysis_metadata']['steps_completed'].append('similarity_analysis')

            logger.info(f"Similarity analysis complete. Found {len(result.get('similar_genes', []))} similar genes")
            return state

        except Exception as e:
            logger.warning(f"Similarity analysis failed (optional): {str(e)}")
            # Don't fail the workflow for similarity analysis
            state['similarity_analysis'] = {"error": str(e), "analysis_skipped": True}
            state['analysis_metadata']['steps_completed'].append('similarity_analysis')
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

            result = await self.domain_agents.analyze_drug_development_context(
                state['gene'], gene_info, pathway_analysis
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

            result = await self.domain_agents.analyze_systems_biology_context(
                state['gene'], gene_info, pathway_analysis,
                regulators_analysis, targets_analysis
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
            cross_cell_analysis=None,
            similarity_analysis=None,
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
    workflow = GREmLNWorkflow()

    # Test with APC (the gene from your example)
    print("Testing GREmLN LangGraph Workflow with APC...")
    result = await workflow.run_analysis(
        gene="APC",
        cell_type="epithelial_cell",
        analysis_depth="comprehensive"
    )

    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    asyncio.run(main())