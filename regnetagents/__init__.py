"""
RegNetAgents - Core Package

Core utilities for gene regulatory network analysis.
"""

from .complete_gene_service import CompleteGeneService, GeneInfo, get_complete_gene_service
from .gene_id_mapper import GeneIDMapper
from .tcga_registry import TCGA_NETWORK_REGISTRY, TCGA_CANCER_TYPES
from .network_loader import load_network

__all__ = [
    'CompleteGeneService',
    'GeneInfo',
    'get_complete_gene_service',
    'GeneIDMapper',
    'TCGA_NETWORK_REGISTRY',
    'TCGA_CANCER_TYPES',
    'load_network',
]

__version__ = '1.2.8'
