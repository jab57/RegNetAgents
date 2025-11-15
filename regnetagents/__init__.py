"""
RegNetAgents - Core Package

Core utilities for gene regulatory network analysis.
"""

from .complete_gene_service import CompleteGeneService, GeneInfo, get_complete_gene_service
from .gene_id_mapper import GeneIDMapper

__all__ = [
    'CompleteGeneService',
    'GeneInfo', 
    'get_complete_gene_service',
    'GeneIDMapper'
]

__version__ = '1.0.0'
