#!/usr/bin/env python3
"""
Data Validation Script for New Cell Types

This script validates single-cell RNA-seq datasets before processing them
through the ARACNe pipeline for RegNetAgents network generation.
"""

import os
import sys
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Optional
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Cell type specific requirements
CELL_TYPE_REQUIREMENTS = {
    'hepatocytes': {
        'min_cells': 3000,
        'min_genes': 15000,
        'max_mito_pct': 25,
        'required_markers': ['ALB', 'CYP3A4', 'CYP2E1', 'G6PC'],
        'tissue_keywords': ['liver', 'hepatic', 'hepatocyte'],
        'description': 'Liver parenchymal cells for drug metabolism research'
    },
    'cardiomyocytes': {
        'min_cells': 2000,
        'min_genes': 12000,
        'max_mito_pct': 30,
        'required_markers': ['TNNT2', 'MYH6', 'MYH7', 'ACTC1'],
        'tissue_keywords': ['heart', 'cardiac', 'cardiomyocyte'],
        'description': 'Heart muscle cells for cardiotoxicity studies'
    },
    'neurons': {
        'min_cells': 5000,
        'min_genes': 18000,
        'max_mito_pct': 20,
        'required_markers': ['MAP2', 'RBFOX3', 'SYP', 'SNAP25'],
        'tissue_keywords': ['brain', 'neural', 'neuron', 'cortex'],
        'description': 'Brain and nervous system cells for neurological research'
    },
    'fibroblasts': {
        'min_cells': 3000,
        'min_genes': 14000,
        'max_mito_pct': 25,
        'required_markers': ['COL1A1', 'VIM', 'FN1', 'THY1'],
        'tissue_keywords': ['skin', 'lung', 'fibroblast', 'connective'],
        'description': 'Connective tissue cells for cancer stroma research'
    },
    'endothelial_cells': {
        'min_cells': 2500,
        'min_genes': 13000,
        'max_mito_pct': 25,
        'required_markers': ['PECAM1', 'VWF', 'CDH5', 'KDR'],
        'tissue_keywords': ['endothelial', 'vascular', 'vessel'],
        'description': 'Blood vessel lining cells for angiogenesis research'
    }
}

class CellTypeDataValidator:
    """Validator for cell type specific single-cell RNA-seq data"""

    def __init__(self, cell_type: str):
        if cell_type not in CELL_TYPE_REQUIREMENTS:
            raise ValueError(f"Unsupported cell type: {cell_type}")

        self.cell_type = cell_type
        self.requirements = CELL_TYPE_REQUIREMENTS[cell_type]
        self.validation_results = {}

    def validate_dataset(self, data_path: str) -> Dict:
        """
        Validate a single-cell dataset for the specified cell type

        Args:
            data_path: Path to H5AD file

        Returns:
            Dictionary containing validation results
        """
        logger.info(f"Validating {self.cell_type} dataset: {data_path}")

        try:
            # Load data
            adata = self._load_data(data_path)
            if adata is None:
                return {'valid': False, 'error': 'Failed to load data'}

            # Run validation checks
            results = {
                'valid': True,
                'cell_type': self.cell_type,
                'data_path': data_path,
                'checks': {}
            }

            # Basic data structure checks
            results['checks']['data_structure'] = self._check_data_structure(adata)

            # Cell count validation
            results['checks']['cell_count'] = self._check_cell_count(adata)

            # Gene count validation
            results['checks']['gene_count'] = self._check_gene_count(adata)

            # Cell type annotation validation
            results['checks']['cell_annotations'] = self._check_cell_annotations(adata)

            # Marker gene validation
            results['checks']['marker_genes'] = self._check_marker_genes(adata)

            # Quality metrics validation
            results['checks']['quality_metrics'] = self._check_quality_metrics(adata)

            # Tissue context validation
            results['checks']['tissue_context'] = self._check_tissue_context(adata)

            # Overall validation status
            results['valid'] = all(check.get('passed', False) for check in results['checks'].values())

            # Summary statistics
            results['summary'] = self._generate_summary(adata, results)

            return results

        except Exception as e:
            logger.error(f"Validation failed: {e}")
            return {'valid': False, 'error': str(e)}

    def _load_data(self, data_path: str):
        """Load single-cell data from file"""
        try:
            import scanpy as sc
            logger.info(f"Loading data from {data_path}")
            adata = sc.read_h5ad(data_path)
            logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
            return adata
        except ImportError:
            logger.error("scanpy not available. Install with: pip install scanpy")
            return None
        except Exception as e:
            logger.error(f"Failed to load data: {e}")
            return None

    def _check_data_structure(self, adata) -> Dict:
        """Check basic data structure requirements"""
        checks = {
            'has_obs': hasattr(adata, 'obs') and len(adata.obs) > 0,
            'has_var': hasattr(adata, 'var') and len(adata.var) > 0,
            'has_X': hasattr(adata, 'X') and adata.X is not None,
            'valid_shape': adata.X.shape == (adata.n_obs, adata.n_vars)
        }

        return {
            'passed': all(checks.values()),
            'details': checks,
            'message': 'Data structure validation' + (' passed' if all(checks.values()) else ' failed')
        }

    def _check_cell_count(self, adata) -> Dict:
        """Check if dataset has sufficient cells"""
        cell_count = adata.n_obs
        min_required = self.requirements['min_cells']

        return {
            'passed': cell_count >= min_required,
            'count': cell_count,
            'required': min_required,
            'message': f'Cell count: {cell_count} (required: {min_required})'
        }

    def _check_gene_count(self, adata) -> Dict:
        """Check if dataset has sufficient genes"""
        gene_count = adata.n_vars
        min_required = self.requirements['min_genes']

        return {
            'passed': gene_count >= min_required,
            'count': gene_count,
            'required': min_required,
            'message': f'Gene count: {gene_count} (required: {min_required})'
        }

    def _check_cell_annotations(self, adata) -> Dict:
        """Check cell type annotations"""
        # Look for common cell type annotation columns
        annotation_columns = ['cell_type', 'celltype', 'Cell_type', 'cluster', 'leiden', 'seurat_clusters']
        found_columns = [col for col in annotation_columns if col in adata.obs.columns]

        if not found_columns:
            return {
                'passed': False,
                'message': 'No cell type annotation columns found',
                'available_columns': list(adata.obs.columns)
            }

        # Check if target cell type is present
        target_cell_type = self.cell_type
        cell_type_found = False
        matching_cells = 0

        for col in found_columns:
            unique_values = adata.obs[col].astype(str).str.lower().unique()

            # Look for exact matches or partial matches
            for value in unique_values:
                if (target_cell_type.lower() in value or
                    any(keyword in value for keyword in self.requirements['tissue_keywords'])):
                    cell_type_found = True
                    matching_cells = sum(adata.obs[col].astype(str).str.lower() == value)
                    break

            if cell_type_found:
                break

        return {
            'passed': cell_type_found and matching_cells >= self.requirements['min_cells'] / 2,
            'cell_type_found': cell_type_found,
            'matching_cells': matching_cells,
            'annotation_columns': found_columns,
            'message': f'Cell type annotations: {matching_cells} {target_cell_type} cells found' if cell_type_found else 'Target cell type not found'
        }

    def _check_marker_genes(self, adata) -> Dict:
        """Check for presence of cell type specific marker genes"""
        required_markers = self.requirements['required_markers']
        gene_names = adata.var_names.str.upper()

        found_markers = []
        missing_markers = []

        for marker in required_markers:
            if marker.upper() in gene_names.values:
                found_markers.append(marker)
            else:
                missing_markers.append(marker)

        marker_score = len(found_markers) / len(required_markers)

        return {
            'passed': marker_score >= 0.5,  # At least 50% of markers should be present
            'score': marker_score,
            'found_markers': found_markers,
            'missing_markers': missing_markers,
            'message': f'Marker genes: {len(found_markers)}/{len(required_markers)} found'
        }

    def _check_quality_metrics(self, adata) -> Dict:
        """Check data quality metrics"""
        try:
            import scanpy as sc

            # Calculate QC metrics if not present
            if 'n_genes_by_counts' not in adata.obs.columns:
                sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)

            # Calculate mitochondrial gene percentage if not present
            if 'pct_counts_mt' not in adata.obs.columns:
                adata.var['mt'] = adata.var_names.str.startswith('MT-')
                sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

            # Quality checks
            mean_genes_per_cell = adata.obs['n_genes_by_counts'].mean()
            mean_counts_per_cell = adata.obs['total_counts'].mean()
            mean_mito_pct = adata.obs.get('pct_counts_mt', pd.Series([0])).mean()

            quality_checks = {
                'sufficient_genes_per_cell': mean_genes_per_cell >= 500,
                'sufficient_counts_per_cell': mean_counts_per_cell >= 1000,
                'acceptable_mito_pct': mean_mito_pct <= self.requirements['max_mito_pct']
            }

            return {
                'passed': all(quality_checks.values()),
                'mean_genes_per_cell': mean_genes_per_cell,
                'mean_counts_per_cell': mean_counts_per_cell,
                'mean_mito_pct': mean_mito_pct,
                'checks': quality_checks,
                'message': f'Quality metrics: {mean_genes_per_cell:.0f} genes/cell, {mean_mito_pct:.1f}% mito'
            }

        except Exception as e:
            return {
                'passed': False,
                'error': str(e),
                'message': 'Failed to calculate quality metrics'
            }

    def _check_tissue_context(self, adata) -> Dict:
        """Check if dataset matches expected tissue context"""
        tissue_keywords = self.requirements['tissue_keywords']

        # Check metadata for tissue context
        metadata_text = ""
        if hasattr(adata, 'uns') and isinstance(adata.uns, dict):
            metadata_text += str(adata.uns).lower()

        # Check observation columns
        for col in adata.obs.columns:
            if 'tissue' in col.lower() or 'organ' in col.lower():
                metadata_text += " " + " ".join(adata.obs[col].astype(str).unique()).lower()

        # Check variable names (gene names might indicate tissue)
        # Add observation metadata
        metadata_text += " " + " ".join(adata.obs.columns).lower()

        # Look for tissue keywords
        found_keywords = [kw for kw in tissue_keywords if kw in metadata_text]

        return {
            'passed': len(found_keywords) > 0,
            'found_keywords': found_keywords,
            'message': f'Tissue context: {found_keywords if found_keywords else "No matching keywords found"}'
        }

    def _generate_summary(self, adata, results: Dict) -> Dict:
        """Generate summary statistics and recommendations"""
        summary = {
            'dataset_size': f"{adata.n_obs} cells × {adata.n_vars} genes",
            'validation_score': sum(1 for check in results['checks'].values() if check.get('passed', False)) / len(results['checks']),
            'ready_for_processing': results['valid'],
            'cell_type': self.cell_type,
            'description': self.requirements['description']
        }

        # Generate recommendations
        recommendations = []
        for check_name, check_result in results['checks'].items():
            if not check_result.get('passed', False):
                if check_name == 'cell_count':
                    recommendations.append(f"Increase cell count to at least {self.requirements['min_cells']}")
                elif check_name == 'marker_genes':
                    recommendations.append("Consider using a different dataset with better cell type markers")
                elif check_name == 'quality_metrics':
                    recommendations.append("Apply quality filtering to improve data quality")
                elif check_name == 'cell_annotations':
                    recommendations.append("Verify cell type annotations or use computational annotation")

        summary['recommendations'] = recommendations

        return summary

def validate_all_cell_types(base_dir: str = "models/networks") -> Dict:
    """Validate all new cell type datasets"""

    results = {}

    for cell_type in CELL_TYPE_REQUIREMENTS.keys():
        logger.info(f"\n{'='*50}")
        logger.info(f"Validating {cell_type}")
        logger.info(f"{'='*50}")

        # Look for data files
        cell_type_dir = os.path.join(base_dir, cell_type, "raw_data")
        data_files = []

        if os.path.exists(cell_type_dir):
            for file in os.listdir(cell_type_dir):
                if file.endswith(('.h5ad', '.h5')):
                    data_files.append(os.path.join(cell_type_dir, file))

        if not data_files:
            results[cell_type] = {
                'valid': False,
                'error': f'No data files found in {cell_type_dir}',
                'status': 'NO_DATA'
            }
            logger.warning(f"No data files found for {cell_type}")
            continue

        # Validate each data file
        validator = CellTypeDataValidator(cell_type)

        for data_file in data_files:
            logger.info(f"Validating file: {data_file}")
            file_results = validator.validate_dataset(data_file)
            results[f"{cell_type}_{os.path.basename(data_file)}"] = file_results

            # Print validation summary
            if file_results['valid']:
                logger.info(f"✓ {cell_type} validation PASSED")
            else:
                logger.warning(f"✗ {cell_type} validation FAILED")

            if 'summary' in file_results:
                summary = file_results['summary']
                logger.info(f"  Dataset: {summary['dataset_size']}")
                logger.info(f"  Score: {summary['validation_score']:.2f}")

                if summary['recommendations']:
                    logger.info("  Recommendations:")
                    for rec in summary['recommendations']:
                        logger.info(f"    - {rec}")

    return results

def main():
    """Main validation function"""
    import argparse

    parser = argparse.ArgumentParser(description="Validate cell type datasets for RegNetAgents processing")
    parser.add_argument('--cell-type', choices=list(CELL_TYPE_REQUIREMENTS.keys()),
                       help='Specific cell type to validate')
    parser.add_argument('--data-file', help='Specific data file to validate')
    parser.add_argument('--base-dir', default='models/networks',
                       help='Base directory containing cell type data')
    parser.add_argument('--report', action='store_true',
                       help='Generate detailed validation report')

    args = parser.parse_args()

    if args.cell_type and args.data_file:
        # Validate specific file
        validator = CellTypeDataValidator(args.cell_type)
        results = validator.validate_dataset(args.data_file)

        print(f"\nValidation Results for {args.cell_type}:")
        print(f"File: {args.data_file}")
        print(f"Valid: {'✓' if results['valid'] else '✗'}")

        if 'summary' in results:
            summary = results['summary']
            print(f"Dataset: {summary['dataset_size']}")
            print(f"Validation Score: {summary['validation_score']:.2f}")

    else:
        # Validate all cell types
        results = validate_all_cell_types(args.base_dir)

        print(f"\n{'='*60}")
        print("VALIDATION SUMMARY")
        print(f"{'='*60}")

        valid_count = 0
        total_count = 0

        for name, result in results.items():
            if 'error' not in result:
                total_count += 1
                if result.get('valid', False):
                    valid_count += 1
                    status = "✓ READY"
                else:
                    status = "✗ NEEDS WORK"

                print(f"{name:<30} {status}")

        print(f"\nOverall: {valid_count}/{total_count} datasets ready for processing")

        if args.report:
            # Save detailed report
            import json
            with open('validation_report.json', 'w') as f:
                json.dump(results, f, indent=2, default=str)
            print(f"\nDetailed report saved to: validation_report.json")

if __name__ == "__main__":
    main()