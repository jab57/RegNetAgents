# Extended Cancer Gene Analysis (Archived)

**Date**: December 6, 2024
**Status**: Exploratory analysis - not included in manuscript

## What This Is

An exploratory analysis testing the RegNetAgents framework on 10 additional cancer genes across 4 cancer types:
- Breast Cancer: BRCA1, BRCA2, ERBB2, ESR1
- Lung Cancer: EGFR, ALK, RET
- Prostate Cancer: AR, PTEN
- Pancreatic Cancer: CDKN2A

## Results Summary

**Success**: 3/10 genes (30%) yielded meaningful therapeutic target data
- BRCA1: 23 regulators, top: ZNF334
- BRCA2: 20 regulators, top: HMGB2
- EGFR: 37 regulators, top: HMGB2

**Limited data**: 7/10 genes had insufficient regulatory network coverage in epithelial cells

## Key Insight

The results revealed that **pre-screening genes for regulatory network coverage** is essential before applying the framework. Random selection of cancer genes yields low success rates because:
1. ARACNe networks are context-specific (cell-type dependent)
2. Not all cancer genes have rich transcriptional regulatory networks
3. Some genes are regulated post-transcriptionally

## Why Archived

This analysis was not included in the bioRxiv manuscript because:
- The 5-gene colorectal cancer validation is already sufficient
- Low success rate (30%) doesn't strengthen the validation
- Reveals network limitations rather than demonstrating broad applicability
- Pre-screening approach (Option A) would be needed for stronger validation

## Files

- `run_extended_cancer_analysis.py` - Analysis script
- `extended_validation_assessment.md` - Detailed assessment and recommendations
- `results/` - Individual gene reports and perturbation analyses

## Future Use

This analysis could be valuable for:
- Responding to reviewer questions about broader applicability
- Developing a "gene pre-screening" feature for the framework
- Future journal submission with proper gene selection methodology
