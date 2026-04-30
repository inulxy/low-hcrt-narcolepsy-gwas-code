# Low HCRT Narcolepsy GWAS Code

This repository contains analysis notebooks and R scripts used for the low hypocretin narcolepsy GWAS research.

## Contents

- `00_Low_HCRT_GWAS_PreQC.ipynb`: Pre-QC workflow
- `01_Low_HCRT_GWAS_QC.ipynb`: Sample and variant QC
- `02_Low_HCRT_GWAS_Imputation.ipynb`: Imputation workflow
- `03_Low_HCRT_GWAS_PCA.ipynb`: PCA and ancestry-related analysis
- `04_Low_HCRT_GWAS_SAIGE.ipynb`: SAIGE association testing
- `05_Low_HCRT_GWAS_Meta_Analysis.ipynb`: Meta-analysis workflow
- `Deming_regression_OR_comparison.R`: Odds ratio comparison and Deming regression analysis

## Input data for Deming regression
The file `data/compare.xlsx` contains summary-level association statistics used for the odds-ratio comparison and Deming regression analysis. 


## Requirements

Main tools used include:
- R
- Jupyter Notebook
- PLINK
- SAIGE
- Tractor
- METAL
