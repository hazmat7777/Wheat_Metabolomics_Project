# Repository README

This repository contains R scripts and results files for a machine learning analysis predicting take-all disease status in wheat using soil microbiome and metabolomic data.

## Analysis Scripts

### Data Processing & Preparation
- **`GT_datawrang.R`** - Data wrangling and preparation of take-all (GT) status data
- **`metabolomics_data_wrang2.R`** - Processing and cleaning of metabolomics LCMS data  
- **`smp_output_singletons_taxaassigned.R`** - Processing and filtering 16S metabarcoding data

### Feature Selection & Stability Analysis
- **`BORUTA_MICROB_AND_METAB.R`** - Combined Boruta feature selection for both microbiome and metabolomic datasets
- **`boruta_microb_fs_juststability.R`** - Boruta feature selection focused on microbial data stability testing
- **`boruta_selection_metab5_stability_cv_acc_se.R`** - Metabolite feature selection with cross-validation and stability assessment
- **`FS_metab_results.R`** - Analysis of feature-selected metabolomic model results

### Machine Learning Models
- **`metab_boruta_linearmodels.R`** - Linear modeling approaches on Boruta-selected metabolites
- **`GLMs.R`** - Generalized linear models for disease prediction
- **`save_metab_rf_results.R`** - Random forest modeling on metabolomic data
- **`metab_pca_then_rf8_pcarf_cv_auc.R`** - PCA-transformed metabolomic data with random forest modeling
- **`microb_rf_fs_visualisn.R`** - Random forest modeling and visualization for microbial data

### Cross-Validation & Model Evaluation
- **`boruta_microbe2_cv_auc.R`** - Cross-validation and AUC assessment for microbial models
- **`metab_mega_results_df_maker6_CV.R`** - Comprehensive cross-validation results compilation for metabolomic models
- **`microb_mega_results_df_maker3_auc.R`** - AUC-based evaluation and results compilation for microbial models
- **`boruta_linear_models_comparison.csv`** - Comparative performance metrics across linear modeling approaches on microbiome/metabolome data

### Correlation & Integration Analysis
- **`metab_165_correl.R`** - Correlation analysis between metabolites and microbial taxa
- **`diversity_microb_metab_vs_gt4_nofs.R`** - Diversity analysis comparing microbiome and metabolome with take-all status

### Figure Generation
- **`FIGURE_1.R`** - Generates Figure 1
- **`FIGURE_2.R`** - Generates Figure 2 
- **`FIGURE_3.R`** - Generates Figure 3  
- **`FIGURE_4.R`** - Generates Figure 4 

