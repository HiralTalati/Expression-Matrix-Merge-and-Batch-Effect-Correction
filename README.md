**Cross-Cohort Transcriptomic Integration & Protected Batch Correction Pipeline**

This repository contains an advanced R-based data science pipeline designed to merge independent gene expression matrices, dynamically recover clinical metadata using NCBI GEO, and perform biological-signal-protected batch effect correction using ComBat.

The pipeline is explicitly optimized for transcriptomic analysis of neurological disorders across heterogeneous studies (specifically combining GSE63060, GSE63061, and GSE140829 formats) to output a harmonized, multi-cohort dataset ready for downstream machine learning or differential expression workflows.

📌 Features

• Automated Metadata Harvesting: Programmatically queries the NCBI Gene Expression Omnibus (GEO) via GEOquery to reconstruct and map clinical metadata factors (status:ch1 / diagnosis:ch1) on the fly, eliminating the need for manual tracking files.

• Phenotype Standardization: Harmonizes categorical raw diagnostic strings into three strict clinical classifications: Controls (CTL), Mild Cognitive Impairment (MCI), and Alzheimer's Disease (AD).

• Unbiased Multi-Cohort Merge: Aggregates variable-width expression data frames using a relational reduce(inner_join) configuration mapped against stable NCBI Entrez IDs (EntrezID).

• Variance Auditing via Linear Modeling: Computes Principal Component Analysis (PCA) and calculates the exact coefficient of determination (R²) from a linear model (lm(PC1 ~ Batch)) to quantify the statistical influence of batch effects before and after correction.

• Biological-Signal-Protected ComBat: Implements an empirical Bayes model matrix using sva::ComBat to shield critical clinical phenotypes (group_clean) from distortion while completely removing technical noise from sample batches.

🛠️ Tech Stack & Dependencies

Ensure you have the following CRAN and Bioconductor packages installed before initiating execution:

# CRAN Packages

install.packages(c("openxlsx", "dplyr", "purrr", "tibble", "ggplot2"))

# Bioconductor Packages

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "sva"))

📂 Data Structure & Pipeline Inputs

The entry phase assumes that individual array preprocessing steps have run successfully and deposited their final assets in designated directories within your active root workspace.

Expected Folder Schema:

├── GSE63060_Final_Results_20260314_171513/

│   ├── [ID]_06_EXPRESSION_MATRIX.xlsx

│   └── [ID]_04_ANNOTATED_ALL.xlsx

├── GSE63061_Final_Results_20260314_171611/

│   ├── [ID]_06_EXPRESSION_MATRIX.xlsx

│   └── [ID]_04_ANNOTATED_ALL.xlsx

└── GSE140829_Final_Results_20260314_180745/
    
    ├── [ID]_06_EXPRESSION_MATRIX.xlsx
    
    └── [ID]_04_ANNOTATED_ALL.xlsx

Note: If your local results subfolder timestamps vary, simply modify the named folders character vector inside Step 8 of your script.

🚀 Execution & Processing Workflow

The pipeline executes sequentially in two distinct automated phases.

Phase 1: Cross-Mapping & Unification (Step 8)

Run the first section of the script to extract expression matrices and map phenotypes. The pipeline automatically connects to the GEO network API to download and parse corresponding sample pData maps.

Rscript step8_master_merge.R

Mid-Stage Outputs Generated:

• 07_MASTER_MERGED_EXPRESSION.xlsx: An integrated matrix indexing all common Entrez IDs (rows) against individual Sample Accession IDs (columns).

• 08_BATCH_METADATA.xlsx: A companion lookup spreadsheet mapping each SampleID to its source dataset (Batch) and verified clinical cohort (group_clean).

Phase 2: Protected Batch Correction & Audit (Step 9)

Once Phase 1 files are available in your root directory, launch the ComBat correction script:


Rscript step9_batch_correction.R

This phase takes the raw merged matrix, evaluates initial technical variance, builds a covariate protection matrix, and applies empirical Bayes corrections.

📊 Evaluation Assets & Output Summary

Upon successful completion, a new production folder titled After_Batch_Correction_Results/ will be generated with the following diagnostic and data logs:

└── After_Batch_Correction_Results/
    
    ├── 09_PCA_Pre_Correction.png        # Scatter plot highlighting severe batch clustering
    
    ├── 10_PCA_Post_Correction.png       # Normalized PCA plot showcasing unified data integration
    
    ├── 11_FINAL_BATCH_CORRECTED_MATRIX.xlsx # The target batch-adjusted matrix for downstream modeling
    
    └── 12_Batch_Audit_Comparison.xlsx   # Audit summary file detailing R² reduction logs

Expected Audit Results (Sample Logs):

The pipeline outputs the statistical variance changes in a format ready for academic inclusion:

• Pre-Correction: Batch noise typically explains a massive majority of the primary Principal Component variance (R² ≈ 80-95%).

• Post-Correction: The remaining batch effect influence is successfully dropped to nominal levels (R² ≈ 0.00%), while true diagnostic cluster separation (CTL vs. MCI vs. AD) is verified via shape geometries on the post-correction plot.

📄 License

This analysis framework is open-source and distributed under the MIT License.
