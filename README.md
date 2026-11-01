# Gene Expression Analysis and Breast Cancer Relapse Prediction

**[View Live Interactive Dashboard](https://gse2034-breast-cancer-relapse.streamlit.app/)**

## Project Overview
This project explores gene expression patterns associated with breast cancer relapse using the GSE2034 microarray dataset. The goal is to investigate differences in gene expression between patients with and without relapse and provide an interactive tool for exploring these patterns.

The project includes data preprocessing, exploratory data analysis, statistical testing, and an interactive Streamlit dashboard. The dashboard allows users to examine gene annotations, compare expression levels across relapse groups, and explore statistically significant differences. <br>
**Project Status:** In Progress

### Dataset Description
Source: Kaggle dataset containing a copy of the GSE2034 breast cancer gene-expression dataset.
- Dataset: GSE2034 - Breast Cancer relapse-free survival
- Primary Repository: NCBI Gene Expression Omnibus (GEO)
- GEO accession: GSE2034
- Publication status: Publicly available as of February 23, 2005
- Original Study: Wang et al., Gene-expression profiles to predict distant metastasis of lymph-node-negative primary breast cancer
- Original Expression Data: 22,283 microarray probe measurements
- Analysis Cohort: 286 primary breast tumour samples after merging expression and clinical data
- Clinical Information: Relapse status, time to relapse or last follow-up, lymph-node status, ER status, and other recurrence-related variables
- License: CC0 1.0 Universal (public domain dedication)

### Background
Gene expression profiling allows researchers to examine patterns of gene activity across thousands of measurements. These patterns may contain information associated with clinical outcomes and prognosis. However, gene expression data is high dimensional, with the dataset in this project containing over 22k probes(features). This poses several analytical challenges, many more features than samples, class imbalance, correlated predictors, and risk of overfitting. This project explores approaches for addressing these challenges while evaluating the potential predictive value of gene expression data.

### Data Source and Licensing
This project uses the GSE2034 breast cancer gene-expression dataset, obtained from the Kaggle dataset Breast Cancer Gene Expression Dataset (GSE2034), which is marked CC0 1.0 Public Domain on Kaggle.
The original study and dataset are available through NCBI Gene Expression Omnibus (GSE2034). Please consult the original GEO record and publication for source information and appropriate attribution.
This project is an independent educational analysis and is not a clinical diagnostic tool.


## To get started:
1. Download the data files provided. (All files that begin with GSE_2034)
3. Run the EDA_Report file to review the Exploratory Analysis, quality control the probe data, and join the necessary columns into one comprehensive dataframe. 
