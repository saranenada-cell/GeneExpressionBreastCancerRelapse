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

### Research Objective
The primary objective of this project is to investigate gene-expression patterns associated with breast cancer relapse and evaluate the potential of machine-learning approaches to predict relapse using gene-expression data. <br>
The project has two main components:
1. Exploratory Data Analysis: This includes examining the dataset and it's distribution of gene expression, exploring clinical characteristics, and identifying gene-expression probes associated with differences in relapse. This project will also apply statistical testing and multiple testing correction to evaluate differential gene expression. An interactive is developed to support exploration of gene-expression patterns and associated clinical outcomes.
2. Predictive Modelling and External Validation: This includes applying feature selection techniques to identify informative gene-expression features, developing a classifier, evaluating the models performance, and asesssing generalizability using an external and independent dataset.


### Data Source and Licensing
This project uses the GSE2034 breast cancer gene-expression dataset, obtained from the Kaggle dataset Breast Cancer Gene Expression Dataset (GSE2034), which is marked CC0 1.0 Public Domain on Kaggle.
The original study and dataset are available through NCBI Gene Expression Omnibus (GSE2034). Please consult the original GEO record and publication for source information and appropriate attribution.
This project is an independent educational analysis and is not a clinical diagnostic tool.


## To get started:
1. Download the data files provided. (All files that begin with GSE_2034)
3. Run the EDA_Report file to review the Exploratory Analysis, quality control the probe data, and join the necessary columns into one comprehensive dataframe. 
