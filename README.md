# GeneExpressionBreastCancerRelapse
Project Status: In Progress
## Project Overview
This project investigates whether gene expression can help distinguish breast cancer patients who experienced relapse from those who did not. Using the publicaly available GSE2034 dataset, the project applies statistical analysis and machine learning techniques to develop and evaluate a predictive model on it's ability to predict relapse versus relapse free based on gene expression values. <br>

This project focuses on data preprocessing, exploratory data analysis, machine learning feature selection, and predictive model evaluation. It also considers the challenge of applying the predictive model to independent patient cohorts.
Repository for analysis and predictive model development for evaluating gene expressions and their relationships to breast cancer relapse.

### Research Question
Can a machine learning model trained on gene expression data from the GSE2034 breast cancer cohort predict relapse status, and how well does this predictive performance generalize to an independent dataset?

### Background
Gene expression profiling allows researchers to examine patterns of gene activity across thousands of measurements. These patterns may contain information associated with clinical outcomes and prognosis. However, gene expression data is high dimensional, with the dataset in this project containing over 22k probes(features). This poses several analytical challenges, many more features than samples, class imbalance, correlated predictors, and risk of overfitting. This project explores approaches for addressing these challenges while evaluating the potential predictive value of gene expression data.

### Data and Provenance
The original dataset comes from a research paper from The Lancet entitled ‘Gene-expression profiles to predict distant metastasis of lymph-node-negative primary breast cancer’. The data used in this project comes from the Kaggle Breast Cancer Gene Expression Dataset (GSE2034). It is CC0 1.0 Universal, meaning the data is public domain. <br>
This project also includes the clinical parameter data obtained from the NCBI Geo Omnibus platform, which is also public domain. 
The clinical variables extracted from the clinical parameter data are the Relapse status and the ER status.

## To get started:
1. Download the data files provided.
2. Run the EDA_Report file to quality control the data, and join the necessary columns into one comprehensive dataframe. 
