# GeneExpressionBreastCancerRelapse
# Project Overview
This project investigates whether gene expression can help distinguish breast cancer patients who experienced relapse from those who did not. Using the publicaly available GSE2034 dataset, the project applies statistical analysis and machine learning techniques to develop and evaluate a predictive model on it's ability to predict relapse versus relapse free based on gene expression values. <br>

This project focuses on data preprocessing, exploratory data analysis, machine learning feature selection, and predictive model evaluation. It also considers the challenge of applying the predictive model to independent patient cohorts.
Repository for analysis and predictive model development for evaluating gene expressions and their relationships to breast cancer relapse.

This project looks into which gene-expression features are most associated with relapse in breast cancer patients, how well a gene-expression profile can differentiate relapse vs relapse free patients, and looks into if separating based on ER status affects the models ability to predict relapse.

The projects works off of the Kaggle breast cancer data set, as well as the NCBI GEO database, containing the clinical relapse parameters.

To get started:
Download the data files provided.
Run the EDA_Report file to quality control the data, and join the necessary columns into one comprehensive dataframe. 
