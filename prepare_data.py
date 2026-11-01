import pandas as pd

# Load clinical data
clinical_data = pd.read_csv("data/clinical_data.txt", sep="\t")

# Rename columns
clinical_parameters_df = clinical_data.rename(
    columns={
        "GEO asscession number": "Sample_Accession",
        "relapse (1=True)": "Relapse",
        "ER Status": "ER_status"
    }
)

print(clinical_parameters_df.head())
print("\nMissing values by column:")
print(clinical_parameters_df.isnull().sum())


# Load gene expression data
gene_seq_df = pd.read_csv("data/GSE2034_family_expression_matrix.csv")

print(gene_seq_df.head())
print("\nMissing values by column:")
print(gene_seq_df.isnull().sum())


# Load ABS_CALL detection data
PAM_data_df = pd.read_csv("data/GSE2034_family_matrix_ABS_CALL.csv")

# Calculate Present detection rate for each probe
pam_value_row = PAM_data_df.iloc[:, 1:]

present_detection_rate = (pam_value_row == "P").mean(axis=1)

present_detection_summary = pd.DataFrame({
    "ID_REF": PAM_data_df.iloc[:, 0],
    "Detection_Rate": present_detection_rate.values
})

print(present_detection_summary.head())

# Retain probes with a Present detection rate of at least 50%
filtered_present_detection_df = present_detection_summary.loc[
    present_detection_summary["Detection_Rate"] >= 0.5,
    ["ID_REF", "Detection_Rate"]
]

print(filtered_present_detection_df.head())

# Filter expression matrix using retained probe IDs
probes_to_keep = filtered_present_detection_df["ID_REF"]

filtered_expression_df = gene_seq_df[
    gene_seq_df["ID_REF"].isin(probes_to_keep)
]

print(filtered_expression_df.head())
print("Filtered expression matrix shape:", filtered_expression_df.shape)

# Set probe IDs as index and transpose expression matrix
expression_by_sample = filtered_expression_df.set_index("ID_REF").T

# Make sample accession numbers a regular column
expression_by_sample.index.name = "Sample_Accession"
expression_by_sample = expression_by_sample.reset_index()

print("Expression matrix shape:", expression_by_sample.shape)
print(expression_by_sample.head())


# Select relevant clinical variables
clinical_info = clinical_parameters_df[
    ["Sample_Accession", "Relapse", "ER_status"]
]

# Merge expression data with clinical information
final_expression_df = expression_by_sample.merge(
    clinical_info,
    on="Sample_Accession",
    how="inner",
    validate="one_to_one"
)

print("Final dataset shape:", final_expression_df.shape)
print(final_expression_df[["Sample_Accession", "Relapse", "ER_status"]].head())

# Save processed dataset
final_expression_df.to_csv(
    "data/dashboard_dataset.csv",
    index=False
)

print("Processed dataset saved successfully!")