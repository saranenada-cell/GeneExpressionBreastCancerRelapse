from pathlib import Path
import pandas as pd

data_folder = Path("data")

clinical_path = data_folder / "clinical_data.txt"

clinical_df = pd.read_csv(clinical_path, sep="\t")
print("Clinical data shape:", clinical_df.shape)
print("\nColumn names:")
print(clinical_df.columns.tolist())

print("\nFirst five rows:")
print(clinical_df.head())

'''Checking the expression data'''

clean_path = data_folder / "final_expression_matrix.csv"

clean_df = pd.read_csv(clean_path)

print("\nCleaned expression data shape:", clean_df.shape)
print("\nFirst five column names:")
print(clean_df.columns[:5].tolist())

print("\nFirst five rows:")
print(clean_df.iloc[:5, :5])

'''Checking the shape of the data because there was an unnamed column'''

print("\nCleaned matrix shape:", clean_df.shape)
print("\nUnique sample accessions:", clean_df["Sample_Accession"].nunique())
print("\nDuplicate sample accessions:", clean_df["Sample_Accession"].duplicated().sum())

print("\nClinical sample accessions:", clinical_df["GEO asscession number"].nunique())

'''Checking to see whether the sample IDs match'''

expression_ids = set(clean_df["Sample_Accession"])
clinical_ids = set(clinical_df["GEO asscession number"])

print("Samples in expression but not clinical:", len(expression_ids - clinical_ids))
print("Samples in clinical but not expression:", len(clinical_ids - expression_ids))

'''Creating a merged dataset'''

# Remove the saved DataFrame index column
clean_df = clean_df.drop(columns=["Unnamed: 0"])

# Rename clinical columns for easier use
clinical_df = clinical_df.rename(columns={
    "GEO asscession number": "Sample_Accession",
    "relapse (1=True)": "Relapse",
    "ER Status": "ER_status"
})

# Keep only the clinical variables needed for this project
clinical_subset = clinical_df[
    ["Sample_Accession", "Relapse", "ER_status"]
]

# Merge expression and clinical data
merged_df = clean_df.merge(
    clinical_subset,
    on="Sample_Accession",
    how="inner",
    validate="one_to_one"
)

print("Merged dataset shape:", merged_df.shape)
print("\nClinical columns:")
#print(merged_df[["Sample_Accession", "Relapse", "ER_status"]].head())
#print("\nMissing clinical values:")
#print(merged_df[["Relapse", "ER_status"]].isna().sum())

print("Columns containing relapse or ER:")
print([
    col for col in merged_df.columns
    if "relapse" in col.lower() or "er_status" in col.lower() or "er status" in col.lower()
])

'''For some reason relapse and status was not findable, turns out it was expored with _x or _y'''
for col in ["Relapse_x", "Relapse_y", "ER_status_x", "ER_status_y"]:
    print(f"\n{col}:")
    print(merged_df[col].value_counts(dropna=False))


'''Comparing the duplicated columns. Both results were true! Nice'''

print(
    "Relapse columns identical:",
    merged_df["Relapse_x"].equals(merged_df["Relapse_y"])
)

print(
    "ER status columns identical:",
    merged_df["ER_status_x"].equals(merged_df["ER_status_y"])
)

'''Cleaning up and saving the verified dataset'''

# Keep one copy of each clinical label
merged_df = merged_df.drop(
    columns=["Relapse_y", "ER_status_y"]
)

# Rename the remaining columns
merged_df = merged_df.rename(columns={
    "Relapse_x": "Relapse",
    "ER_status_x": "ER_status"
})

# Save the verified dataset
output_path = data_folder / "dashboard_dataset.csv"
merged_df.to_csv(output_path, index=False)

print("Saved dataset:", output_path)
print("Final shape:", merged_df.shape)
print("Missing values in clinical labels:")
print(merged_df[["Relapse", "ER_status"]].isna().sum())