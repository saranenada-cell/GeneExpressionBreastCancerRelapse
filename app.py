from pathlib import Path

import pandas as pd
import streamlit as st
import altair as alt
import numpy as np

st.set_page_config(
    page_title="GSE2034 Research Explorer",
    page_icon="🧬",
    layout="wide"
)

# Load the verified dataset
DATA_PATH = Path(__file__).parent / "data" / "dashboard_dataset.csv"
df = pd.read_csv(DATA_PATH)

# Dashboard title
st.title("GSE2034 Research Explorer")

st.write(
    "An interactive exploration of breast cancer "
    "gene expression and relapse prediction."
)

# Cohort overview
st.header("Cohort Overview")

col1, col2, col3 = st.columns(3)

col1.metric("Total Samples", len(df))
col2.metric("Relapse Positive", int(df["Relapse"].sum()))
col3.metric("Gene Expression Probes", df.shape[1] - 3)

st.subheader("Clinical Data Preview")
st.dataframe(
    df[["Sample_Accession", "Relapse", "ER_status"]].head(10),
    use_container_width=True
)

'''Relapse distribution Chart for visualization'''
st.subheader("Relapse Status Distribution")

relapse_counts = (
    df["Relapse"]
    .map({0: "Relapse-free", 1: "Relapse"})
    .value_counts()
    .reindex(["Relapse-free", "Relapse"], fill_value=0)
)

chart_col, summary_col = st.columns([2, 1])

with chart_col:
    st.bar_chart(relapse_counts)

with summary_col:
    st.markdown("**Cohort breakdown**")

    total = len(df)

    for label in ["Relapse-free", "Relapse"]:
        count = int(relapse_counts[label])
        percentage = count / total * 100 if total else 0

        st.metric(label, f"{percentage:.1f}%")
        st.caption(f"{count:,} of {total:,} samples")


'''ER Status Distribution'''
st.subheader("Estrogen Receptor (ER) Status")

er_counts = (
    df["ER_status"]
    .value_counts()
    .reindex(["ER-", "ER+"], fill_value=0)
)

chart_col, summary_col = st.columns([2, 1])

with chart_col:
    st.bar_chart(er_counts)

with summary_col:
    st.markdown("**Cohort breakdown**")

    total = len(df)

    for label in ["ER-", "ER+"]:
        count = int(er_counts[label])
        percentage = count / total * 100 if total else 0

        st.metric(label, f"{percentage:.1f}%")
        st.caption(f"{count:,} of {total:,} samples")

'''Proportional Comparison of ER vs Relapse'''

st.subheader("Relapse Status by ER Status")

# Count samples in each ER and relapse group
relapse_by_er = (
    df.assign(
        Relapse_Label=df["Relapse"].map(
            {0: "No relapse", 1: "Relapse"}
        )
    )
    .groupby(["ER_status", "Relapse_Label"])
    .size()
    .reset_index(name="Sample Count")
)

# Grouped bar chart
chart = alt.Chart(relapse_by_er).mark_bar().encode(
    x=alt.X("ER_status:N", title="ER Status"),
    xOffset=alt.XOffset("Relapse_Label:N"),
    y=alt.Y("Sample Count:Q", title="Number of Samples"),
    color=alt.Color("Relapse_Label:N", title="Relapse Status"),
    tooltip=["ER_status", "Relapse_Label", "Sample Count"]
)

# Calculate relapse rates within each ER group
er_rates = df.groupby("ER_status")["Relapse"].mean() * 100

er_negative_rate = er_rates.get("ER-", 0)
er_positive_rate = er_rates.get("ER+", 0)
difference = er_positive_rate - er_negative_rate

# Display chart and proportional comparison side by side
chart_col, summary_col = st.columns([2, 1])

with chart_col:
    st.altair_chart(chart, use_container_width=True)

with summary_col:
    st.markdown("**Relapse rate within each ER group**")
    st.metric("ER-negative", f"{er_negative_rate:.1f}%")
    st.metric("ER-positive", f"{er_positive_rate:.1f}%")
    st.metric(
        "Difference (ER+ minus ER−)",
        f"{difference:+.1f} percentage points"
    )

st.caption(
    "Descriptive cohort comparison only; the difference does not "
    "establish causation or statistical significance."
)


'''GENE EXPRESSION ANALYSIS SECTION'''
'''Mann-Whitney U test for every expression probe and adjusting the p-value using benjamini-hochberg FDR'''


from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# --- Differential expression analysis ---

# Identify expression probe columns
metadata_columns = [
    "Sample_Accession",
    "Relapse",
    "ER_status"
]

probe_columns = [
    col for col in df.columns
    if col not in metadata_columns
]

# Separate the two relapse groups
relapse_free = df[df["Relapse"] == 0]
relapse_positive = df[df["Relapse"] == 1]

results = []

for probe in probe_columns:
    group_0 = relapse_free[probe].dropna()
    group_1 = relapse_positive[probe].dropna()

    if len(group_0) == 0 or len(group_1) == 0:
        continue

    test = mannwhitneyu(
        group_0,
        group_1,
        alternative="two-sided"
    )

    results.append({
        "Probe ID": probe,
        "Relapse-free Median": group_0.median(),
        "Relapse Median": group_1.median(),
        "Median Difference": group_1.median() - group_0.median(),
        "P-value": test.pvalue
    })

# Create results table
results_df = pd.DataFrame(results)

# Adjust p-values for multiple testing
results_df["FDR-adjusted P-value"] = multipletests(
    results_df["P-value"],
    method="fdr_bh"
)[1]

results_df["Significant (FDR < 0.05)"] = (
    results_df["FDR-adjusted P-value"] < 0.05
)

results_df = results_df.sort_values(
    "FDR-adjusted P-value"
).reset_index(drop=True)

#Display the results

results_df = results_df.sort_values(
    "FDR-adjusted P-value"
).reset_index(drop=True)

# --- Gene Expression Explorer ---

'''Gene Annotation Preview Table'''

# --- Prepare gene annotations ---
annotation_path = (
    Path(__file__).parent
    / "data"
    / "GSE2034_family_platform_GPL96_annotation.csv"
)


# --- Prepare gene annotations ---
annotation_df = pd.read_csv(annotation_path)

gene_annotations = annotation_df[
    [
        "ID",
        "Gene Symbol",
        "Gene Title",
        "Gene Ontology Biological Process"
    ]
].copy()

# Clean text fields
for col in ["Gene Symbol", "Gene Title"]:
    gene_annotations[col] = (
        gene_annotations[col]
        .fillna("")
        .astype(str)
        .str.strip()
    )

# Remove probes without a usable gene symbol
gene_annotations = gene_annotations[
    (gene_annotations["Gene Symbol"] != "")
    & (gene_annotations["Gene Symbol"] != "---")
    & (gene_annotations["Gene Symbol"].str.lower() != "nan")
]

st.subheader("Gene Annotation Preview")
st.write(f"Annotated probe records: {len(gene_annotations):,}")
st.dataframe(gene_annotations.head(10), width="stretch")

'''Searchable Gene Selector'''


# --- Gene Explorer ---
st.divider()
st.header("Gene Explorer")

# Match annotated probe IDs to expression columns
available_probes = [
    probe for probe in gene_annotations["ID"].unique()
    if probe in df.columns
]

explorer_genes = gene_annotations[
    gene_annotations["ID"].isin(available_probes)
].copy()

# Create readable labels for the dropdown
explorer_genes["Display Name"] = (
    explorer_genes["Gene Symbol"]
    + " — "
    + explorer_genes["Gene Title"]
    + " ["
    + explorer_genes["ID"]
    + "]"
)

explorer_genes = explorer_genes.drop_duplicates(
    subset=["ID"]
).sort_values("Display Name")


# Filter gene explorer by statistical significance
explorer_results = explorer_genes.merge(
    results_df[
        ["Probe ID", "FDR-adjusted P-value",
         "Significant (FDR < 0.05)"]
    ],
    left_on="ID",
    right_on="Probe ID",
    how="left"
)

significance_filter = st.radio(
    "Which probes would you like to explore?",
    ["All annotated probes", "Significant probes only"],
    horizontal=True
)

if significance_filter == "Significant probes only":
    filtered_genes = explorer_results[
        explorer_results["Significant (FDR < 0.05)"] == True
    ]
else:
    filtered_genes = explorer_results

if filtered_genes.empty:
    st.warning("No probes match this filter.")
    st.stop()

selected_display = st.selectbox(
    "Search for a gene",
    options=filtered_genes["Display Name"].tolist()
)


selected_gene = filtered_genes[
    filtered_genes["Display Name"] == selected_display
].iloc[0]

selected_probe = selected_gene["ID"]

st.markdown(f"**Gene symbol:** {selected_gene['Gene Symbol']}")
st.markdown(f"**Gene title:** {selected_gene['Gene Title']}")

st.caption(f"Expression probe: {selected_probe}")

biological_process = selected_gene["Gene Ontology Biological Process"]

if biological_process:
    st.markdown("**Known biological processes (GO annotation)**")
    st.write(biological_process)
else:
    st.info("No biological process annotation is available for this probe.")



# --- Expression comparison ---
st.subheader("Expression by Relapse Status")

plot_df = df[["Relapse", selected_probe]].copy()

plot_df["Relapse Status"] = plot_df["Relapse"].map(
    {0: "Relapse-free", 1: "Relapse"}
)

chart = alt.Chart(plot_df).mark_boxplot().encode(
    x=alt.X(
        "Relapse Status:N",
        title="Relapse Status",
        sort=["Relapse-free", "Relapse"]
    ),
    y=alt.Y(
        f"{selected_probe}:Q",
        title="Expression Value"
    ),
    color=alt.Color(
        "Relapse Status:N",
        legend=None
    ),
    tooltip=["Relapse Status", selected_probe]
)

st.altair_chart(chart, width="stretch")


# --- Statistical results for the selected gene probe ---
probe_result = results_df[
    results_df["Probe ID"] == selected_probe
]

if not probe_result.empty:
    result = probe_result.iloc[0]

    col1, col2, col3, col4 = st.columns(4)

    col1.metric(
        "Relapse-free median",
        f"{result['Relapse-free Median']:.3f}"
    )
    col2.metric(
        "Relapse median",
        f"{result['Relapse Median']:.3f}"
    )
    col3.metric(
        "Mann–Whitney U p-value",
        f"{result['P-value']:.3g}"
    )
    col4.metric(
        "FDR-adjusted p-value",
        f"{result['FDR-adjusted P-value']:.3g}"
    )

    if result["Significant (FDR < 0.05)"]:
        st.success("Significant after FDR correction (adjusted p < 0.05).")
    else:
        st.info("Not significant after FDR correction (adjusted p ≥ 0.05).")
else:
    st.warning("Statistical results were not found for this probe.")


st.caption(
    "This chart compares expression values for the selected probe. "
    "Differences are exploratory and do not establish statistical "
    "significance or predictive performance."
)




# --- Differential Expression Results ---
st.header("Differential Expression Results")

st.write(
    "Genes are ranked by their FDR-adjusted p-values "
    "from the Mann–Whitney U test."
)

st.write(f"Probes tested: {len(results_df):,}")

significant_results = results_df[
    results_df["FDR-adjusted P-value"] < 0.05
]

st.write(
    f"Probes significant after FDR correction: "
    f"{len(significant_results):,}"
)


# Add gene annotations to the statistical results
#st.write("Results columns:", results_df.columns.tolist())
#st.write("Annotation columns:", gene_annotations.columns.tolist())

annotated_results = results_df.merge(
    gene_annotations[["ID", "Gene Symbol", "Gene Title"]],
    left_on="Probe ID",
    right_on="ID",
    how="left"
)


# Put gene information first
annotated_results = annotated_results[
    [
        "Probe ID",
        "Gene Symbol",
        "Gene Title",
        "Relapse-free Median",
        "Relapse Median",
        "Median Difference",
        "P-value",
        "FDR-adjusted P-value",
        "Significant (FDR < 0.05)"
    ]
]

# Filter option
show_significant_only = st.checkbox(
    "Show only statistically significant probes (FDR < 0.05)",
    value=False
)

if show_significant_only:
    display_results = annotated_results[
        annotated_results["Significant (FDR < 0.05)"]
    ]
else:
    display_results = annotated_results

st.write(f"Showing {len(display_results):,} probes")

# Highlight significant results in light green
def highlight_significance(row):
    if row["Significant (FDR < 0.05)"]:
        return [
            "background-color: #d9ead3"
            if col == "Significant (FDR < 0.05)"
            else ""
            for col in row.index
        ]
    return [""] * len(row)

styled_results = display_results.style.apply(
    highlight_significance,
    axis=1
)

st.dataframe(
    styled_results,
    width="stretch",
    hide_index=True
)





# --- Top 20 Most Significant Probes ---
st.divider()
st.header("Top 20 Most Significant Probes")

st.write(
    "Probes are ranked by FDR-adjusted p-value. "
    "The median difference indicates the direction and size "
    "of the difference between the two groups' median expression."
)

top20 = results_df.nsmallest(
    20,
    "FDR-adjusted P-value"
).copy()

# Add gene symbols
top20 = top20.merge(
    gene_annotations[["ID", "Gene Symbol"]],
    left_on="Probe ID",
    right_on="ID",
    how="left"
)

top20["Gene Symbol"] = top20["Gene Symbol"].fillna("Unknown")

top20["Chart Label"] = (
    top20["Gene Symbol"]
    + " ("
    + top20["Probe ID"]
    + ")"
)

# Direction of median expression difference
top20["Expression Direction"] = np.where(
    top20["Median Difference"] > 0,
    "Higher in relapse",
    np.where(
        top20["Median Difference"] < 0,
        "Lower in relapse",
        "Equal medians"
    )
)

# Plot median difference, with zero as the reference point
top20_chart = alt.Chart(top20).mark_bar().encode(
    x=alt.X(
        "Median Difference:Q",
        title="Median Difference (relapse − relapse-free)"
    ),
    y=alt.Y(
        "Chart Label:N",
        sort=alt.EncodingSortField(
            field="FDR-adjusted P-value",
            order="ascending"
        ),
        title="Gene symbol (probe ID)"
    ),
    color=alt.Color(
        "Expression Direction:N",
        scale=alt.Scale(
            domain=[
                "Higher in relapse",
                "Lower in relapse",
                "Equal medians"
            ],
            range=["#d95f02", "#1b9e77", "#999999"]
        ),
        title="Direction"
    ),
    tooltip=[
        "Probe ID",
        "Gene Symbol",
        "Median Difference",
        "Expression Direction",
        "P-value",
        "FDR-adjusted P-value"
    ]
).properties(height=500)

st.altair_chart(top20_chart, width="stretch")
