import streamlit as st
import pandas as pd

import plotly.express as px
import numpy as np


# Page title
st.title("Gene Expression Patterns")

st.write(
    "Explore gene expression differences between "
    "patients who experienced relapse and those who remained relapse-free."
)

# Load differential expression results
de_results = pd.read_csv(
    "data/differential_expression_results.csv"
)

# ----------------------------------
# Page Development Notice
# ----------------------------------

with st.container(border=True):
    st.markdown("### 🚧 This page is a work in progress!")

    st.write(
        "Currently, you can explore the volcano plot to see "
        "which genes are associated with breast cancer relapse "
        "and view the differential expression results in the table below."
    )

    st.markdown("**Coming soon:**")
    st.write(
        "- Explore gene expression by other patient characteristics, "
        "including ER status."
    )
    st.write(
        "- Investigate associations with time to relapse and survival."
    )
    st.write(
        "- Explore the biological functions and pathways of significant genes."
    )
# ----------------------------------
# Adding gene-expression summary
# ----------------------------------

# Identify significant probes
significant = de_results[de_results["Significant"] == True]

# Count significant probes by expression direction
higher_expression = (
    significant["Mean_Difference"] > 0
).sum()

lower_expression = (
    significant["Mean_Difference"] < 0
).sum()

# Overview metrics
st.subheader("Differential Expression Overview")

col1, col2, col3, col4 = st.columns(4)

col1.metric("Probes Analyzed", f"{len(de_results):,}")
col2.metric("Significant Probes", f"{len(significant):,}")
col3.metric("Higher in Relapse", f"{higher_expression:,}")
col4.metric("Lower in Relapse", f"{lower_expression:,}")


# ----------------------------------
# Volcano Plot
# ----------------------------------

# Volcano plot
st.subheader("Differential Expression Volcano Plot")

volcano_df = de_results.copy()

# Calculate negative log10 adjusted p-value
volcano_df["-log10(FDR)"] = -np.log10(
    volcano_df["Adjusted_P_Value"].clip(lower=1e-300)
)

# Assign significance categories
volcano_df["Expression"] = np.where(
    volcano_df["Significant"],
    np.where(
        volcano_df["Mean_Difference"] > 0,
        "Higher in relapse",
        "Lower in relapse"
    ),
    "Not significant"
)

fig = px.scatter(
    volcano_df,
    x="Mean_Difference",
    y="-log10(FDR)",
    color="Expression",
    hover_data=[
        "Probe_ID",
        "Gene Symbol",
        "Gene Title",
        "Adjusted_P_Value"
    ],
    color_discrete_map={
        "Higher in relapse": "#D96C6C",
        "Lower in relapse": "#397CA5",
        "Not significant": "#C8CDD2"
    },
    labels={
        "Mean_Difference": "Mean expression difference (log2)",
        "-log10(FDR)": "-log10 (adjusted p-value)"
    }
)

fig.add_hline(
    y=-np.log10(0.05),
    line_dash="dash",
    line_color="black",
    annotation_text="FDR = 0.05"
)

fig.add_vline(x=0, line_dash="dash", line_color="black")

fig.update_layout(
    height=600,
    legend_title_text="Expression difference",
    xaxis_title="Mean expression difference (Relapse − No relapse)",
    yaxis_title="-log10 (FDR-adjusted p-value)"
)

st.plotly_chart(fig, use_container_width=True)

# ----------------------------------
# Significant Genes Table
# ----------------------------------

# Gene expression results table
st.subheader("Gene Expression Results")

# Filter option
significant_only = st.checkbox(
    "View significant probes only",
    value=True
)

# Apply filter
if significant_only:
    display_df = de_results[
        de_results["Significant"] == True
    ].copy()
else:
    display_df = de_results.copy()

# Select and organize columns
display_columns = [
    "Probe_ID",
    "Gene Symbol",
    "Gene Title",
    "Mean_Difference",
    "P_Value",
    "Adjusted_P_Value"
]

display_df = display_df[display_columns]

display_df = display_df.rename(columns={
    "Probe_ID": "Probe ID",
    "Gene Symbol": "Gene Symbol",
    "Gene Title": "Gene Title",
    "Mean_Difference": "Mean Difference (log2)",
    "P_Value": "Raw P-value",
    "Adjusted_P_Value": "Adjusted P-value"
})

display_df = display_df.sort_values("Adjusted P-value")

st.caption(f"Showing {len(display_df):,} probes")

st.dataframe(
    display_df,
    use_container_width=True,
    hide_index=True,
    column_config={
        "Mean Difference (log2)": st.column_config.NumberColumn(format="%.3f"),
        "Raw P-value": st.column_config.NumberColumn(format="%.2e"),
        "Adjusted P-value": st.column_config.NumberColumn(format="%.4f")
    }
)


# ----------------------------------
# Gene Explorer
# ----------------------------------

dashboard_dataset = pd.read_csv(
    "data/dashboard_dataset.csv"
)

st.subheader("Individual Gene Explorer")

st.write(
    "Select a gene to explore its expression levels "
    "in patients who experienced relapse versus those who did not."
)

# Filter option
significant_genes_only = st.checkbox(
    "Show significant genes only",
    value=False
)

# Keep genes with an available gene symbol
available_genes = de_results.dropna(
    subset=["Gene Symbol"]
).copy()

# Filter to significant probes if selected
if significant_genes_only:
    available_genes = available_genes[
        available_genes["Significant"] == True
    ]

# Create sorted gene options
gene_options = sorted(
    available_genes["Gene Symbol"].unique().tolist()
)

if not gene_options:
    st.warning("No genes are available with the selected filter.")

else:
    # Gene selector
    selected_gene = st.selectbox(
        "Select a gene",
        options=gene_options
    )

    # Find probes associated with selected gene
    gene_probes = available_genes[
        available_genes["Gene Symbol"] == selected_gene
    ].sort_values("Adjusted_P_Value")

    # Use probe with lowest adjusted p-value
    selected_probe = gene_probes.iloc[0]["Probe_ID"]

    # Gene information
    gene_title = gene_probes.iloc[0]["Gene Title"]

    st.markdown(f"**Gene:** {selected_gene}")
    st.caption(gene_title)
    st.caption(f"Representative probe: {selected_probe}")

    # Statistical results
    adjusted_p = gene_probes.iloc[0]["Adjusted_P_Value"]
    mean_difference = gene_probes.iloc[0]["Mean_Difference"]
    is_significant = gene_probes.iloc[0]["Significant"]

    st.divider()
    st.subheader("Statistical Results")

    if is_significant:
        st.success(
            f"Statistically significant "
            f"(FDR-adjusted p-value = {adjusted_p:.4e})"
        )
    else:
        st.info(
            f"Not statistically significant "
            f"(FDR-adjusted p-value = {adjusted_p:.4e})"
        )

    # Direction of expression difference
    if mean_difference > 0:
        direction = "Higher expression in the relapse group"
    elif mean_difference < 0:
        direction = "Higher expression in the no-relapse group"
    else:
        direction = "Equal mean expression between groups"

    st.markdown(f"**Direction:** {direction}")

    st.write(
        f"**Mean expression difference:** {mean_difference:.3f} "
        "(Relapse − No relapse)"
    )

    st.caption(
        "Statistical significance is defined as a Benjamini–Hochberg "
        "adjusted p-value below 0.05. Direction describes the difference "
        "in group means, not a causal relationship."
    )

    # ----------------------------------
    # Relapse Expression Plot
    # ----------------------------------

    st.divider()
    st.subheader("Gene Expression by Relapse Status")

    expression_df = dashboard_dataset[
        ["Relapse", selected_probe]
    ].copy()

    expression_df["Relapse Group"] = expression_df["Relapse"].map({
        0: "No relapse",
        1: "Relapse"
    })

    fig = px.box(
        expression_df,
        x="Relapse Group",
        y=selected_probe,
        color="Relapse Group",
        color_discrete_map={
            "No relapse": "#397CA5",
            "Relapse": "#D96C6C"
        },
        points="all",
        labels={
            selected_probe: "Gene expression (log2)",
            "Relapse Group": "Patient group"
        }
    )

    fig.update_layout(
        height=450,
        showlegend=False,
        xaxis_title="",
        yaxis_title="Gene expression (log2)"
    )

    st.plotly_chart(
        fig,
        use_container_width=True,
        key="relapse_expression_chart"
    )


