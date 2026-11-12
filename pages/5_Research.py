
import streamlit as st
import pandas as pd
import plotly.express as px
import numpy as np



# --------------------------------------------------
# Reading required datasets
# --------------------------------------------------
@st.cache_data
def load_de_results():
    er_positive = pd.read_csv(
        "data/ER_positive_differential_expression.csv"
    )
    er_negative = pd.read_csv(
        "data/ER_negative_differential_expression.csv"
    )

    return er_positive, er_negative

interaction_results = pd.read_csv("data/interaction_results.csv")

gene_results = pd.read_csv("data/gene_level_results_clean.csv")

# --------------------------------------------------
# Required plots
# --------------------------------------------------
def create_volcano_plot(df, title):
    fig = px.scatter(
        df,
        x="Mean_Difference",
        y="Neg_Log10_P",
        color="Significant",
        hover_data=[
            "Probe_ID",
            "Mean_Difference",
            "P_Value",
            "Adjusted_P_Value"
        ],
        color_discrete_map={
            True: "#D95F59",
            False: "#B8C4CC"
        },
        title=title,
        labels={
            "Mean_Difference": "Mean Difference (Relapse − No Relapse)",
            "Neg_Log10_P": "−log₁₀(P-value)",
            "Significant": "FDR Significant"
        }
    )

    fig.update_layout(
        template="plotly_white",
        height=420,
        legend_title_text="FDR Significant",
        margin=dict(l=20, r=20, t=55, b=20)
    )

    return fig

def create_interaction_plot(df):
    fig = px.scatter(
        df,
        x="Interaction_Coefficient",
        y=-np.log10(df["P_Value"].clip(lower=1e-300)),
        color="Significant",
        hover_data=[
            "Probe_ID",
            "Interaction_Coefficient",
            "P_Value",
            "Adjusted_P_Value"
        ],
        color_discrete_map={
            True: "#D95F59",
            False: "#B8C4CC"
        },
        labels={
            "Interaction_Coefficient": "Interaction Coefficient",
            "y": "−log₁₀(P-value)",
            "Significant": "FDR Significant"
        }
    )

    fig.update_layout(
        template="plotly_white",
        height=420,
        margin=dict(l=20, r=20, t=45, b=20),
        showlegend=True
    )

    return fig

def create_chromosome_plot(df):

    plot_df = df.copy()

    plot_df["Position_Mb"] = (
        (plot_df["start"] + plot_df["end"]) / 2
    ) / 1_000_000

    fig = px.scatter(
        plot_df,
        x="Position_Mb",
        y="Chromosome",
        color="Expression_Category",
        hover_name="Gene_Symbol",
        hover_data={
            "Entrez_ID": True,
            "Mean_Difference": ":.3f",
            "Adjusted_P_Value": ":.3e",
            "Position_Mb": ":.2f",
            "Chromosome": True
        },
        color_discrete_map={
            "Higher in relapse": "#D95F59",
            "Higher in no relapse": "#3182BD",
            "Conflicting": "#9467BD",
            "Not significant": "#D9DEE3"
        },
        labels={
            "Position_Mb": "Genomic Position (Mb)",
            "Chromosome": "Chromosome",
            "Expression_Category": "Expression Category"
        }
    )

    fig.update_layout(
        template="plotly_white",
        height=750,
        margin=dict(l=20, r=20, t=30, b=20),
        legend_title_text="Expression Category",
        yaxis=dict(
            categoryorder="array",
            categoryarray=[
                str(i) for i in range(1, 23)
            ] + ["X"]
        )
    )

    return fig


st.title("Research Questions & Findings")

st.markdown("""
This page presents the research questions guiding the investigation
of gene expression patterns and breast cancer relapse in the GSE2034 cohort.

Each section describes the research question, analytical approach,
key findings, and interpretation.
""")

st.divider()

# --------------------------------------------------
# 1. ER Subgroup Analysis
# --------------------------------------------------

st.header("1. ER Subgroup Analysis")

st.subheader(
    "Does gene expression associated with relapse differ by ER status?"
)
er_positive, er_negative = load_de_results()

col1, col2 = st.columns(2)

with col1:
    st.markdown("**ER-positive patients (n=209)**")
    fig_positive = create_volcano_plot(
        er_positive,
        "ER-positive"
    )
    st.plotly_chart(fig_positive, use_container_width=True)

with col2:
    st.markdown("**ER-negative patients (n=77)**")
    fig_negative = create_volcano_plot(
        er_negative,
        "ER-negative"
    )
    st.plotly_chart(fig_negative, use_container_width=True)
st.markdown("""
**Analytical approach**

Differential expression analysis was performed separately in
ER-positive and ER-negative patients using Welch's t-tests
and Benjamini-Hochberg FDR correction.

**Key findings**

- ER-positive: 134 FDR-significant probes.
- ER-negative: 0 FDR-significant probes.
- ER-positive patients: 209.
- ER-negative patients: 77.

**Interpretation**

Relapse-associated differential expression was detected in the
ER-positive subgroup, but not in the ER-negative subgroup.
The smaller ER-negative sample size may have limited statistical power.
""")

st.divider()

# --------------------------------------------------
# 2. ER Interaction Analysis
# --------------------------------------------------

st.header("2. ER Interaction Analysis")

st.subheader(
    "Does ER status modify the relationship between gene expression and relapse?"
)

col1, col2 = st.columns([1.5, 1])

with col1:
    fig_interaction = create_interaction_plot(interaction_results)
    st.plotly_chart(fig_interaction, use_container_width=True)

with col2:
    st.markdown("""
    **Analytical approach**

    Gene-level interaction models were used to investigate whether
    the association between expression and relapse differed by ER status.
    
    **Findings**

    - Probes tested: 11,130
    - FDR-significant interactions: 0

    **Interpretation**

    No statistically significant evidence was found that ER status modifies the relationship between gene expression and relapse.

    The smaller ER-negative subgroup may limit the power to detect interaction effects.
    """)


st.divider()



# --------------------------------------------------
# 3. Chromosomal Distribution
# --------------------------------------------------

st.header("3. Chromosomal Distribution")

st.subheader(
    "Where are relapse-associated genes located across the genome?"
)
col1, col2 = st.columns([1.7, 1])

with col1:
    fig_chromosome = create_chromosome_plot(gene_results)
    st.plotly_chart(
        fig_chromosome,
        use_container_width=True
    )

with col2:
    st.markdown("""
    **Analytical approach**

    Significant differential-expression probes were mapped to
    GRCh38 genomic coordinates using Entrez Gene IDs.
    Multiple probes mapping to the same gene were consolidated
    for visualization.
    
    **Findings**

    - Genes mapped: 7,397
    - Higher expression in relapse: 141
    - Higher expression in no relapse: 52
    - Conflicting probe directions: 1

    **Interpretation**

    Relapse-associated genes are distributed across multiple chromosomes.

    The visualization allows individual genes to be explored by genomic position, expression direction, and statistical significance.
    """)


st.divider()

# --------------------------------------------------
# 4. ER and Predictive Performance
# --------------------------------------------------

st.header("4. ER Status and Relapse Prediction")

st.subheader(
    "Does including ER status improve our ability to predict relapse?"
)

st.info("This investigation is currently in progress.")

st.markdown("""
**Planned approach**

Compare machine-learning models using:

- Gene expression alone.
- Gene expression combined with ER status.
- ER status alone.

Model performance will be evaluated using ROC-AUC,
balanced accuracy, sensitivity, and specificity.
""")
