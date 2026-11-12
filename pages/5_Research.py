
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


def create_er_correlation_plot(er_positive, er_negative):

    # Select the 134 significant ER-positive probes
    positive_sig = er_positive[
        er_positive["Significant"].astype(str).str.lower() == "true"
    ][["Probe_ID", "Mean_Difference"]].copy()

    # Match those same probes in the ER-negative results
    negative_results = er_negative[
        ["Probe_ID", "Mean_Difference"]
    ].copy()

    # Merge using Probe_ID
    comparison = positive_sig.merge(
        negative_results,
        on="Probe_ID",
        suffixes=("_ER_Positive", "_ER_Negative")
    )

    # Calculate Pearson correlation
    correlation = comparison[
        "Mean_Difference_ER_Positive"
    ].corr(
        comparison["Mean_Difference_ER_Negative"]
    )

    # Determine common axis limits for the reference line
    min_value = comparison[
        ["Mean_Difference_ER_Positive", "Mean_Difference_ER_Negative"]
    ].min().min()

    max_value = comparison[
        ["Mean_Difference_ER_Positive", "Mean_Difference_ER_Negative"]
    ].max().max()

    # Create interactive scatter plot
    fig = px.scatter(
        comparison,
        x="Mean_Difference_ER_Positive",
        y="Mean_Difference_ER_Negative",
        hover_data=["Probe_ID"],
        title=f"ER+ vs ER− Expression Differences (r = {correlation:.3f})",
        labels={
            "Mean_Difference_ER_Positive": "Mean Difference (ER-positive)",
            "Mean_Difference_ER_Negative": "Mean Difference (ER-negative)"
        }
    )

    # Add diagonal reference line (identical effect sizes)
    fig.add_shape(
        type="line",
        x0=min_value,
        y0=min_value,
        x1=max_value,
        y1=max_value,
        line=dict(color="red", dash="dash", width=2)
    )

    fig.update_layout(
        template="plotly_white",
        height=450,
        margin=dict(l=20, r=20, t=60, b=20)
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

st.markdown("#### Consistency of expression patterns across ER subgroups")

col1, col2 = st.columns([1.3, 1])

with col1:
    fig_correlation = create_er_correlation_plot(
        er_positive,
        er_negative
    )

    st.plotly_chart(
        fig_correlation,
        use_container_width=True
    )

with col2:
    st.markdown("""
    
    **Approach**
    There are no significant probes in the ER- subgroup, however it would be interesting to see how the signficant probes in the ER+ group behave in the ER- group.

    **Findings**

    - Probes compared: 134
    - Pearson correlation: 0.763
    - Same direction in ER−:
        - 87 higher in relapse
        - 47 higher in no relapse
    - Opposite directions: 0

    **Interpretation**

    The positive correlation suggests that relapse-associated expression patterns are broadly consistent across ER subgroups. AKA when that probe is high in expression in the ER+ group, it also seems to be high in the ER- group.

    However, this does not establish statistical significance in ER-negative patients or prove identical biological effects. 
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

st.header("4. ER Status and Predictive Performance")

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

st.subheader(
    "How can we improve the Models Performance"
)
st.markdown("""
**Potential Approaches**

Compare machine-learning models using:

- ER status as predictor.
- RFE for feature selection.
- Elastic Net Logistic Regression

Model performance will be evaluated using ROC-AUC,
balanced accuracy, sensitivity, and specificity.
""")
