from pathlib import Path

import pandas as pd
import streamlit as st
import altair as alt
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

st.markdown("""
<style>
div[data-testid="stVerticalBlockBorderWrapper"] {
    border-radius: 12px;
    border: 1px solid rgba(128, 128, 128, 0.18);
    box-shadow: 0px 2px 8px rgba(0, 0, 0, 0.07);
}
</style>
""", unsafe_allow_html=True)

# ----------------------------------
# Load dataset
# ----------------------------------

DATA_PATH = (
    Path(__file__).parent.parent
    / "data"
    / "dashboard_dataset.csv"
)

@st.cache_data
def load_data(path):
    return pd.read_csv(path)

st.cache_data.clear()


df = load_data(DATA_PATH)

# ----------------------------------
# Page introduction
# ----------------------------------

st.title("Breast Cancer Relapse Prediction Project")

st.write(
    "This is an overview of the GSE2034 breast cancer cohort, " \
    "published in ‘Gene-expression profiles to predict distant metastasis " \
    "of lymph-node-negative primary breast cancer’ paper in The Lancet, " \
    "providing a snapshot into the data presented in this report. " \
    "On this page you will find descriptive information on the total samples, " \
    "relapses cases, ER status of the patients, and the distribution of relapse " \
    "showcased in this study."
)

st.divider()


# ----------------------------------
# 1. Cohort Summary
# ----------------------------------

with st.container(border=True):

    st.subheader("Cohort Summary")

    # Calculate cohort counts
    total_samples = len(df)
    probe_count = df.shape[1] - 3

    relapse_counts = (
        df["Relapse"]
        .map({0: "Relapse-free", 1: "Relapse"})
        .value_counts()
        .reindex(["Relapse-free", "Relapse"], fill_value=0)
    )

    er_counts = (
        df["ER_status"]
        .value_counts()
        .reindex(["ER-", "ER+"], fill_value=0)
    )


    # Create chart data
    relapse_chart_data = pd.DataFrame({
        "Category": relapse_counts.index,
        "Count": relapse_counts.values
    })

    er_chart_data = pd.DataFrame({
        "Category": er_counts.index,
        "Count": er_counts.values
    })


    # Calculate percentages for tooltips
    relapse_chart_data["Percentage"] = (
        relapse_chart_data["Count"] / total_samples * 100
    )

    er_chart_data["Percentage"] = (
        er_chart_data["Count"] / total_samples * 100
    )


    # Layout
    col1, col2, col3 = st.columns(3)

    with col1:
        st.metric("Total Samples", f"{total_samples:,}")
        st.caption(f"{probe_count:,} gene expression probes")

    with col2:
        st.markdown("**Relapse Status**")

        relapse_chart = (
            alt.Chart(relapse_chart_data)
            .mark_arc(outerRadius=85)
            .encode(
                theta=alt.Theta("Count:Q"),
                color=alt.Color(
                    "Category:N",
                    title="Relapse Status",
                    scale=alt.Scale(
                        domain=["Relapse-free", "Relapse"],
                        range=["#00436C", "#91C0CB"]
                    )
                ),
                tooltip=[
                    alt.Tooltip("Category:N", title="Status"),
                    alt.Tooltip("Count:Q", title="Samples"),
                    alt.Tooltip(
                        "Percentage:Q",
                        title="Percentage",
                        format=".1f"
                    )
                ]
            )
            .properties(height=230)
        )

        st.altair_chart(relapse_chart, width="stretch")

    with col3:
        st.markdown("**ER Status**")

        er_chart = (
            alt.Chart(er_chart_data)
            .mark_arc(outerRadius=85)
            .encode(
                theta=alt.Theta("Count:Q"),
                color=alt.Color(
                    "Category:N",
                    title="ER Status",
                    scale=alt.Scale(
                        domain=["ER+", "ER-"],
                        range=["#00436C", "#91C0CB"]
                    )
                ),
                tooltip=[
                    alt.Tooltip("Category:N", title="Status"),
                    alt.Tooltip("Count:Q", title="Samples"),
                    alt.Tooltip(
                        "Percentage:Q",
                        title="Percentage",
                        format=".1f"
                    )
                ]
            )
            .properties(height=230)
        )

        st.altair_chart(er_chart, width="stretch")




# ----------------------------------
# 4. Relapse by ER Status
# ----------------------------------

#st.subheader("Relapse Outcomes by ER Status")


with st.container(border=True):

    st.subheader("Relapse Outcomes by ER Status")

    st.caption(
        "This comparison examines relapse rates among ER-positive "
        "and ER-negative patients."
    )

    # Calculate relapse counts and rates
    er_relapse = df.groupby("ER_status")["Relapse"].agg(
        ["sum", "count"]
    )

    er_relapse["Relapse-free"] = (
        er_relapse["count"] - er_relapse["sum"]
    )

    er_relapse["Relapse Rate"] = (
        er_relapse["sum"] / er_relapse["count"] * 100
    )

    # Chi-square test
    contingency = pd.crosstab(df["ER_status"], df["Relapse"])
    contingency = contingency.reindex(
        index=["ER-", "ER+"],
        columns=[0, 1],
        fill_value=0
    )

    chi2, p_value, dof, expected = chi2_contingency(contingency)

    # Calculate relapse percentages for chart
    chart_data = pd.DataFrame({
        "ER Status": ["ER-", "ER+"],
        "Relapse Rate": [
            er_relapse.loc["ER-", "Relapse Rate"],
            er_relapse.loc["ER+", "Relapse Rate"]
        ]
    })

    # Chart and relapse counts side by side
    chart_col, summary_col = st.columns([2, 1])

    with chart_col:

        relapse_chart = (
            alt.Chart(chart_data)
            .mark_bar()
            .encode(
                x=alt.X(
                    "ER Status:N",
                    title="ER Status",
                    sort=["ER-", "ER+"]
                ),
                y=alt.Y(
                    "Relapse Rate:Q",
                    title="Relapse Rate (%)",
                    scale=alt.Scale(domain=[0, 100])
                ),
                color=alt.Color(
                    "ER Status:N",
                    scale=alt.Scale(
                        domain=["ER-", "ER+"],
                        range=["#00436C", "#91C0CB"]
                    ),
                    legend=None
                ),
                tooltip=[
                    alt.Tooltip("ER Status:N"),
                    alt.Tooltip(
                        "Relapse Rate:Q",
                        title="Relapse Rate",
                        format=".1f"
                    )
                ]
            )
            .properties(height=350)
        )

        st.altair_chart(relapse_chart, width="stretch")

    with summary_col:

        st.markdown("**Relapse Counts**")

        st.metric(
            "ER-negative patients",
            f"{int(er_relapse.loc['ER-', 'sum'])} / "
            f"{int(er_relapse.loc['ER-', 'count'])}"
        )

        st.metric(
            "ER-positive patients",
            f"{int(er_relapse.loc['ER+', 'sum'])} / "
            f"{int(er_relapse.loc['ER+', 'count'])}"
        )

    # Full-width statistical significance summary
    st.divider()

    st.markdown("**Statistical Significance**")

    if p_value < 0.05:
        st.success(
            f"ER status was significantly associated with relapse "
            f"(p = {p_value:.3f})."
        )
    else:
        st.info(
            f"ER status was not significantly associated with "
            f"relapse rates in this cohort (p = {p_value:.3f})."
        )

# ----------------------------------
# 5. Relapse-Free Survival: Survival Analysis
# ----------------------------------

with st.container(border=True):

    st.subheader("Relapse-Free Survival")

    st.caption(
        "Kaplan–Meier estimation of relapse-free survival "
        "over the observed follow-up period."
    )

    # Validate survival data
    survival_columns = ["Relapse", "Follow_up_months"]

    if df[survival_columns].isnull().any().any():
        st.error("Missing values detected in the survival data.")
        st.stop()

    if not df["Relapse"].isin([0, 1]).all():
        st.error("Relapse must contain only 0 and 1.")
        st.stop()

    # Fit Kaplan-Meier estimator
    kmf = KaplanMeierFitter()

    kmf.fit(
        durations=df["Follow_up_months"],
        event_observed=df["Relapse"],
        label="Relapse-free survival"
    )

    # Calculate cohort statistics
    total_patients = len(df)
    relapse_events = int(df["Relapse"].sum())
    censored_patients = total_patients - relapse_events

    median_rfs = kmf.median_survival_time_

    # Summary metrics
    col1, col2, col3 = st.columns(3)

    col1.metric(
        "Patients",
        f"{total_patients:,}"
    )

    col2.metric(
        "Relapse Events",
        f"{relapse_events:,}"
    )

    col3.metric(
        "Median Relapse-Free Survival",
        (
            f"{median_rfs:.1f} months"
            if pd.notna(median_rfs) and median_rfs != float("inf")
            else "Not reached"
        )
    )

    st.divider()


    # Prepare Kaplan-Meier estimates and confidence intervals
    survival_plot_data = pd.DataFrame({
        "Follow_up_months": kmf.survival_function_.index,
        "Survival_probability": kmf.survival_function_.iloc[:, 0].values,
        "CI_lower": kmf.confidence_interval_.iloc[:, 0].values,
        "CI_upper": kmf.confidence_interval_.iloc[:, 1].values
    })

    # Identify censored observations
    censor_times = kmf.event_table[
        kmf.event_table["censored"] > 0
    ].index

    censor_data = pd.DataFrame({
        "Follow_up_months": censor_times,
        "Survival_probability": kmf.predict(censor_times).values
    })

    # Confidence interval
    confidence_band = (
        alt.Chart(survival_plot_data)
        .mark_area(
            interpolate="step-after",
            color="#91C0CB",
            opacity=0.25
        )
        .encode(
            x=alt.X(
                "Follow_up_months:Q",
                title="Follow-up Time (Months)",
                scale=alt.Scale(
                    domain=[0, df["Follow_up_months"].max() + 5]
                )
            ),
            y=alt.Y(
                "CI_lower:Q",
                title="Probability of Remaining Relapse-Free",
                scale=alt.Scale(domain=[0, 1])
            ),
            y2="CI_upper:Q"
        )
    )

    # Kaplan-Meier survival curve
    survival_line = (
        alt.Chart(survival_plot_data)
        .mark_line(
            interpolate="step-after",
            color="#00436C",
            strokeWidth=2.5
        )
        .encode(
            x="Follow_up_months:Q",
            y="Survival_probability:Q"
        )
    )

    # Invisible hover targets along the survival curve
    hover_points = (
        alt.Chart(survival_plot_data)
        .mark_circle(
            size=80,
            opacity=0
        )
        .encode(
            x="Follow_up_months:Q",
            y="Survival_probability:Q",
            tooltip=[
                alt.Tooltip(
                    "Follow_up_months:Q",
                    title="Follow-up",
                    format=".0f"
                ),
                alt.Tooltip(
                    "Survival_probability:Q",
                    title="Relapse-free survival",
                    format=".1%"
                ),
                alt.Tooltip(
                    "CI_lower:Q",
                    title="95% CI lower",
                    format=".1%"
                ),
                alt.Tooltip(
                    "CI_upper:Q",
                    title="95% CI upper",
                    format=".1%"
                )
            ]
        )
    )

    # Censoring indicators
    censor_marks = (
        alt.Chart(censor_data)
        .mark_tick(
            color="#E89583",
            thickness=1.5,
            size=90
        )
        .encode(
            x="Follow_up_months:Q",
            y="Survival_probability:Q",
            tooltip=[
                alt.Tooltip(
                    "Follow_up_months:Q",
                    title="Censoring time",
                    format=".0f"
                ),
                alt.Tooltip(
                    "Survival_probability:Q",
                    title="Relapse-free survival",
                    format=".1%"
                )
            ]
        )
    )

    # Combine interactive chart layers
    survival_chart = (
        confidence_band +
        survival_line +
        hover_points +
        censor_marks
    ).properties(
        height=400,
        title="Kaplan–Meier Relapse-Free Survival Curve"
    ).interactive()

    st.altair_chart(
        survival_chart,
        width="stretch"
    )


    # Survival estimates at selected time points
    st.markdown("**Estimated Relapse-Free Survival**")

    time_points = [36, 60, 84, 120]

    estimate_cols = st.columns(4)

    for col, month in zip(estimate_cols, time_points):

        survival_probability = float(kmf.predict(month))

        col.metric(
            f"{month} Months",
            f"{survival_probability * 100:.1f}%"
        )

    st.caption(
        f"{censored_patients:,} patients had no recorded relapse "
        "during their observed follow-up and were treated as "
        "right-censored observations. The curve displays the "
        "estimated probability of remaining relapse-free."
    )
# ----------------------------------
# 6. Data Quality
# ----------------------------------

with st.container(border=True):

    st.subheader("Data Quality")

    st.caption(
        "Review the missing data in this project"
    )

    total_missing = int(df.isna().sum().sum())
    missing_percentage = (
        total_missing / df.size * 100 if df.size else 0
    )

    complete_rows = int(df.notna().all(axis=1).sum())

    col1, col2, col3 = st.columns(3)

    col1.metric("Total Missing Values", f"{total_missing:,}")

    col2.metric(
        "Missing Data",
        f"{missing_percentage:.2f}%"
    )

    col3.metric(
        "Complete Sample Records",
        f"{complete_rows:,} / {total_samples:,}"
    )

    with st.expander("View missing values by column"):

        missing_by_column = (
            df.isna()
            .sum()
            .reset_index()
        )

        missing_by_column.columns = [
            "Column",
            "Missing Values"
        ]

        missing_by_column = (
            missing_by_column[
                missing_by_column["Missing Values"] > 0
            ]
            .sort_values("Missing Values", ascending=False)
        )

        if missing_by_column.empty:
            st.success("No missing values detected.")
        else:
            st.dataframe(
                missing_by_column,
                width="stretch",
                hide_index=True
            )




# ----------------------------------
# 6. Clinical Data Preview
# ----------------------------------

#with st.container(border=True):

    #st.subheader("Clinical Data Preview")

    #st.caption(
     #   "A preview of the clinical variables available "
     #       "for the samples in this cohort."
   # )

   # st.dataframe(
   #     df[[
   #         "Sample_Accession",
    #        "Relapse",
     #       "ER_status"
      #  ]].head(10),
       # width="stretch",
        #hide_index=True
    #)

    #from scipy.stats import chi2_contingency
    #import numpy as np

