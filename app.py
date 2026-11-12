import streamlit as st

st.set_page_config(
    page_title="GSE2034 Research Explorer",
    page_icon="🧬",
    layout="wide"
)

# Define dashboard pages
cohort_page = st.Page(
    "pages/1_Cohort_Overview.py",
    title="Cohort Overview",
    icon="👥",
    default=True
)

expression_page = st.Page(
    "pages/2_Gene_Expression_Patterns.py",
    title="Gene Expression Patterns",
    icon="🧬"
)

survival_page = st.Page(
    "pages/3_Survival_Analysis.py",
    title="Survival Analysis",
    icon="📈"
)

prediction_page = st.Page(
    "pages/4_Predictive_Modeling.py",
    title="Predictive Modeling",
    icon="🤖"
)

research_page = st.Page(
    "pages/5_Research.py",
    title="Research Questions & Findings",
    icon="🔬"
)

# Create navigation
navigation = st.navigation([
    cohort_page,
    expression_page,
    survival_page,
    prediction_page,
    research_page
])

navigation.run()
