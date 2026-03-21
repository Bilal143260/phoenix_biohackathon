"""
GEM Media Optimizer — Genome-to-media optimization for microbial strains.

A polished, hackathon-ready Streamlit UI for non-specialists and bioprocess engineers.
Helps users upload genome/protein input, generate a GEM, design media compositions,
run growth optimization, and view structured results.

All biological outputs are mocked for MVP demonstration purposes.
"""

import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import json
import time
import io
from datetime import datetime

from bio_utils import parse_fasta, generate_example_fasta, SequenceSummary

# ─────────────────────────────────────────────
# Page config & custom CSS
# ─────────────────────────────────────────────
st.set_page_config(
    page_title="GEM Media Optimizer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# Premium custom CSS
st.markdown("""
<style>
/* ── Global ── */
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

html, body, [class*="css"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
}

.stApp {
    background: linear-gradient(180deg, #f8f9fc 0%, #f0f2f8 100%);
}

/* Hide default Streamlit branding */
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}

/* ── Section containers ── */
.section-card {
    background: #ffffff;
    border-radius: 16px;
    padding: 2rem 2.5rem;
    margin-bottom: 1.5rem;
    box-shadow: 0 1px 3px rgba(0,0,0,0.04), 0 4px 12px rgba(0,0,0,0.03);
    border: 1px solid rgba(0,0,0,0.04);
}

.section-title {
    font-size: 1.35rem;
    font-weight: 600;
    color: #1a1a2e;
    margin-bottom: 0.25rem;
    letter-spacing: -0.01em;
}

.section-subtitle {
    font-size: 0.88rem;
    color: #6b7280;
    margin-bottom: 1.5rem;
    line-height: 1.5;
}

/* ── Hero ── */
.hero-container {
    background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
    border-radius: 20px;
    padding: 3rem 3.5rem 2.5rem;
    margin-bottom: 2rem;
    position: relative;
    overflow: hidden;
}

.hero-container::before {
    content: '';
    position: absolute;
    top: -50%;
    right: -20%;
    width: 500px;
    height: 500px;
    background: radial-gradient(circle, rgba(83,178,166,0.12) 0%, transparent 70%);
    border-radius: 50%;
}

.hero-title {
    font-size: 2.2rem;
    font-weight: 700;
    color: #ffffff;
    margin-bottom: 0.4rem;
    letter-spacing: -0.02em;
}

.hero-tagline {
    font-size: 1.1rem;
    font-weight: 400;
    color: rgba(255,255,255,0.75);
    margin-bottom: 0.5rem;
}

.hero-subtext {
    font-size: 0.85rem;
    color: rgba(255,255,255,0.5);
    max-width: 520px;
}

/* ── Progress stepper ── */
.stepper-container {
    display: flex;
    align-items: center;
    gap: 0;
    margin-top: 1.8rem;
    flex-wrap: wrap;
}

.step {
    display: flex;
    align-items: center;
    gap: 0.5rem;
    padding: 0.45rem 1rem;
    border-radius: 100px;
    font-size: 0.78rem;
    font-weight: 500;
    white-space: nowrap;
}

.step-active {
    background: rgba(83,178,166,0.2);
    color: #53b2a6;
    border: 1px solid rgba(83,178,166,0.3);
}

.step-completed {
    background: rgba(83,178,166,0.15);
    color: #53b2a6;
}

.step-pending {
    color: rgba(255,255,255,0.3);
}

.step-arrow {
    color: rgba(255,255,255,0.15);
    font-size: 0.9rem;
    margin: 0 0.15rem;
}

/* ── KPI cards ── */
.kpi-card {
    background: #ffffff;
    border-radius: 14px;
    padding: 1.4rem 1.6rem;
    border: 1px solid rgba(0,0,0,0.05);
    box-shadow: 0 1px 3px rgba(0,0,0,0.03);
    text-align: center;
}

.kpi-label {
    font-size: 0.75rem;
    font-weight: 500;
    color: #9ca3af;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    margin-bottom: 0.5rem;
}

.kpi-value {
    font-size: 1.75rem;
    font-weight: 700;
    color: #1a1a2e;
    line-height: 1.1;
}

.kpi-unit {
    font-size: 0.8rem;
    font-weight: 400;
    color: #6b7280;
    margin-top: 0.2rem;
}

/* ── Media component cards ── */
.media-card {
    background: #f9fafb;
    border-radius: 12px;
    padding: 1rem 1.25rem;
    margin-bottom: 0.6rem;
    border: 1px solid #e5e7eb;
    display: flex;
    align-items: center;
    justify-content: space-between;
    transition: all 0.15s ease;
}

.media-card:hover {
    border-color: #53b2a6;
    box-shadow: 0 2px 8px rgba(83,178,166,0.08);
}

.media-name {
    font-size: 0.92rem;
    font-weight: 600;
    color: #1a1a2e;
}

/* ── Badges ── */
.badge {
    display: inline-block;
    padding: 0.2rem 0.65rem;
    border-radius: 100px;
    font-size: 0.7rem;
    font-weight: 600;
    letter-spacing: 0.03em;
}

.badge-teal {
    background: rgba(83,178,166,0.12);
    color: #3a9e92;
}

.badge-blue {
    background: rgba(59,130,246,0.1);
    color: #3b82f6;
}

.badge-amber {
    background: rgba(245,158,11,0.1);
    color: #d97706;
}

.badge-green {
    background: rgba(34,197,94,0.1);
    color: #16a34a;
}

/* ── Upload area ── */
.upload-zone {
    border: 2px dashed #d1d5db;
    border-radius: 14px;
    padding: 2.5rem 2rem;
    text-align: center;
    transition: all 0.2s ease;
    background: #fafbfc;
}

.upload-zone:hover {
    border-color: #53b2a6;
    background: rgba(83,178,166,0.02);
}

.upload-icon {
    font-size: 2.2rem;
    margin-bottom: 0.6rem;
}

.upload-text {
    font-size: 0.92rem;
    color: #6b7280;
    margin-bottom: 0.3rem;
}

.upload-hint {
    font-size: 0.78rem;
    color: #9ca3af;
}

/* ── File summary card ── */
.file-summary {
    background: rgba(83,178,166,0.05);
    border: 1px solid rgba(83,178,166,0.15);
    border-radius: 12px;
    padding: 1.2rem 1.5rem;
}

/* ── Buttons ── */
.stButton > button {
    border-radius: 10px !important;
    font-weight: 500 !important;
    padding: 0.5rem 1.5rem !important;
    font-size: 0.88rem !important;
    transition: all 0.2s ease !important;
    border: none !important;
}

div[data-testid="stButton"] > button[kind="primary"] {
    background: linear-gradient(135deg, #53b2a6 0%, #3a9e92 100%) !important;
    color: white !important;
    box-shadow: 0 2px 8px rgba(83,178,166,0.25) !important;
}

div[data-testid="stButton"] > button[kind="primary"]:hover {
    box-shadow: 0 4px 16px rgba(83,178,166,0.35) !important;
    transform: translateY(-1px) !important;
}

/* ── Range cards ── */
.range-card {
    background: #ffffff;
    border-radius: 12px;
    padding: 1rem 1.2rem;
    border: 1px solid #e5e7eb;
    text-align: center;
}

.range-label {
    font-size: 0.82rem;
    font-weight: 600;
    color: #1a1a2e;
    margin-bottom: 0.35rem;
}

.range-value {
    font-size: 0.78rem;
    font-weight: 500;
    color: #53b2a6;
    background: rgba(83,178,166,0.08);
    padding: 0.2rem 0.6rem;
    border-radius: 100px;
    display: inline-block;
}

/* ── Recipe card ── */
.recipe-card {
    background: linear-gradient(135deg, #f0fdf9 0%, #f0f9ff 100%);
    border-radius: 14px;
    padding: 1.5rem 1.8rem;
    border: 1px solid rgba(83,178,166,0.15);
}

/* ── Proxy metric card ── */
.proxy-card {
    background: #ffffff;
    border-radius: 12px;
    padding: 1.2rem 1.4rem;
    border: 1px solid #e5e7eb;
    text-align: center;
}

.proxy-label {
    font-size: 0.75rem;
    font-weight: 500;
    color: #9ca3af;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    margin-bottom: 0.4rem;
}

.proxy-value {
    font-size: 1.4rem;
    font-weight: 700;
    color: #1a1a2e;
}

.proxy-note {
    font-size: 0.7rem;
    color: #9ca3af;
    margin-top: 0.25rem;
}

/* ── Divider ── */
.section-divider {
    height: 1px;
    background: linear-gradient(90deg, transparent 0%, #e5e7eb 50%, transparent 100%);
    margin: 2.5rem 0;
}

/* ── Success state ── */
.success-banner {
    background: rgba(34,197,94,0.06);
    border: 1px solid rgba(34,197,94,0.15);
    border-radius: 12px;
    padding: 1rem 1.4rem;
    display: flex;
    align-items: center;
    gap: 0.6rem;
}

.success-icon {
    color: #16a34a;
    font-size: 1.1rem;
}

.success-text {
    font-size: 0.88rem;
    font-weight: 500;
    color: #15803d;
}

/* ── Streamlit overrides ── */
.stSelectbox > div > div {
    border-radius: 10px !important;
}

.stNumberInput > div > div > input {
    border-radius: 10px !important;
}

div[data-testid="stExpander"] {
    border-radius: 14px !important;
    border: 1px solid #e5e7eb !important;
    box-shadow: none !important;
}

.stProgress > div > div > div {
    background: linear-gradient(90deg, #53b2a6, #3a9e92) !important;
    border-radius: 10px !important;
}

/* ── Tabs ── */
.stTabs [data-baseweb="tab-list"] {
    gap: 0.5rem;
}

.stTabs [data-baseweb="tab"] {
    border-radius: 10px;
    padding: 0.5rem 1rem;
    font-weight: 500;
}

/* ── Empty state ── */
.empty-state {
    text-align: center;
    padding: 3rem 2rem;
    color: #9ca3af;
}

.empty-state-icon {
    font-size: 2.5rem;
    margin-bottom: 0.8rem;
    opacity: 0.4;
}

.empty-state-text {
    font-size: 0.92rem;
    margin-bottom: 0.3rem;
}

.empty-state-hint {
    font-size: 0.8rem;
    color: #c4c9d4;
}
</style>
""", unsafe_allow_html=True)


# ─────────────────────────────────────────────
# Session state initialization
# ─────────────────────────────────────────────
DEFAULTS = {
    "file_uploaded": False,
    "file_name": "",
    "seq_summary": None,  # SequenceSummary from bio_utils
    "gem_generated": False,
    "media_components": [],
    "optimization_run": False,
    "current_step": 0,  # 0=upload, 1=GEM, 2=media, 3=optimize, 4=results
}

for key, val in DEFAULTS.items():
    if key not in st.session_state:
        st.session_state[key] = val

# Available media components (realistic set)
AVAILABLE_COMPONENTS = {
    "Glucose": {"unit": "mM", "default": 20.0, "category": "Carbon source"},
    "Ammonium": {"unit": "mM", "default": 10.0, "category": "Nitrogen source"},
    "Phosphate": {"unit": "mM", "default": 5.0, "category": "Phosphorus source"},
    "Sulfate": {"unit": "mM", "default": 2.0, "category": "Sulfur source"},
    "Oxygen": {"unit": "mmol/gDW/h", "default": 20.0, "category": "Electron acceptor"},
    "Magnesium": {"unit": "mM", "default": 1.0, "category": "Cofactor"},
    "Acetate": {"unit": "mM", "default": 5.0, "category": "Carbon source"},
    "Iron": {"unit": "µM", "default": 50.0, "category": "Trace element"},
    "Potassium": {"unit": "mM", "default": 5.0, "category": "Electrolyte"},
    "Calcium": {"unit": "mM", "default": 0.5, "category": "Cofactor"},
    "L-Glutamine": {"unit": "mM", "default": 2.0, "category": "Amino acid"},
    "Sodium Chloride": {"unit": "mM", "default": 100.0, "category": "Osmolyte"},
}

RECOMMENDED_SET = ["Glucose", "Ammonium", "Phosphate", "Sulfate", "Oxygen", "Magnesium"]


# ─────────────────────────────────────────────
# Helper: compute current workflow step
# ─────────────────────────────────────────────
def get_current_step():
    if st.session_state.optimization_run:
        return 4
    if st.session_state.media_components:
        return 3
    if st.session_state.gem_generated:
        return 2
    if st.session_state.file_uploaded:
        return 1
    return 0


# ─────────────────────────────────────────────
# Mock data generators
# ─────────────────────────────────────────────
MOCK_MODEL_STATS = {
    "reactions": 2583,
    "metabolites": 1877,
    "genes": 1366,
    "compartments": 3,
    "biomass_objective": "Defined",
}

MOCK_OPTIMIZATION_RESULTS = {
    "growth_rate": 0.873,
    "growth_unit": "h⁻¹",
    "active_components": 6,
    "top_driver": "Glucose",
    "space_reduction": "62%",
    "drivers": {
        "Glucose": 0.91,
        "Ammonium": 0.74,
        "Oxygen": 0.68,
        "Phosphate": 0.52,
        "Sulfate": 0.31,
    },
    "ranges": {
        "Glucose": "Medium–High",
        "Ammonium": "Medium",
        "Phosphate": "Low–Medium",
        "Oxygen": "High",
        "Sulfate": "Low",
    },
    "try_benchmarks": {
        "Growth Rate": {"value": "0.873 h⁻¹", "confidence": "High"},
        "Yield Proxy": {"value": "0.41 g/g", "confidence": "Medium"},
        "Product Rate Proxy": {"value": "—", "confidence": "N/A"},
    },
}


def generate_structured_output():
    """Generate realistic structured output JSON for download."""
    components = st.session_state.media_components
    recipe = {}
    for c in components:
        info = AVAILABLE_COMPONENTS.get(c["name"], {})
        recipe[c["name"]] = {
            "concentration": c["quantity"],
            "unit": info.get("unit", "mM"),
        }

    return {
        "gem_media_optimizer": {
            "version": "1.0.0-mvp",
            "timestamp": datetime.now().isoformat(),
            "input": {
                "file": st.session_state.file_name or "sample_proteome.faa",
                "type": (st.session_state.seq_summary.format
                         if st.session_state.seq_summary else "Protein FASTA (.faa)"),
                "sequences": (st.session_state.seq_summary.sequence_count
                              if st.session_state.seq_summary else 0),
            },
            "model": {
                "status": "Draft GEM generated",
                "reactions": MOCK_MODEL_STATS["reactions"],
                "metabolites": MOCK_MODEL_STATS["metabolites"],
                "genes": MOCK_MODEL_STATS["genes"],
                "biomass_objective": "Defined",
            },
            "optimization": {
                "objective": "Maximize Growth",
                "baseline_growth": 0.42,
                "optimized_growth": MOCK_OPTIMIZATION_RESULTS["growth_rate"],
                "improvement": "+107.9%",
            },
            "top_variables": list(MOCK_OPTIMIZATION_RESULTS["drivers"].keys()),
            "recommended_media": recipe,
            "recommended_ranges": MOCK_OPTIMIZATION_RESULTS["ranges"],
            "performance_proxies": MOCK_OPTIMIZATION_RESULTS["try_benchmarks"],
            "experimental_space_reduction": MOCK_OPTIMIZATION_RESULTS["space_reduction"],
        }
    }


# ─────────────────────────────────────────────
# 1. HERO / HEADER
# ─────────────────────────────────────────────
step = get_current_step()
step_labels = ["Upload", "GEM Generation", "Media Design", "Optimization", "Results"]

stepper_html = ""
for i, label in enumerate(step_labels):
    if i > 0:
        stepper_html += '<span class="step-arrow">›</span>'
    if i < step:
        stepper_html += f'<span class="step step-completed">✓ {label}</span>'
    elif i == step:
        stepper_html += f'<span class="step step-active">● {label}</span>'
    else:
        stepper_html += f'<span class="step step-pending">○ {label}</span>'

st.markdown(f"""
<div class="hero-container">
    <div class="hero-title">GEM Media Optimizer</div>
    <div class="hero-tagline">Genome-to-media optimization for microbial strains</div>
    <div class="hero-subtext">
        Upload your protein sequences, generate a metabolic model, and discover optimal
        media compositions — no modeling expertise required.
    </div>
    <div class="stepper-container">
        {stepper_html}
    </div>
</div>
""", unsafe_allow_html=True)


# ─────────────────────────────────────────────
# Sidebar: Project controls
# ─────────────────────────────────────────────
with st.sidebar:
    st.markdown("### Project Controls")

    if st.button("🔄 Reset Project", use_container_width=True):
        for key in DEFAULTS:
            st.session_state[key] = DEFAULTS[key] if not isinstance(DEFAULTS[key], list) else []
        st.rerun()

    st.markdown("---")

    st.markdown("##### Workflow Status")
    for i, label in enumerate(step_labels):
        if i < step:
            st.markdown(f"✅ {label}")
        elif i == step:
            st.markdown(f"🔵 {label}")
        else:
            st.markdown(f"⬜ {label}")

    st.markdown("---")
    st.markdown(
        '<span class="badge badge-teal">MVP</span> '
        '<span class="badge badge-blue">Growth Mode</span>',
        unsafe_allow_html=True,
    )
    st.caption("GEM Media Optimizer v1.0 — Hackathon MVP")


# ─────────────────────────────────────────────
# 2. INPUT / UPLOAD SECTION
# ─────────────────────────────────────────────
st.markdown("""
<div class="section-card">
    <div class="section-title">Project Input</div>
    <div class="section-subtitle">
        Upload your organism's protein sequences to begin. We'll use them to build a
        genome-scale metabolic model automatically.
    </div>
</div>
""", unsafe_allow_html=True)

col_upload, col_summary = st.columns([3, 2], gap="large")

with col_upload:
    if not st.session_state.file_uploaded:
        uploaded_file = st.file_uploader(
            "Upload protein FASTA file",
            type=["faa", "fasta", "fa"],
            help="Preferred format: .faa (protein FASTA). GenBank (.gbk) support coming soon.",
            label_visibility="collapsed",
        )

        if uploaded_file is not None:
            try:
                summary = parse_fasta(uploaded_file, uploaded_file.name)
                st.session_state.file_uploaded = True
                st.session_state.file_name = uploaded_file.name
                st.session_state.seq_summary = summary
                st.rerun()
            except ValueError as exc:
                st.error(str(exc))
        else:
            st.markdown("""
            <div class="upload-zone">
                <div class="upload-icon">📄</div>
                <div class="upload-text">Drag & drop your <strong>.faa</strong> file here</div>
                <div class="upload-hint">or click to browse · Protein FASTA format · GenBank coming soon</div>
            </div>
            """, unsafe_allow_html=True)

            if st.button("📋 Load Example Dataset", type="secondary", use_container_width=True):
                example_bytes, example_name = generate_example_fasta()
                summary = parse_fasta(io.BytesIO(example_bytes), example_name)
                st.session_state.file_uploaded = True
                st.session_state.file_name = example_name
                st.session_state.seq_summary = summary
                st.rerun()
    else:
        st.markdown(f"""
        <div class="success-banner">
            <span class="success-icon">✓</span>
            <span class="success-text">File uploaded successfully</span>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("")

        # File summary card — uses real parsed data
        summ = st.session_state.seq_summary
        fname = st.session_state.file_name
        seq_count = f"{summ.sequence_count:,}" if summ else "—"
        fmt_label = summ.format if summ else "FASTA"
        size_label = summ.size_label if summ else "—"
        avg_len = f"{summ.avg_length:.0f}" if summ else "—"

        st.markdown(f"""
        <div class="file-summary">
            <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:0.6rem;">
                <span style="font-weight:600; color:#1a1a2e; font-size:0.92rem;">📄 {fname}</span>
                <span class="badge badge-green">Ready</span>
            </div>
            <div style="display:flex; gap:2rem; font-size:0.82rem; color:#6b7280;">
                <span><strong>Type:</strong> {fmt_label}</span>
                <span><strong>Sequences:</strong> {seq_count}</span>
                <span><strong>Total:</strong> {size_label}</span>
            </div>
            <div style="display:flex; gap:2rem; font-size:0.82rem; color:#6b7280; margin-top:0.3rem;">
                <span><strong>Avg length:</strong> {avg_len} residues</span>
                <span><strong>Range:</strong> {summ.min_length:,}–{summ.max_length:,}</span>
            </div>
        </div>
        """, unsafe_allow_html=True)

with col_summary:
    if st.session_state.file_uploaded and not st.session_state.gem_generated:
        st.markdown("")
        st.markdown("")
        if st.button("⚡ Generate GEM", type="primary", use_container_width=True):
            progress = st.progress(0, text="Initializing model reconstruction...")
            steps_text = [
                "Annotating protein functions...",
                "Mapping reactions to pathways...",
                "Building stoichiometric matrix...",
                "Adding transport reactions...",
                "Defining biomass objective...",
                "Validating draft model...",
            ]
            for i, txt in enumerate(steps_text):
                time.sleep(0.4)
                progress.progress((i + 1) / len(steps_text), text=txt)
            time.sleep(0.3)
            progress.empty()
            st.session_state.gem_generated = True
            st.rerun()

    elif st.session_state.gem_generated:
        st.markdown(f"""
        <div class="success-banner" style="margin-bottom:1rem;">
            <span class="success-icon">✓</span>
            <span class="success-text">Draft GEM generated</span>
        </div>
        """, unsafe_allow_html=True)

        stats = MOCK_MODEL_STATS
        c1, c2 = st.columns(2)
        with c1:
            st.markdown(f"""
            <div class="kpi-card" style="margin-bottom:0.6rem;">
                <div class="kpi-label">Reactions</div>
                <div class="kpi-value">{stats['reactions']:,}</div>
            </div>
            """, unsafe_allow_html=True)
            st.markdown(f"""
            <div class="kpi-card">
                <div class="kpi-label">Genes</div>
                <div class="kpi-value">{stats['genes']:,}</div>
            </div>
            """, unsafe_allow_html=True)
        with c2:
            st.markdown(f"""
            <div class="kpi-card" style="margin-bottom:0.6rem;">
                <div class="kpi-label">Metabolites</div>
                <div class="kpi-value">{stats['metabolites']:,}</div>
            </div>
            """, unsafe_allow_html=True)
            st.markdown(f"""
            <div class="kpi-card">
                <div class="kpi-label">Biomass Obj.</div>
                <div class="kpi-value" style="font-size:1.1rem; color:#16a34a;">Defined ✓</div>
            </div>
            """, unsafe_allow_html=True)

    elif not st.session_state.file_uploaded:
        st.markdown("""
        <div class="empty-state" style="padding:2rem;">
            <div class="empty-state-icon">🧬</div>
            <div class="empty-state-text">Upload a file to get started</div>
            <div class="empty-state-hint">Model stats will appear here</div>
        </div>
        """, unsafe_allow_html=True)


# ─────────────────────────────────────────────
# 3. MEDIA BUILDER SECTION
# ─────────────────────────────────────────────
if st.session_state.gem_generated:
    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)

    st.markdown("""
    <div class="section-card" style="padding-bottom:1rem;">
        <div class="section-title">Build Your Media</div>
        <div class="section-subtitle">
            Select the nutrients and compounds for your growth medium. Each component
            affects how well your organism can grow. Start with the recommended set
            or customize your own.
        </div>
    </div>
    """, unsafe_allow_html=True)

    # Controls row
    ctrl_col1, ctrl_col2, ctrl_col3 = st.columns([3, 1.5, 1.5], gap="medium")

    existing_names = [c["name"] for c in st.session_state.media_components]
    available = [n for n in AVAILABLE_COMPONENTS if n not in existing_names]

    with ctrl_col1:
        selected_component = st.selectbox(
            "Available components",
            options=available if available else ["— All components added —"],
            label_visibility="collapsed",
            placeholder="Select a media component to add...",
        )

    with ctrl_col2:
        add_clicked = st.button(
            "➕ Add Component",
            type="primary",
            use_container_width=True,
            disabled=(not available),
        )

    with ctrl_col3:
        rec_clicked = st.button(
            "⚡ Recommended Set",
            type="secondary",
            use_container_width=True,
        )

    # Handle add
    if add_clicked and available and selected_component != "— All components added —":
        info = AVAILABLE_COMPONENTS[selected_component]
        st.session_state.media_components.append({
            "name": selected_component,
            "quantity": info["default"],
            "level": "Medium",
        })
        st.rerun()

    # Handle recommended set
    if rec_clicked:
        for name in RECOMMENDED_SET:
            if name not in existing_names:
                info = AVAILABLE_COMPONENTS[name]
                st.session_state.media_components.append({
                    "name": name,
                    "quantity": info["default"],
                    "level": "Medium",
                })
        st.rerun()

    # Display added components
    if st.session_state.media_components:
        st.markdown(
            f"<div style='font-size:0.82rem; color:#6b7280; margin: 0.5rem 0 0.8rem;'>"
            f"{len(st.session_state.media_components)} component(s) in your media</div>",
            unsafe_allow_html=True,
        )

        to_remove = None
        for idx, comp in enumerate(st.session_state.media_components):
            info = AVAILABLE_COMPONENTS.get(comp["name"], {"unit": "mM", "category": ""})
            c1, c2, c3, c4, c5 = st.columns([2.5, 1.5, 1.5, 1.5, 0.5], gap="small")

            with c1:
                st.markdown(
                    f"<div style='padding-top:0.5rem;'>"
                    f"<span style='font-weight:600; color:#1a1a2e;'>{comp['name']}</span>"
                    f"<br/><span style='font-size:0.75rem; color:#9ca3af;'>{info['category']}</span>"
                    f"</div>",
                    unsafe_allow_html=True,
                )

            with c2:
                new_qty = st.number_input(
                    f"qty_{idx}",
                    value=comp["quantity"],
                    min_value=0.0,
                    step=1.0,
                    format="%.1f",
                    label_visibility="collapsed",
                    key=f"qty_{idx}",
                )
                st.session_state.media_components[idx]["quantity"] = new_qty

            with c3:
                st.markdown(
                    f"<div style='padding-top:0.6rem; font-size:0.82rem; color:#6b7280;'>"
                    f"{info['unit']}</div>",
                    unsafe_allow_html=True,
                )

            with c4:
                level = st.selectbox(
                    f"level_{idx}",
                    options=["Low", "Medium", "High"],
                    index=["Low", "Medium", "High"].index(comp["level"]),
                    label_visibility="collapsed",
                    key=f"level_{idx}",
                )
                st.session_state.media_components[idx]["level"] = level

            with c5:
                if st.button("✕", key=f"remove_{idx}", help="Remove component"):
                    to_remove = idx

        if to_remove is not None:
            st.session_state.media_components.pop(to_remove)
            st.rerun()
    else:
        st.markdown("""
        <div class="empty-state">
            <div class="empty-state-icon">🧪</div>
            <div class="empty-state-text">No media components added yet</div>
            <div class="empty-state-hint">
                Add components above or click "Recommended Set" to get started quickly
            </div>
        </div>
        """, unsafe_allow_html=True)


# ─────────────────────────────────────────────
# 4. OPTIMIZATION GOAL + RUN
# ─────────────────────────────────────────────
if st.session_state.gem_generated and st.session_state.media_components:
    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)

    st.markdown("""
    <div class="section-card" style="padding-bottom:1rem;">
        <div class="section-title">Optimization</div>
        <div class="section-subtitle">
            Choose what to optimize for and run the analysis. The optimizer will find
            the best media composition for your selected goal.
        </div>
    </div>
    """, unsafe_allow_html=True)

    opt_col1, opt_col2 = st.columns([2, 3], gap="large")

    with opt_col1:
        st.markdown(
            "<div style='font-size:0.82rem; font-weight:500; color:#6b7280; margin-bottom:0.3rem;'>"
            "Optimization Goal</div>",
            unsafe_allow_html=True,
        )
        goal = st.selectbox(
            "goal",
            options=["Maximize Growth", "Product Rate (coming soon)", "Yield Proxy (coming soon)"],
            label_visibility="collapsed",
        )
        if "coming soon" in goal:
            st.info("This objective will be available in a future release. Growth mode is active.", icon="ℹ️")

        st.markdown(
            '<span class="badge badge-blue">Growth Mode</span> '
            '<span class="badge badge-teal">MVP</span>',
            unsafe_allow_html=True,
        )

    with opt_col2:
        st.markdown("")
        st.markdown("")

        comp_count = len(st.session_state.media_components)
        st.markdown(
            f"<div style='font-size:0.85rem; color:#6b7280; margin-bottom:1rem;'>"
            f"Ready to optimize with <strong>{comp_count}</strong> media component(s)</div>",
            unsafe_allow_html=True,
        )

        if not st.session_state.optimization_run:
            if st.button("🚀 Run Optimization", type="primary", use_container_width=True):
                progress = st.progress(0, text="Preparing optimization...")
                opt_steps = [
                    "Setting flux constraints...",
                    "Running flux balance analysis...",
                    "Performing sensitivity analysis...",
                    "Ranking media drivers...",
                    "Computing recommended ranges...",
                    "Generating results...",
                ]
                for i, txt in enumerate(opt_steps):
                    time.sleep(0.45)
                    progress.progress((i + 1) / len(opt_steps), text=txt)
                time.sleep(0.3)
                progress.empty()
                st.session_state.optimization_run = True
                st.rerun()
        else:
            st.markdown("""
            <div class="success-banner">
                <span class="success-icon">✓</span>
                <span class="success-text">Optimization complete — scroll down for results</span>
            </div>
            """, unsafe_allow_html=True)


# ─────────────────────────────────────────────
# 5. RESULTS SECTION
# ─────────────────────────────────────────────
if st.session_state.optimization_run:
    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)

    results = MOCK_OPTIMIZATION_RESULTS

    # ── A. KPI Cards ──
    st.markdown("""
    <div class="section-card" style="padding-bottom:1.5rem;">
        <div class="section-title">Optimization Results</div>
        <div class="section-subtitle">
            Key performance indicators from the metabolic model optimization.
        </div>
    </div>
    """, unsafe_allow_html=True)

    k1, k2, k3, k4 = st.columns(4, gap="medium")
    with k1:
        st.markdown(f"""
        <div class="kpi-card">
            <div class="kpi-label">Predicted Growth</div>
            <div class="kpi-value" style="color:#16a34a;">{results['growth_rate']}</div>
            <div class="kpi-unit">{results['growth_unit']}</div>
        </div>
        """, unsafe_allow_html=True)
    with k2:
        st.markdown(f"""
        <div class="kpi-card">
            <div class="kpi-label">Active Components</div>
            <div class="kpi-value">{results['active_components']}</div>
            <div class="kpi-unit">in media</div>
        </div>
        """, unsafe_allow_html=True)
    with k3:
        st.markdown(f"""
        <div class="kpi-card">
            <div class="kpi-label">Top Driver</div>
            <div class="kpi-value" style="font-size:1.3rem;">{results['top_driver']}</div>
            <div class="kpi-unit">most influential</div>
        </div>
        """, unsafe_allow_html=True)
    with k4:
        st.markdown(f"""
        <div class="kpi-card">
            <div class="kpi-label">Space Reduction</div>
            <div class="kpi-value" style="color:#3b82f6;">{results['space_reduction']}</div>
            <div class="kpi-unit">experimental space</div>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("")

    # ── B. Key Media Drivers Chart ──
    chart_col, range_col = st.columns([3, 2], gap="large")

    with chart_col:
        st.markdown("""
        <div class="section-title" style="font-size:1.1rem;">Key Media Drivers</div>
        <div class="section-subtitle" style="margin-bottom:1rem;">
            Components ranked by their influence on growth rate.
        </div>
        """, unsafe_allow_html=True)

        drivers = results["drivers"]
        driver_names = list(drivers.keys())[::-1]
        driver_values = [drivers[n] for n in driver_names]

        # Color gradient based on values
        colors = [
            f"rgba(83,178,166,{0.3 + 0.7 * (v / max(driver_values))})"
            for v in driver_values
        ]

        fig = go.Figure(go.Bar(
            x=driver_values,
            y=driver_names,
            orientation="h",
            marker=dict(
                color=colors,
                line=dict(width=0),
                cornerradius=6,
            ),
            text=[f"{v:.2f}" for v in driver_values],
            textposition="outside",
            textfont=dict(size=12, color="#374151", family="Inter"),
        ))
        fig.update_layout(
            height=260,
            margin=dict(l=0, r=40, t=10, b=10),
            xaxis=dict(
                title=dict(text="Sensitivity Score", font=dict(size=11, color="#9ca3af")),
                range=[0, 1.1],
                showgrid=True,
                gridcolor="rgba(0,0,0,0.04)",
                zeroline=False,
                tickfont=dict(size=10, color="#9ca3af"),
            ),
            yaxis=dict(
                tickfont=dict(size=12, color="#374151", family="Inter"),
                showgrid=False,
            ),
            plot_bgcolor="rgba(0,0,0,0)",
            paper_bgcolor="rgba(0,0,0,0)",
            font=dict(family="Inter"),
        )
        st.plotly_chart(fig, use_container_width=True, config={"displayModeBar": False})

    # ── C. Recommended Media Window ──
    with range_col:
        st.markdown("""
        <div class="section-title" style="font-size:1.1rem;">Recommended Ranges</div>
        <div class="section-subtitle" style="margin-bottom:1rem;">
            Suggested concentration levels for experimental design.
        </div>
        """, unsafe_allow_html=True)

        ranges = results["ranges"]
        for comp_name, range_val in ranges.items():
            level_color = {
                "Low": "#ef4444",
                "Low–Medium": "#f97316",
                "Medium": "#eab308",
                "Medium–High": "#22c55e",
                "High": "#3b82f6",
            }.get(range_val, "#53b2a6")

            st.markdown(f"""
            <div class="range-card" style="margin-bottom:0.5rem; text-align:left; display:flex; justify-content:space-between; align-items:center;">
                <span class="range-label" style="margin-bottom:0;">{comp_name}</span>
                <span class="range-value" style="background:rgba({','.join(str(int(level_color.lstrip('#')[i:i+2], 16)) for i in (0, 2, 4))},0.1); color:{level_color};">{range_val}</span>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("")

    # ── D. Recommended Recipe ──
    rec_col, try_col = st.columns([3, 2], gap="large")

    with rec_col:
        st.markdown("""
        <div class="section-title" style="font-size:1.1rem;">Recommended Recipe</div>
        <div class="section-subtitle">
            Optimized media composition based on the model analysis.
        </div>
        """, unsafe_allow_html=True)

        st.markdown('<div class="recipe-card">', unsafe_allow_html=True)

        recipe_data = []
        for comp in st.session_state.media_components:
            info = AVAILABLE_COMPONENTS.get(comp["name"], {"unit": "mM"})
            # Simulate optimized quantities (slightly adjusted from input)
            optimized_qty = round(comp["quantity"] * (0.8 + 0.5 * (results["drivers"].get(comp["name"], 0.3))), 1)
            recipe_data.append({
                "Component": comp["name"],
                "Optimized": f"{optimized_qty} {info['unit']}",
                "Original": f"{comp['quantity']} {info['unit']}",
                "Range": results["ranges"].get(comp["name"], "Medium"),
            })

        if recipe_data:
            for item in recipe_data:
                st.markdown(
                    f"<div style='display:flex; justify-content:space-between; align-items:center; "
                    f"padding:0.5rem 0; border-bottom:1px solid rgba(0,0,0,0.04);'>"
                    f"<span style='font-weight:600; color:#1a1a2e; font-size:0.88rem;'>{item['Component']}</span>"
                    f"<div style='display:flex; gap:1rem; align-items:center;'>"
                    f"<span style='font-size:0.82rem; color:#6b7280;'>{item['Original']} → </span>"
                    f"<span style='font-size:0.88rem; font-weight:600; color:#16a34a;'>{item['Optimized']}</span>"
                    f"<span class='badge badge-teal' style='font-size:0.68rem;'>{item['Range']}</span>"
                    f"</div></div>",
                    unsafe_allow_html=True,
                )

        st.markdown("</div>", unsafe_allow_html=True)

    # ── E. TRY Benchmarking (Performance Proxies) ──
    with try_col:
        st.markdown("""
        <div class="section-title" style="font-size:1.1rem;">Performance Proxies</div>
        <div class="section-subtitle">
            Model-based estimates — not wet-lab measurements.
        </div>
        """, unsafe_allow_html=True)

        benchmarks = results["try_benchmarks"]
        for metric, data in benchmarks.items():
            confidence_color = {
                "High": "#16a34a",
                "Medium": "#eab308",
                "N/A": "#9ca3af",
            }.get(data["confidence"], "#6b7280")

            st.markdown(f"""
            <div class="proxy-card" style="margin-bottom:0.6rem;">
                <div class="proxy-label">{metric}</div>
                <div class="proxy-value">{data['value']}</div>
                <div class="proxy-note">
                    Confidence: <span style="color:{confidence_color}; font-weight:600;">{data['confidence']}</span>
                </div>
            </div>
            """, unsafe_allow_html=True)

        st.markdown(
            "<div style='font-size:0.72rem; color:#9ca3af; margin-top:0.6rem; padding:0 0.2rem; line-height:1.5;'>"
            "These are simplified proxies derived from the metabolic model. "
            "They are not predictions of wet-lab titers.</div>",
            unsafe_allow_html=True,
        )

    # ── F. Structured Output ──
    st.markdown("")
    st.markdown('<div class="section-divider"></div>', unsafe_allow_html=True)

    st.markdown("""
    <div class="section-title" style="font-size:1.1rem;">Structured Output</div>
    <div class="section-subtitle">
        Download or inspect the full analysis results.
    </div>
    """, unsafe_allow_html=True)

    structured_data = generate_structured_output()
    json_str = json.dumps(structured_data, indent=2)

    # Cache download data in session state to prevent stale media file errors
    st.session_state["_dl_json"] = json_str

    csv_rows = []
    for comp in st.session_state.media_components:
        info = AVAILABLE_COMPONENTS.get(comp["name"], {"unit": "mM"})
        optimized_qty = round(comp["quantity"] * (0.8 + 0.5 * (results["drivers"].get(comp["name"], 0.3))), 1)
        csv_rows.append({
            "Component": comp["name"],
            "Original_Concentration": comp["quantity"],
            "Optimized_Concentration": optimized_qty,
            "Unit": info["unit"],
            "Recommended_Range": results["ranges"].get(comp["name"], "Medium"),
            "Sensitivity_Score": results["drivers"].get(comp["name"], 0.0),
        })
    csv_df = pd.DataFrame(csv_rows)
    st.session_state["_dl_csv"] = csv_df.to_csv(index=False)

    # Download buttons
    dl_col1, dl_col2, dl_col3 = st.columns([1, 1, 3], gap="medium")

    with dl_col1:
        st.download_button(
            "📥 Download JSON",
            data=st.session_state["_dl_json"],
            file_name="gem_media_optimizer_results.json",
            mime="application/json",
            use_container_width=True,
        )

    with dl_col2:
        st.download_button(
            "📥 Download CSV",
            data=st.session_state["_dl_csv"],
            file_name="gem_media_optimizer_results.csv",
            mime="text/csv",
            use_container_width=True,
        )

    # Collapsible JSON view
    with st.expander("View full JSON output", expanded=False):
        st.json(structured_data)


# ─────────────────────────────────────────────
# Footer
# ─────────────────────────────────────────────
st.markdown("")
st.markdown("")
st.markdown(
    "<div style='text-align:center; font-size:0.75rem; color:#c4c9d4; padding:2rem 0 1rem;'>"
    "GEM Media Optimizer · Hackathon MVP · Built with Streamlit & COBRA"
    "</div>",
    unsafe_allow_html=True,
)
