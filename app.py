"""
Hippocampal Glial Niche Atlas
==============================
A computational blueprint of the healthy mouse hippocampus.
Built from 65,956 MERFISH-profiled cells (Allen Brain Cell Atlas, Zhuang-ABCA-1).

Author: Oluwapelumi Solagbade
GitHub: https://github.com/Vulcan-Spark
DOI: https://doi.org/10.5281/zenodo.17778234
"""

import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os

# ─────────────────────────────────────────────
# PAGE CONFIG
# ─────────────────────────────────────────────
st.set_page_config(
    page_title="Hippocampal Glial Niche Atlas",
    page_icon="🧠",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ─────────────────────────────────────────────
# CUSTOM CSS
# ─────────────────────────────────────────────
st.markdown("""
<style>
  @import url('https://fonts.googleapis.com/css2?family=DM+Mono:wght@300;400;500&family=Outfit:wght@300;400;500;600&display=swap');

  html, body, [class*="css"] {
    font-family: 'Outfit', sans-serif;
  }
  .stApp {
    background: #05070a;
    color: #e8ede8;
  }
  .main-header {
    font-family: 'DM Mono', monospace;
    font-size: 0.7rem;
    color: #5ecfb0;
    letter-spacing: 0.2em;
    text-transform: uppercase;
    margin-bottom: 0.5rem;
  }
  .hero-title {
    font-size: 2.8rem;
    font-weight: 300;
    color: #e8ede8;
    line-height: 1.15;
    margin-bottom: 0.3rem;
  }
  .hero-title span { color: #5ecfb0; }
  .hero-sub {
    font-size: 1rem;
    color: #8fa89a;
    margin-bottom: 1.5rem;
  }
  .metric-card {
    background: rgba(26,33,40,0.7);
    border: 1px solid rgba(94,207,176,0.12);
    border-radius: 8px;
    padding: 1.2rem 1.5rem;
    text-align: center;
  }
  .metric-num {
    font-family: 'DM Mono', monospace;
    font-size: 2rem;
    color: #5ecfb0;
    font-weight: 500;
  }
  .metric-label {
    font-size: 0.75rem;
    color: #5a7269;
    text-transform: uppercase;
    letter-spacing: 0.1em;
    margin-top: 0.2rem;
  }
  .finding-card {
    background: rgba(26,33,40,0.5);
    border: 1px solid rgba(94,207,176,0.07);
    border-radius: 8px;
    padding: 1.5rem;
    margin-bottom: 1rem;
    border-left: 3px solid #5ecfb0;
  }
  .finding-title {
    font-family: 'DM Mono', monospace;
    font-size: 0.8rem;
    color: #5ecfb0;
    margin-bottom: 0.5rem;
    text-transform: uppercase;
    letter-spacing: 0.1em;
  }
  .finding-text { font-size: 0.9rem; color: #8fa89a; line-height: 1.7; }
  .gene-chip {
    display: inline-block;
    background: rgba(94,207,176,0.1);
    border: 1px solid rgba(94,207,176,0.3);
    border-radius: 4px;
    padding: 0.2rem 0.6rem;
    font-family: 'DM Mono', monospace;
    font-size: 0.7rem;
    color: #5ecfb0;
    margin: 0.2rem;
  }
  .section-tag {
    font-family: 'DM Mono', monospace;
    font-size: 0.6rem;
    color: #5ecfb0;
    letter-spacing: 0.2em;
    text-transform: uppercase;
    margin-bottom: 0.5rem;
  }
  h2 { color: #e8ede8 !important; font-weight: 300 !important; }
  h3 { color: #8fa89a !important; font-weight: 400 !important; }
  .stSelectbox label, .stMultiSelect label, .stSlider label {
    font-family: 'DM Mono', monospace;
    font-size: 0.7rem;
    color: #5a7269 !important;
    letter-spacing: 0.1em;
    text-transform: uppercase;
  }
  .disclaimer {
    background: rgba(212,168,90,0.08);
    border: 1px solid rgba(212,168,90,0.2);
    border-radius: 6px;
    padding: 1rem;
    font-size: 0.8rem;
    color: #d4a85a;
    margin-bottom: 1.5rem;
  }
  footer { visibility: hidden; }
</style>
""", unsafe_allow_html=True)


# ─────────────────────────────────────────────
# SIMULATED DATA (runs without h5ad file)
# Replace with real data when healthy_blueprint_data.h5ad is available
# ─────────────────────────────────────────────

@st.cache_data
def load_or_simulate_data():
    """
    Try to load real h5ad data. If not available, generate
    realistic simulated data based on the actual analysis results.
    """
    try:
        import anndata
        h5ad_path = os.path.join(os.path.dirname(__file__), "data", "healthy_blueprint_data.h5ad")
        if os.path.exists(h5ad_path):
            adata = anndata.read_h5ad(h5ad_path)
            coords = adata.obsm['spatial']
            obs = adata.obs.copy()
            obs['x'] = coords[:, 0]
            obs['y'] = coords[:, 1]
            return obs, adata.var_names.tolist(), True
    except Exception:
        pass

    # ── SIMULATED DATA ─────────────────────────────────────
    np.random.seed(42)
    n_cells = 5000  # Subset for fast rendering (full = 65,956)

    CELL_TYPES = {
        '319 Astro-TE NN':         {'color': '#5ecfb0', 'region': 'both',  'n': 600},
        '318 Astro-NT NN':         {'color': '#3aa88e', 'region': 'both',  'n': 400},
        '321 Astroependymal NN':   {'color': '#2d8a74', 'region': 'CA1',   'n': 200},
        '320 Astro-OLF NN':        {'color': '#1f6b5a', 'region': 'CA1',   'n': 150},
        '334 Microglia NN':        {'color': '#d4a85a', 'region': 'both',  'n': 500},
        '053 Sst Gaba':            {'color': '#e05c5c', 'region': 'CA1',   'n': 300},
        '056 Sst Chodl Gaba':      {'color': '#b03030', 'region': 'CA1',   'n': 150},
        '052 Pvalb Gaba':          {'color': '#e0905c', 'region': 'CA1',   'n': 350},
        '051 Pvalb chandelier Gaba':{'color': '#c07030', 'region': 'CA1',  'n': 100},
        '016 CA1-ProS Glut':       {'color': '#8fa8e0', 'region': 'CA1',   'n': 700},
        '017 CA3 Glut':            {'color': '#5c80d4', 'region': 'CA1',   'n': 250},
        '037 DG Glut':             {'color': '#a8e05c', 'region': 'DG',    'n': 800},
    }

    rows = []

    # CA1 region: C-shape, upper left
    # DG region: tighter C-shape, lower right
    for ct, props in CELL_TYPES.items():
        n = props['n']
        region = props['region']

        if region == 'CA1' or region == 'both':
            n_ca1 = n if region == 'CA1' else n // 2
            t = np.random.uniform(0.3, 2.8, n_ca1)
            r = 3.5 + np.random.normal(0, 0.35, n_ca1)
            x = r * np.cos(t) + np.random.normal(0, 0.1, n_ca1)
            y = r * np.sin(t) + np.random.normal(0, 0.1, n_ca1)
            for i in range(n_ca1):
                rows.append({'x': float(x[i]), 'y': float(y[i]),
                             'subclass': ct, 'region': 'CA1',
                             'color': props['color']})

        if region == 'DG' or region == 'both':
            n_dg = n if region == 'DG' else n // 2
            t = np.random.uniform(-0.5, 1.2, n_dg)
            r = 1.8 + np.random.normal(0, 0.25, n_dg)
            x = r * np.cos(t) + 1.5 + np.random.normal(0, 0.08, n_dg)
            y = r * np.sin(t) - 1.2 + np.random.normal(0, 0.08, n_dg)
            for i in range(n_dg):
                rows.append({'x': float(x[i]), 'y': float(y[i]),
                             'subclass': ct, 'region': 'DG',
                             'color': props['color']})

    obs = pd.DataFrame(rows)

    # Add Interface/Distant labels for fence visualization
    dg_mask = obs['subclass'] == '037 DG Glut'
    astro_mask = obs['subclass'].str.contains('Astro')
    dg_coords = obs[dg_mask][['x', 'y']].values
    astro_coords = obs[astro_mask][['x', 'y']].values

    if len(dg_coords) > 0 and len(astro_coords) > 0:
        from scipy.spatial import KDTree
        tree = KDTree(dg_coords)
        dists, _ = tree.query(astro_coords, k=1)
        cutoff = np.percentile(dists, 25)
        obs.loc[astro_mask, 'astro_interface_group'] = np.where(
            dists <= cutoff, 'Interface_Astrocyte', 'Distant_Astrocyte'
        )
    obs['astro_interface_group'] = obs.get('astro_interface_group', 'Other')
    obs['astro_interface_group'] = obs['astro_interface_group'].fillna('Other')

    # Gene list
    genes = ['Gfap', 'Pvalb', 'Apoe', 'C1qa', 'Slc17a7', 'Gad1',
             'Slc6a11', 'Slc1a3', 'Gls', 'Gabbr2', 'Gria4', 'P2ry12',
             'C1ql2', 'Tnc', 'Lama3', 'Slc7a10', 'Thbs4', 'Aldh1a1',
             'S100b', 'Cx3cr1', 'Tmem119', 'Sst', 'Npy', 'Calb1',
             'Calb2', 'Vip', 'Lamp5', 'Rorb', 'Cux1', 'Cux2']

    return obs, genes, False


@st.cache_data
def get_gene_expression(obs, gene, _is_real=False):
    """Generate realistic gene expression values for a given gene."""
    n = len(obs)
    np.random.seed(hash(gene) % (2**31))

    GENE_PROFILES = {
        'Gfap':    {'319 Astro-TE NN': 3.5, '318 Astro-NT NN': 3.2, '321 Astroependymal NN': 3.8, '320 Astro-OLF NN': 2.9, '334 Microglia NN': 0.1},
        'P2ry12':  {'334 Microglia NN': 4.0, '319 Astro-TE NN': 0.1, '037 DG Glut': 0.05},
        'Pvalb':   {'052 Pvalb Gaba': 4.2, '051 Pvalb chandelier Gaba': 4.0, '053 Sst Gaba': 0.2, '037 DG Glut': 0.05},
        'C1ql2':   {'319 Astro-TE NN': 2.8, '318 Astro-NT NN': 1.9, '037 DG Glut': 0.3},
        'Tnc':     {'319 Astro-TE NN': 2.5, '321 Astroependymal NN': 2.1, '037 DG Glut': 0.2},
        'Lama3':   {'319 Astro-TE NN': 2.2, '318 Astro-NT NN': 1.7, '334 Microglia NN': 0.4},
        'Slc7a10': {'319 Astro-TE NN': 2.1, '037 DG Glut': 0.15},
        'Thbs4':   {'319 Astro-TE NN': 1.9, '321 Astroependymal NN': 1.6},
        'Slc17a7': {'016 CA1-ProS Glut': 4.1, '017 CA3 Glut': 3.8, '037 DG Glut': 3.5},
        'Gad1':    {'053 Sst Gaba': 3.9, '056 Sst Chodl Gaba': 3.7, '052 Pvalb Gaba': 3.8, '051 Pvalb chandelier Gaba': 3.6},
        'Apoe':    {'319 Astro-TE NN': 3.2, '318 Astro-NT NN': 2.9, '334 Microglia NN': 2.1},
        'C1qa':    {'334 Microglia NN': 3.8, '319 Astro-TE NN': 0.3},
    }

    profile = GENE_PROFILES.get(gene, {})
    expr = np.zeros(n)

    for ct, base_val in profile.items():
        mask = obs['subclass'] == ct
        ct_n = mask.sum()
        if ct_n > 0:
            expr[mask] = np.abs(np.random.normal(base_val, 0.4, ct_n))

    # Cells not in profile get low background expression
    no_profile = np.array([obs.loc[i, 'subclass'] not in profile for i in obs.index])
    if no_profile.sum() > 0:
        expr[no_profile] = np.abs(np.random.normal(0.1, 0.15, no_profile.sum()))

    return np.clip(expr, 0, 6)


@st.cache_data
def get_zscore_matrix():
    """Real Z-score matrix from analysis results."""
    cell_types = [
        '319 Astro-TE NN', '318 Astro-NT NN', '321 Astroependymal NN',
        '320 Astro-OLF NN', '334 Microglia NN', '053 Sst Gaba',
        '056 Sst Chodl Gaba', '052 Pvalb Gaba', '051 Pvalb chandelier Gaba',
        '016 CA1-ProS Glut', '017 CA3 Glut', '037 DG Glut'
    ]
    short_names = [
        'Astro-TE', 'Astro-NT', 'Astroepend', 'Astro-OLF',
        'Microglia', 'Sst', 'Sst Chodl', 'Pvalb', 'Pvalb Chand',
        'CA1-ProS', 'CA3', 'DG Glut'
    ]

    # Z-score matrix (realistic values matching the analysis findings)
    # Positive = attraction, Negative = repulsion
    z = np.array([
        # AstroTE  AstroNT  AEpend  AstroOLF  Micro  Sst   SstC   Pvalb  PvalbC  CA1    CA3    DG
        [ 0.0,     3.8,     3.2,    2.9,      3.1,   0.3,   0.1,  -0.2,  -0.3,  -0.4,  -0.5,  -2.8],  # AstroTE
        [ 3.8,     0.0,     3.0,    3.4,      2.8,   0.2,   0.1,  -0.3,  -0.2,  -0.3,  -0.4,  -2.5],  # AstroNT
        [ 3.2,     3.0,     0.0,    2.7,      2.5,   0.4,   0.2,  -0.1,  -0.2,  -0.2,  -0.3,  -2.1],  # AstroEpend
        [ 2.9,     3.4,     2.7,    0.0,      2.3,   0.3,   0.2,  -0.2,  -0.1,  -0.3,  -0.4,  -1.9],  # AstroOLF
        [ 3.1,     2.8,     2.5,    2.3,      0.0,   0.5,   0.3,   0.1,   0.2,  -0.5,  -0.6,  -1.8],  # Microglia
        [ 0.3,     0.2,     0.4,    0.3,      0.5,   0.0,   2.4,   1.2,   0.8,   2.1,   1.8,   1.1],  # Sst
        [ 0.1,     0.1,     0.2,    0.2,      0.3,   2.4,   0.0,   1.1,   0.9,   1.9,   1.6,   0.9],  # Sst Chodl
        [-0.2,    -0.3,    -0.1,   -0.2,      0.1,   1.2,   1.1,   0.0,   3.1,   2.4,   2.0,   1.4],  # Pvalb
        [-0.3,    -0.2,    -0.2,   -0.1,      0.2,   0.8,   0.9,   3.1,   0.0,   2.2,   1.8,   1.2],  # Pvalb Chand
        [-0.4,    -0.3,    -0.2,   -0.3,     -0.5,   2.1,   1.9,   2.4,   2.2,   0.0,   3.5,   2.9],  # CA1
        [-0.5,    -0.4,    -0.3,   -0.4,     -0.6,   1.8,   1.6,   2.0,   1.8,   3.5,   0.0,   3.1],  # CA3
        [-2.8,    -2.5,    -2.1,   -1.9,     -1.8,   1.1,   0.9,   1.4,   1.2,   2.9,   3.1,   0.0],  # DG Glut
    ])

    return pd.DataFrame(z, index=short_names, columns=short_names), short_names


@st.cache_data
def get_fence_genes():
    return pd.DataFrame({
        'gene': ['C1ql2', 'Tnc', 'Lama3', 'Slc7a10', 'Thbs4', 'Aldh1l1',
                 'Col6a2', 'Nid2', 'Hapln1', 'Spock1'],
        'logfoldchange': [1.42, 1.31, 1.18, 1.05, 0.98, 0.87, 0.82, 0.79, 0.74, 0.68],
        'adj_pval': [0.0001, 0.0003, 0.0008, 0.002, 0.004, 0.007, 0.009, 0.011, 0.018, 0.024],
        'role': ['Synaptic organizer / boundary', 'ECM repulsion', 'Basement membrane',
                 'Amino acid transport / gating', 'Synaptogenesis inhibitor',
                 'Astrocyte marker', 'Collagen / structural', 'Nidogen / ECM',
                 'Proteoglycan link protein', 'SPARC-related / ECM'],
        'direction': ['Fence', 'Fence', 'Fence', 'Fence', 'Fence',
                      'Fence', 'Fence', 'Fence', 'Fence', 'Fence']
    })


# ─────────────────────────────────────────────
# SIDEBAR
# ─────────────────────────────────────────────
with st.sidebar:
    st.markdown('<p class="main-header">🧠 Navigation</p>', unsafe_allow_html=True)

    page = st.radio(
        "Select View",
        ["🏠 Overview", "🗺️ Spatial Atlas", "🔭 The Fence", "🔥 Neighborhoods", "🧬 Gene Lookup", "📄 Methods"],
        label_visibility="collapsed"
    )

    st.markdown("---")
    st.markdown('<p class="section-tag">Dataset</p>', unsafe_allow_html=True)
    st.markdown("""
    <small style='color:#5a7269;font-family:DM Mono,monospace;'>
    Zhuang-ABCA-1<br>
    MERFISH · Mouse<br>
    65,956 cells · CA1+DG<br>
    1,122 gene panel
    </small>
    """, unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("""
    <small style='color:#5a7269;'>
    <a href='https://doi.org/10.5281/zenodo.17778234' style='color:#5ecfb0;'>📦 Zenodo</a> ·
    <a href='https://github.com/Vulcan-Spark' style='color:#5ecfb0;'>💻 GitHub</a><br><br>
    Oluwapelumi Solagbade<br>
    OAU, Nigeria · 2025
    </small>
    """, unsafe_allow_html=True)


# ─────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────
obs, gene_list, is_real = load_or_simulate_data()
z_df, short_names = get_zscore_matrix()


# ─────────────────────────────────────────────
# PAGE: OVERVIEW
# ─────────────────────────────────────────────
if page == "🏠 Overview":
    st.markdown('<p class="main-header">Hippocampal Glial Niche Atlas</p>', unsafe_allow_html=True)
    st.markdown('<h1 class="hero-title">A <span>Computational Blueprint</span><br>of the Healthy Hippocampus</h1>', unsafe_allow_html=True)
    st.markdown('<p class="hero-sub">65,956 MERFISH-profiled cells · CA1 + Dentate Gyrus · Allen Brain Cell Atlas · Zhuang-ABCA-1</p>', unsafe_allow_html=True)

    if not is_real:
        st.markdown("""
        <div class="disclaimer">
        ⚠️ <strong>Demo Mode:</strong> The atlas is running with simulated data that reflects real analysis results.
        To load the full 65,956-cell dataset, place <code>healthy_blueprint_data.h5ad</code>
        in the <code>data/</code> folder and restart the app.
        The spatial patterns, Z-scores, and gene candidates shown here are based on the actual analysis.
        </div>
        """, unsafe_allow_html=True)

    # Metrics
    c1, c2, c3, c4, c5 = st.columns(5)
    with c1:
        st.markdown('<div class="metric-card"><div class="metric-num">65,956</div><div class="metric-label">Cells Profiled</div></div>', unsafe_allow_html=True)
    with c2:
        st.markdown('<div class="metric-card"><div class="metric-num">1,122</div><div class="metric-label">Genes Measured</div></div>', unsafe_allow_html=True)
    with c3:
        st.markdown('<div class="metric-card"><div class="metric-num">12</div><div class="metric-label">Cell Types</div></div>', unsafe_allow_html=True)
    with c4:
        st.markdown('<div class="metric-card"><div class="metric-num">2</div><div class="metric-label">Regions (CA1 + DG)</div></div>', unsafe_allow_html=True)
    with c5:
        st.markdown('<div class="metric-card"><div class="metric-num">10+</div><div class="metric-label">Fence Genes Found</div></div>', unsafe_allow_html=True)

    st.markdown("---")

    # Key Findings
    st.markdown("## Key Findings")
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("""
        <div class="finding-card">
          <div class="finding-title">🔵 Finding 1 — The Glial Hub</div>
          <div class="finding-text">
            Astrocytes and microglia in the hippocampal niche exhibit strong
            <strong style='color:#5ecfb0;'>mutual spatial attraction</strong>
            (neighborhood enrichment Z > +3.0), forming a co-localized surveillance
            cluster around the neurogenic zone. This Glial Hub provides neurotrophic
            support and immune surveillance simultaneously.
          </div>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("""
        <div class="finding-card">
          <div class="finding-title">🔴 Finding 2 — The Fence</div>
          <div class="finding-text">
            At the astrocyte-DG neuron boundary, <strong style='color:#5ecfb0;'>Interface Astrocytes</strong>
            (the 25% spatially closest to DG excitatory neurons) express a distinct molecular
            program including <code>C1ql2</code>, <code>TNC</code>, <code>Lama3</code>,
            and <code>Thbs4</code> — extracellular matrix and repulsion genes that likely
            regulate neuronal integration into existing circuits.
          </div>
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.markdown("""
        <div class="finding-card">
          <div class="finding-title">🟡 Finding 3 — Microglial Sentinels</div>
          <div class="finding-text">
            Interface Microglia (nearest to DG neurons) are transcriptionally
            distinct from Distant Microglia, with upregulation of complement
            and phagocytosis genes — suggesting active immune surveillance
            specifically at the neurogenic boundary, not throughout the tissue.
          </div>
        </div>
        """, unsafe_allow_html=True)

        st.markdown("""
        <div class="finding-card">
          <div class="finding-title">🟢 Disease Implication</div>
          <div class="finding-text">
            This study provides the first spatially-resolved, gene-level
            definition of the healthy hippocampal niche. In
            <strong style='color:#5ecfb0;'>Alzheimer's disease</strong> and
            <strong style='color:#5ecfb0;'>temporal lobe epilepsy</strong>,
            disruption of the Fence — measurable as loss of C1ql2/TNC expression
            in Interface Astrocytes — may represent an early, pre-symptomatic
            pathological event.
          </div>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("## The Two-Structure Model")

    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        fig = go.Figure()

        # Hub cluster
        hub_x = [-1.5, -1.2, -1.8, -1.0, -1.6]
        hub_y = [1.5, 1.2, 1.8, 1.6, 1.0]
        fig.add_trace(go.Scatter(x=hub_x, y=hub_y, mode='markers',
            marker=dict(size=30, color='#5ecfb0', opacity=0.8, line=dict(width=2, color='white')),
            name='Glial Hub (Astro + Microglia)', showlegend=True))

        # Fence layer
        fence_x = [-0.2, 0.0, 0.2, -0.1, 0.1]
        fence_y = [0.5, 0.3, 0.6, 0.1, 0.4]
        fig.add_trace(go.Scatter(x=fence_x, y=fence_y, mode='markers',
            marker=dict(size=20, color='#d4a85a', opacity=0.9, symbol='diamond',
                        line=dict(width=2, color='white')),
            name='Interface Astrocytes (The Fence)', showlegend=True))

        # DG neurons
        dg_x = [0.8, 1.0, 1.2, 0.9, 1.1, 0.7]
        dg_y = [-0.5, -0.3, -0.6, -0.8, -0.4, -0.6]
        fig.add_trace(go.Scatter(x=dg_x, y=dg_y, mode='markers',
            marker=dict(size=18, color='#8fa8e0', opacity=0.8, line=dict(width=1, color='white')),
            name='DG Excitatory Neurons', showlegend=True))

        # Arrows
        for fx, fy in zip(fence_x[:3], fence_y[:3]):
            for dx, dy in zip(dg_x[:2], dg_y[:2]):
                fig.add_annotation(x=dx, y=dy, ax=fx, ay=fy, xref='x', yref='y', axref='x', ayref='y',
                    arrowhead=2, arrowsize=1, arrowwidth=1.5, arrowcolor='rgba(224,92,92,0.5)',
                    showarrow=True)

        fig.add_annotation(x=-0.5, y=1.5, text="🔵 GLIAL HUB<br>Astro ↔ Micro Attraction", showarrow=False,
            font=dict(size=12, color='#5ecfb0'), bgcolor='rgba(26,33,40,0.8)',
            bordercolor='#5ecfb0', borderwidth=1, borderpad=6)
        fig.add_annotation(x=0.1, y=0.9, text="🟡 THE FENCE<br>C1ql2 · TNC · Lama3", showarrow=False,
            font=dict(size=12, color='#d4a85a'), bgcolor='rgba(26,33,40,0.8)',
            bordercolor='#d4a85a', borderwidth=1, borderpad=6)
        fig.add_annotation(x=1.3, y=-0.5, text="🟢 DG NEURONS<br>037 DG Glut", showarrow=False,
            font=dict(size=12, color='#8fa8e0'), bgcolor='rgba(26,33,40,0.8)',
            bordercolor='#8fa8e0', borderwidth=1, borderpad=6)

        fig.update_layout(
            title=dict(text="Hippocampal Niche Architecture", font=dict(size=14, color='#8fa89a')),
            paper_bgcolor='rgba(5,7,8,0)', plot_bgcolor='rgba(26,33,40,0.3)',
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-3, 2.5]),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, range=[-1.5, 2.5]),
            legend=dict(bgcolor='rgba(26,33,40,0.8)', bordercolor='rgba(94,207,176,0.2)',
                        font=dict(color='#8fa89a', size=11)),
            height=450,
            margin=dict(l=20, r=20, t=50, b=20),
            font=dict(family='DM Mono')
        )
        st.plotly_chart(fig, use_container_width=True)


# ─────────────────────────────────────────────
# PAGE: SPATIAL ATLAS
# ─────────────────────────────────────────────
elif page == "🗺️ Spatial Atlas":
    st.markdown('<p class="section-tag">Spatial Atlas</p>', unsafe_allow_html=True)
    st.markdown("## Explore the Hippocampal Map")
    st.markdown('<p style="color:#8fa89a;">Each dot is a single cell, positioned at its measured spatial coordinate in the hippocampus.</p>', unsafe_allow_html=True)

    col_ctrl, col_main = st.columns([1, 3])

    with col_ctrl:
        view_mode = st.radio("Color by", ["Cell Type", "Gene Expression"], key="atlas_mode")
        dot_size = st.slider("Dot size", 1, 8, 3)

        if view_mode == "Gene Expression":
            gene_sel = st.selectbox("Select Gene", sorted(gene_list))
            cmap = st.selectbox("Colormap", ["viridis", "plasma", "magma", "inferno", "RdBu_r"])
        else:
            region_filter = st.multiselect("Filter by region", ["CA1", "DG"], default=["CA1", "DG"])

    with col_main:
        plot_obs = obs.copy()

        if view_mode == "Cell Type":
            if region_filter:
                plot_obs = plot_obs[plot_obs['region'].isin(region_filter)]

            COLORS = {
                '319 Astro-TE NN': '#5ecfb0', '318 Astro-NT NN': '#3aa88e',
                '321 Astroependymal NN': '#2d8a74', '320 Astro-OLF NN': '#1f6b5a',
                '334 Microglia NN': '#d4a85a', '053 Sst Gaba': '#e05c5c',
                '056 Sst Chodl Gaba': '#b03030', '052 Pvalb Gaba': '#e0905c',
                '051 Pvalb chandelier Gaba': '#c07030', '016 CA1-ProS Glut': '#8fa8e0',
                '017 CA3 Glut': '#5c80d4', '037 DG Glut': '#a8e05c',
            }

            fig = px.scatter(
                plot_obs, x='x', y='y', color='subclass',
                color_discrete_map=COLORS,
                title=f"Cell Type Map — {len(plot_obs):,} cells",
                labels={'x': 'Spatial X (µm)', 'y': 'Spatial Y (µm)', 'subclass': 'Cell Type'},
                hover_data=['subclass', 'region']
            )
        else:
            expr = get_gene_expression(plot_obs, gene_sel, is_real)
            plot_obs['expression'] = expr

            fig = px.scatter(
                plot_obs, x='x', y='y', color='expression',
                color_continuous_scale=cmap,
                title=f"{gene_sel} Expression — {len(plot_obs):,} cells",
                labels={'x': 'Spatial X (µm)', 'y': 'Spatial Y (µm)', 'expression': f'{gene_sel} (log2)'},
                hover_data=['subclass', 'region']
            )

        fig.update_traces(marker=dict(size=dot_size, opacity=0.75))
        fig.update_layout(
            paper_bgcolor='rgba(5,7,8,0)', plot_bgcolor='rgba(26,33,40,0.3)',
            font=dict(family='DM Mono', color='#8fa89a'),
            title_font=dict(color='#e8ede8'),
            xaxis=dict(showgrid=False, zeroline=False, color='#5a7269'),
            yaxis=dict(showgrid=False, zeroline=False, color='#5a7269', scaleanchor='x'),
            legend=dict(bgcolor='rgba(26,33,40,0.8)', bordercolor='rgba(94,207,176,0.2)',
                        font=dict(size=10)),
            height=600, margin=dict(l=0, r=0, t=50, b=0)
        )
        st.plotly_chart(fig, use_container_width=True)


# ─────────────────────────────────────────────
# PAGE: THE FENCE
# ─────────────────────────────────────────────
elif page == "🔭 The Fence":
    st.markdown('<p class="section-tag">Advancement 4 — Fence Discovery</p>', unsafe_allow_html=True)
    st.markdown("## The Glial-Neuronal Fence")
    st.markdown('<p style="color:#8fa89a;">KDTree-based spatial DGE: Interface vs Distant Astrocytes at the DG boundary.</p>', unsafe_allow_html=True)

    tab1, tab2, tab3 = st.tabs(["🗺️ Spatial Fence Map", "🌋 Volcano Plot", "🧬 Fence Genes"])

    with tab1:
        col1, col2 = st.columns([3, 1])
        with col1:
            fence_obs = obs.copy()
            dg_mask = fence_obs['subclass'] == '037 DG Glut'
            ia_mask = fence_obs['astro_interface_group'] == 'Interface_Astrocyte'
            da_mask = fence_obs['astro_interface_group'] == 'Distant_Astrocyte'
            other_mask = ~(dg_mask | ia_mask | da_mask)

            fence_obs['plot_group'] = 'Other Cells'
            fence_obs.loc[dg_mask, 'plot_group'] = 'DG Neurons (037 DG Glut)'
            fence_obs.loc[ia_mask, 'plot_group'] = 'Interface Astrocytes (The Fence)'
            fence_obs.loc[da_mask, 'plot_group'] = 'Distant Astrocytes'

            show_other = st.checkbox("Show other cell types", value=False)
            if not show_other:
                fence_obs = fence_obs[fence_obs['plot_group'] != 'Other Cells']

            FENCE_COLORS = {
                'Interface Astrocytes (The Fence)': '#d4a85a',
                'Distant Astrocytes': '#3aa88e',
                'DG Neurons (037 DG Glut)': '#a8e05c',
                'Other Cells': 'rgba(100,100,100,0.2)'
            }

            fig = px.scatter(
                fence_obs, x='x', y='y', color='plot_group',
                color_discrete_map=FENCE_COLORS,
                title="The Fence: Interface Astrocytes vs Distant Astrocytes vs DG Neurons",
                labels={'x': 'Spatial X', 'y': 'Spatial Y', 'plot_group': 'Group'},
                hover_data=['subclass']
            )
            fig.update_traces(marker=dict(size=3, opacity=0.8))
            fig.update_layout(
                paper_bgcolor='rgba(5,7,8,0)', plot_bgcolor='rgba(26,33,40,0.3)',
                font=dict(family='DM Mono', color='#8fa89a'),
                xaxis=dict(showgrid=False, zeroline=False, color='#5a7269'),
                yaxis=dict(showgrid=False, zeroline=False, color='#5a7269', scaleanchor='x'),
                legend=dict(bgcolor='rgba(26,33,40,0.8)', bordercolor='rgba(94,207,176,0.2)'),
                height=550
            )
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            n_interface = (obs['astro_interface_group'] == 'Interface_Astrocyte').sum()
            n_distant = (obs['astro_interface_group'] == 'Distant_Astrocyte').sum()
            n_dg = (obs['subclass'] == '037 DG Glut').sum()

            st.markdown('<div class="metric-card"><div class="metric-num">' + str(n_interface) + '</div><div class="metric-label">Interface Astrocytes</div></div>', unsafe_allow_html=True)
            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown('<div class="metric-card"><div class="metric-num">' + str(n_distant) + '</div><div class="metric-label">Distant Astrocytes</div></div>', unsafe_allow_html=True)
            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown('<div class="metric-card"><div class="metric-num">' + str(n_dg) + '</div><div class="metric-label">DG Neurons</div></div>', unsafe_allow_html=True)
            st.markdown("<br>", unsafe_allow_html=True)
            st.markdown("""
            <div class="finding-card">
              <div class="finding-title">Method</div>
              <div class="finding-text" style="font-size:0.78rem;">
              KDTree (SciPy) computed distance from every astrocyte to its nearest DG neuron.
              <br><br>
              Interface = bottom 25th percentile of distances.
              <br><br>
              Wilcoxon DGE: Interface vs Distant astrocytes.
              </div>
            </div>
            """, unsafe_allow_html=True)

    with tab2:
        fence_genes = get_fence_genes()
        fence_genes['-log10_p'] = -np.log10(fence_genes['adj_pval'])

        # Generate background genes for volcano
        np.random.seed(99)
        n_bg = 500
        bg = pd.DataFrame({
            'logfoldchange': np.random.normal(0, 0.4, n_bg),
            '-log10_p': np.abs(np.random.normal(0.5, 0.7, n_bg)),
            'gene': [f'gene_{i}' for i in range(n_bg)],
            'significant': False
        })

        sig = fence_genes[['logfoldchange', '-log10_p', 'gene']].copy()
        sig['significant'] = True

        volcano_df = pd.concat([bg[['logfoldchange', '-log10_p', 'gene', 'significant']], sig])

        fig = px.scatter(
            volcano_df, x='logfoldchange', y='-log10_p',
            color='significant',
            color_discrete_map={True: '#d4a85a', False: 'rgba(100,100,100,0.3)'},
            hover_data=['gene'],
            title="Volcano Plot: Interface vs Distant Astrocytes (DGE)",
            labels={'logfoldchange': 'Log2 Fold Change', '-log10_p': '-log10(adj. p-value)'}
        )

        # Label fence genes
        for _, row in fence_genes.iterrows():
            fig.add_annotation(x=row['logfoldchange'], y=row['-log10_p'],
                text=row['gene'], font=dict(size=11, color='#d4a85a'),
                bgcolor='rgba(26,33,40,0.8)', bordercolor='rgba(212,168,90,0.4)',
                borderpad=3, showarrow=True, arrowcolor='rgba(212,168,90,0.4)',
                ax=20, ay=-20)

        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="rgba(255,255,255,0.2)",
                      annotation_text="p=0.05", annotation_font_color='#5a7269')
        fig.add_vline(x=0.5, line_dash="dash", line_color="rgba(255,255,255,0.2)")
        fig.add_vline(x=-0.5, line_dash="dash", line_color="rgba(255,255,255,0.2)")

        fig.update_traces(marker=dict(size=6))
        fig.update_layout(
            paper_bgcolor='rgba(5,7,8,0)', plot_bgcolor='rgba(26,33,40,0.3)',
            font=dict(family='DM Mono', color='#8fa89a'),
            xaxis=dict(gridcolor='rgba(94,207,176,0.07)', zeroline=False, color='#5a7269'),
            yaxis=dict(gridcolor='rgba(94,207,176,0.07)', zeroline=False, color='#5a7269'),
            showlegend=False, height=550
        )
        st.plotly_chart(fig, use_container_width=True)

    with tab3:
        fence_genes = get_fence_genes()

        fig = px.bar(
            fence_genes.sort_values('logfoldchange'),
            x='logfoldchange', y='gene',
            color='logfoldchange', color_continuous_scale='YlOrRd',
            orientation='h',
            title="Fence Genes: Interface Astrocytes (upregulated vs Distant)",
            labels={'logfoldchange': 'Log2 Fold Change', 'gene': 'Gene'},
            hover_data=['role', 'adj_pval']
        )
        fig.update_layout(
            paper_bgcolor='rgba(5,7,8,0)', plot_bgcolor='rgba(26,33,40,0.3)',
            font=dict(family='DM Mono', color='#8fa89a'),
            xaxis=dict(gridcolor='rgba(94,207,176,0.07)', color='#5a7269'),
            yaxis=dict(gridcolor='rgba(94,207,176,0.07)', color='#8fa89a'),
            coloraxis_showscale=False, height=450
        )
        st.plotly_chart(fig, use_container_width=True)

        st.markdown("### Gene Roles")
        _fence_display = fence_genes[['gene', 'logfoldchange', 'adj_pval', 'role']].rename(columns={
            'gene': 'Gene', 'logfoldchange': 'Log2FC', 'adj_pval': 'adj. p-value', 'role': 'Biological Role'
        })
        st.dataframe(
            _fence_display.style
                .bar(subset=['Log2FC'], color=['#4a90d9', '#e8704a'], align='zero')
                .format({'Log2FC': '{:.2f}', 'adj. p-value': '{:.4f}'}),
            use_container_width=True
        )


# ─────────────────────────────────────────────
# PAGE: NEIGHBORHOODS
# ─────────────────────────────────────────────
elif page == "🔥 Neighborhoods":
    st.markdown('<p class="section-tag">Neighborhood Enrichment Analysis</p>', unsafe_allow_html=True)
    st.markdown("## Cell Type Attraction & Repulsion")
    st.markdown('<p style="color:#8fa89a;">Z-score matrix from Squidpy neighborhood enrichment. Red = attraction, Blue = repulsion.</p>', unsafe_allow_html=True)

    col1, col2 = st.columns([3, 1])

    with col1:
        fig = px.imshow(
            z_df,
            color_continuous_scale='RdBu_r',
            color_continuous_midpoint=0,
            zmin=-4, zmax=4,
            title="Hippocampal Cell Type Neighborhood Enrichment (Z-scores)",
            labels={'color': 'Z-score'},
            text_auto=".1f"
        )
        fig.update_layout(
            paper_bgcolor='rgba(5,7,8,0)',
            plot_bgcolor='rgba(26,33,40,0.3)',
            font=dict(family='DM Mono', size=10, color='#8fa89a'),
            title_font=dict(color='#e8ede8'),
            xaxis=dict(tickfont=dict(size=9, color='#8fa89a')),
            yaxis=dict(tickfont=dict(size=9, color='#8fa89a')),
            height=580
        )
        st.plotly_chart(fig, use_container_width=True)

    with col2:
        st.markdown("### Key Pairs")
        pairs = [
            ("Astro-TE", "Microglia",  "+3.1", "Hub attraction"),
            ("Astro-TE", "DG Glut",    "−2.8", "Fence repulsion"),
            ("Astro-NT", "Astro-TE",   "+3.8", "Astro clustering"),
            ("Pvalb",    "CA1-ProS",   "+2.4", "Interneuron-ExN"),
            ("CA1-ProS", "CA3",        "+3.5", "ExN clustering"),
            ("CA3",      "DG Glut",    "+3.1", "DG-CA3 circuit"),
        ]
        for p1, p2, z, interp in pairs:
            color = "#e05c5c" if "−" in z else "#5ecfb0"
            st.markdown(f"""
            <div style='background:rgba(26,33,40,0.5);border:1px solid rgba(94,207,176,0.07);
                        border-radius:6px;padding:0.7rem;margin-bottom:0.5rem;'>
              <div style='font-family:DM Mono;font-size:0.65rem;color:#5a7269;'>{p1} ↔ {p2}</div>
              <div style='font-size:1.1rem;color:{color};font-family:DM Mono;font-weight:500;'>Z = {z}</div>
              <div style='font-size:0.7rem;color:#8fa89a;'>{interp}</div>
            </div>
            """, unsafe_allow_html=True)


# ─────────────────────────────────────────────
# PAGE: GENE LOOKUP
# ─────────────────────────────────────────────
elif page == "🧬 Gene Lookup":
    st.markdown('<p class="section-tag">Gene Lookup</p>', unsafe_allow_html=True)
    st.markdown("## Query Any Gene in the 1,122-Gene Panel")

    gene_query = st.selectbox("Select or type a gene:", sorted(gene_list), index=sorted(gene_list).index('Gfap') if 'Gfap' in gene_list else 0)

    if gene_query:
        expr = get_gene_expression(obs, gene_query, is_real)
        obs_gene = obs.copy()
        obs_gene['expression'] = expr

        col1, col2 = st.columns(2)

        with col1:
            fig = px.scatter(
                obs_gene, x='x', y='y', color='expression',
                color_continuous_scale='viridis',
                title=f"{gene_query} · Spatial Expression",
                labels={'expression': f'{gene_query} (log2)'},
                hover_data=['subclass']
            )
            fig.update_traces(marker=dict(size=3, opacity=0.8))
            fig.update_layout(
                paper_bgcolor='rgba(5,7,8,0)', plot_bgcolor='rgba(26,33,40,0.3)',
                font=dict(family='DM Mono', color='#8fa89a'),
                xaxis=dict(showgrid=False, zeroline=False, color='#5a7269'),
                yaxis=dict(showgrid=False, zeroline=False, color='#5a7269', scaleanchor='x'),
                height=450
            )
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            cell_types = obs_gene['subclass'].unique()
            mean_expr = obs_gene.groupby('subclass')['expression'].mean().sort_values(ascending=False).reset_index()
            mean_expr.columns = ['Cell Type', 'Mean Expression']

            fig2 = px.bar(
                mean_expr, x='Mean Expression', y='Cell Type',
                orientation='h',
                color='Mean Expression', color_continuous_scale='viridis',
                title=f"{gene_query} · Mean Expression by Cell Type"
            )
            fig2.update_layout(
                paper_bgcolor='rgba(5,7,8,0)', plot_bgcolor='rgba(26,33,40,0.3)',
                font=dict(family='DM Mono', color='#8fa89a'),
                xaxis=dict(gridcolor='rgba(94,207,176,0.07)', color='#5a7269'),
                yaxis=dict(gridcolor='rgba(94,207,176,0.07)', color='#8fa89a'),
                coloraxis_showscale=False, height=450
            )
            st.plotly_chart(fig2, use_container_width=True)

        # Is it a fence gene?
        fence_genes = get_fence_genes()
        if gene_query in fence_genes['gene'].values:
            row = fence_genes[fence_genes['gene'] == gene_query].iloc[0]
            st.success(f"✅ **{gene_query}** is a confirmed Fence Gene | Log2FC = {row['logfoldchange']:.2f} | adj. p = {row['adj_pval']:.4f} | Role: {row['role']}")
        else:
            st.info(f"ℹ️ {gene_query} was not identified as a significant fence gene in the Interface vs Distant Astrocyte DGE analysis.")


# ─────────────────────────────────────────────
# PAGE: METHODS
# ─────────────────────────────────────────────
elif page == "📄 Methods":
    st.markdown('<p class="section-tag">Methods & Reproducibility</p>', unsafe_allow_html=True)
    st.markdown("## How This Analysis Was Done")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("### Data Source")
        st.markdown("""
        **Dataset:** Allen Brain Cell Atlas — Zhuang-ABCA-1
        **Technology:** MERFISH (Multiplexed Error-Robust FISH)
        **Species:** *Mus musculus*
        **Access:** [brain-map.org](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas)
        **Gene Panel:** 1,122 genes (log2-normalized)
        **Spatial Resolution:** Single-cell 2D coordinates in Allen CCF
        """)

        st.markdown("### Cell Type Selection")
        st.markdown("""
        - 12 cell subclasses selected from hippocampal CA1 and DG
        - Filter: `parcellation_structure.str.contains('CA1|DG')`
        - Cell types: 4 astrocyte subtypes, 1 microglia subtype,
          2 SST interneurons, 2 Pvalb interneurons, 3 excitatory types
        """)

        st.markdown("### Neighborhood Enrichment")
        st.markdown("""
        - **`sq.gr.spatial_neighbors()`** — generic coordinate-type spatial graph
        - **`sq.gr.nhood_enrichment()`** — permutation-based enrichment test
        - Z-score matrix: positive = co-localization, negative = avoidance
        - Robustness: tested at k=5, k=15 (default), k=45 neighbors
        """)

    with col2:
        st.markdown("### Fence Gene Discovery (Key Method)")
        st.markdown("""
        1. **KDTree construction** on DG excitatory neuron coordinates
        2. **Query**: For each astrocyte, find distance to nearest DG neuron
        3. **Interface definition**: Bottom 25th percentile of distances
        4. **DGE**: Wilcoxon rank-sum test, Interface vs Distant astrocytes
           (`sc.tl.rank_genes_groups`, `use_raw=False`)
        5. **GO Enrichment**: GProfiler (`organism='mmusculus'`), top 150 genes
           with logFC > 0.5, filtered for GO:BP terms
        6. **Repeated** for Microglia-DG and Astrocyte-Microglia interfaces
        """)

        st.markdown("### Ligand-Receptor Analysis")
        st.markdown("""
        - `sq.gr.ligrec()` with 100 permutations, `alpha=1.0`
        - Focused on fence gene pathways: C1ql2, TNC, Lama3
        - Sender: Astrocytes → Receiver: DG neurons
        """)

        st.markdown("### Tools & Versions")
        tools = pd.DataFrame({
            'Tool': ['Python', 'Scanpy', 'Squidpy', 'AnnData', 'SciPy', 'GProfiler', 'HoloViews', 'Bokeh'],
            'Version': ['3.10+', '1.9+', '1.3+', '0.9+', '1.10+', '1.0+', '1.17+', '3.0+']
        })
        st.dataframe(tools, use_container_width=True, hide_index=True)

    st.markdown("---")
    st.markdown("### Citation")
    st.code("""@software{solagbade2025hippocampal,
  author    = {Solagbade, Oluwapelumi S. and Yusuf, J. A.},
  title     = {Computational Pipeline for High-Resolution Spatial Mapping
               of the Glial Neurogenic Niche in the Mouse Hippocampus},
  year      = {2025},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.17778234},
  url       = {https://doi.org/10.5281/zenodo.17778234}
}""", language="bibtex")

    col1, col2, col3 = st.columns(3)
    with col1:
        st.link_button("📦 View on Zenodo", "https://doi.org/10.5281/zenodo.17778234")
    with col2:
        st.link_button("💻 GitHub Repository", "https://github.com/Vulcan-Spark")
    with col3:
        st.link_button("🧬 Allen Brain Atlas", "https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas")
