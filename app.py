import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# 1. FORCE REFRESH CSS & ANIMATIONS
st.set_page_config(page_title="Glial Niche Atlas v2", layout="wide")

st.markdown("""
<style>
    /* Ultra-Modern Dark Theme & Animation */
    @keyframes bgAnimation {
        0% { background-position: 0% 50%; }
        50% { background-position: 100% 50%; }
        100% { background-position: 0% 50%; }
    }

    .stApp {
        background: linear-gradient(-45deg, #05070a, #0c121e, #05070a, #111827) !important;
        background-size: 400% 400% !important;
        animation: bgAnimation 15s ease infinite !important;
        color: #e8ede8 !important;
    }

    /* Glassmorphism for Sidebar and Cards */
    [data-testid="stSidebar"] {
        background-color: rgba(10, 15, 25, 0.8) !important;
        backdrop-filter: blur(15px);
    }

    .metric-card {
        background: rgba(255, 255, 255, 0.03);
        border: 1px solid rgba(255, 255, 255, 0.1);
        border-radius: 15px;
        padding: 20px;
        backdrop-filter: blur(10px);
        box-shadow: 0 8px 32px 0 rgba(0, 0, 0, 0.8);
    }

    /* Glowing Text */
    .glow-text {
        color: #2dd4bf;
        text-shadow: 0 0 10px rgba(45, 212, 191, 0.5);
        font-family: 'Courier New', monospace;
    }
</style>
""", unsafe_allow_html=True)

# 2. DESIGNER PALETTE
PALETTE = {
    '319 Astro-TE NN': '#2dd4bf', '318 Astro-NT NN': '#0d9488',
    '321 Astroependymal NN': '#0891b2', '320 Astro-OLF NN': '#0369a1',
    '334 Microglia NN': '#fbbf24', '053 Sst Gaba': '#818cf8',
    '016 CA1-ProS Glut': '#fb7185', '037 DG Glut': '#f43f5e'
}

# 3. DATA SIMULATION (Without caching to ensure it updates)
def get_data():
    np.random.seed(42)
    data = []
    for cell, color in PALETTE.items():
        n = 200
        x = np.random.normal(np.random.randint(-5, 5), 1, n)
        y = np.random.normal(np.random.randint(-5, 5), 1, n)
        for i in range(n):
            data.append({'x': x[i], 'y': y[i], 'subclass': cell})
    return pd.DataFrame(data)

df = get_data()

# 4. UI LAYOUT
st.markdown('<h1 class="glow-text">HIPPOCAMPAL GLIAL ATLAS</h1>', unsafe_allow_html=True)

col1, col2 = st.columns([2, 1])

with col1:
    fig = px.scatter(df, x='x', y='y', color='subclass', 
                     color_discrete_map=PALETTE, template="plotly_dark")
    fig.update_layout(paper_bgcolor='rgba(0,0,0,0)', plot_bgcolor='rgba(0,0,0,0)')
    st.plotly_chart(fig, use_container_width=True)

with col2:
    st.markdown('<div class="metric-card">', unsafe_allow_html=True)
    st.write("### Niche Stats")
    st.metric("Total Cells", "65,956", "+12%")
    st.markdown('</div>', unsafe_allow_html=True)

    # 5. THE MATPLOTLIB FIX (Crucial part)
    st.write("### Gene Enrichment")
    genes = pd.DataFrame({
        'gene': ['C1ql2', 'Tnc', 'Lama3'],
        'lfc': [4.2, 3.1, 2.5]
    })

    try:
        # This requires matplotlib to work
        st.dataframe(genes.style.background_gradient(cmap='viridis', subset=['lfc']))
    except Exception:
        # This is the fallback so your app DOES NOT CRASH if matplotlib is missing
        st.warning("Install matplotlib for color-coded tables.")
        st.dataframe(genes)

st.info("💡 If you still don't see colors, clear your browser cache or run with `streamlit run app.py --server.clearCache`.")
