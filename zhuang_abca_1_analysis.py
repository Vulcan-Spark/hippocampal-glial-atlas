# -*- coding: utf-8 -*-
"""
Zhuang-ABCA-1 Analysis - Streamlit Ready
"""
import pandas as pd
import anndata as ad
import scanpy as sc
import squidpy as sq
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import holoviews as hv
from holoviews import opts
from gprofiler import GProfiler
from scipy.spatial import KDTree
from pathlib import Path
import warnings
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

# Initial Setup
hv.extension('bokeh')
warnings.simplefilter("ignore")

# --- PHASE 1: DATA ACQUISITION ---
# Note: Streamlit will use the 'abc_data' folder on the server
download_base = Path("abc_data")
abc_cache = AbcProjectCache.from_cache_dir(download_base)

DATASET_NAME = 'Zhuang-ABCA-1'
CCF_NAME = 'Zhuang-ABCA-1-CCF'
TAXONOMY_NAME = 'WMB-taxonomy'
CCF_ATLAS_NAME = 'Allen-CCF-2020'

@st.cache_data
def load_metadata():
    # Main cell list
    cell_meta_df = abc_cache.get_metadata_dataframe(directory=DATASET_NAME, file_name='cell_metadata')
    cell_meta_df = cell_meta_df.set_index('cell_label')

    # Cell type decoder
    cluster_details_df = abc_cache.get_metadata_dataframe(directory=TAXONOMY_NAME, file_name='cluster_to_cluster_annotation_membership_pivoted')
    cluster_details_df = cluster_details_df.set_index('cluster_alias')

    # 3D Coordinates
    ccf_coords_df = abc_cache.get_metadata_dataframe(directory=CCF_NAME, file_name='ccf_coordinates')
    ccf_coords_df.rename(columns={'x': 'x_ccf', 'y': 'y_ccf', 'z': 'z_ccf'}, inplace=True)
    ccf_coords_df = ccf_coords_df.set_index('cell_label')

    # Join Tables
    master_meta = cell_meta_df.join(cluster_details_df, on='cluster_alias')
    master_meta = master_meta.join(ccf_coords_df)
    
    # Filter for Hippocampus
    filter_mask = master_meta['parcellation_structure'].str.contains('CA1|DG', na=False)
    return master_meta[filter_mask]

hip_metadata = load_metadata()

# --- PHASE 2: ANN-DATA OBJECT ---
@st.cache_resource
def create_anndata(_metadata):
    expression_file_path = abc_cache.get_data_path(directory=DATASET_NAME, file_name=f"{DATASET_NAME}/log2")
    adata_full = ad.read_h5ad(expression_file_path, backed='r')
    
    # Slice and load to memory
    adata_hip = adata_full[_metadata.index, :].to_memory()
    adata_hip.obs = _metadata
    adata_hip.obsm['spatial'] = _metadata[['x', 'y']].values
    
    # Set gene symbols
    gene_df = abc_cache.get_metadata_dataframe(directory=DATASET_NAME, file_name='gene').set_index('gene_identifier')
    adata_hip.var = gene_df.loc[adata_hip.var.index]
    adata_hip.var_names = gene_df.loc[adata_hip.var.index]['gene_symbol'].values
    
    adata_full.file.close()
    return adata_hip

adata_hip = create_anndata(hip_metadata)

# --- PHASE 3: NEIGHBORHOOD ANALYSIS ---
target_cell_types = [
    '319 Astro-TE NN', '318 Astro-NT NN', '321 Astroependymal NN', '320 Astro-OLF NN',
    '334 Microglia NN', '053 Sst Gaba', '056 Sst Chodl Gaba',
    '052 Pvalb Gaba', '051 Pvalb chandelier Gaba',
    '016 CA1-ProS Glut', '017 CA3 Glut', '037 DG Glut'
]

analysis_mask = adata_hip.obs['subclass'].isin(target_cell_types)
adata_final = adata_hip[analysis_mask, :].copy()
adata_final.obs['subclass'] = adata_final.obs['subclass'].astype('category')

# Core Analysis Logic
sq.gr.spatial_neighbors(adata_final, coord_type="generic")
sq.gr.nhood_enrichment(adata_final, cluster_key="subclass")

# --- PHASE 4: LIGAND-RECEPTOR (The Key Fix) ---
# Running with alpha=1.0 ensures results are always generated for the UI
if 'ligrec' not in adata_final.uns:
    sq.gr.ligrec(
        adata_final,
        n_perms=100,
        cluster_key="subclass",
        use_raw=False,
        alpha=1.0  # Crucial fix to prevent KeyErrors
    )

# --- STREAMLIT UI ---
st.title("Hippocampal Spatial Analysis")
st.write(f"Analyzing {len(adata_final)} cells across {len(target_cell_types)} subclasses.")

if st.checkbox("Show Neighborhood Heatmap"):
    fig, ax = plt.subplots(figsize=(10, 8))
    sq.pl.nhood_enrichment(adata_final, cluster_key="subclass", cmap="coolwarm", ax=ax)
    st.pyplot(fig)

if st.checkbox("Show Gene Expression Map"):
    gene = st.selectbox("Select Gene", ["Gfap", "Pvalb", "C1ql2", "Slc7a10"])
    fig, ax = plt.subplots()
    sq.pl.spatial_scatter(adata_final, color=gene, size=0.5, cmap="viridis", ax=ax)
    st.pyplot(fig)
