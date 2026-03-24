# 🧠 Hippocampal Glial Niche Atlas
### *A Computational Blueprint of the Healthy Hippocampus*

> **"Before we can understand how the brain breaks, we must first precisely define how it is built."**

[![Python](https://img.shields.io/badge/Python-3.10+-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![Scanpy](https://img.shields.io/badge/Scanpy-1.9+-brightgreen?style=flat-square)](https://scanpy.readthedocs.io/)
[![Squidpy](https://img.shields.io/badge/Squidpy-1.3+-orange?style=flat-square)](https://squidpy.readthedocs.io/)
[![Data](https://img.shields.io/badge/Data-Allen%20Brain%20Cell%20Atlas-blue?style=flat-square)](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.17778234-teal?style=flat-square)](https://doi.org/10.5281/zenodo.17778234)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](LICENSE)
[![Presented](https://img.shields.io/badge/Presented-AfInMic%202025-purple?style=flat-square)](https://doi.org/10.5281/zenodo.17778234)

---

## The Scientific Question

The hippocampus is one of the only regions in the adult mammalian brain where neurogenesis continues throughout life — new neurons are born in the **Dentate Gyrus (DG)** and must integrate into existing circuits. This process depends on the precise spatial organization of glial cells (astrocytes, microglia) that act as both **gatekeepers and supporters** of this neurogenic niche.

**But what does "healthy" actually look like at single-cell resolution?**

Without a quantitative spatial baseline, we cannot detect the *earliest* glial disruptions that precede Alzheimer's disease, epilepsy, or age-related cognitive decline. This project provides that baseline.

---

## What Was Discovered

Using **MERFISH spatial transcriptomics** on 65,956 hippocampal cells from the Allen Brain Cell Atlas (Zhuang-ABCA-1 dataset), this study mapped the spatial organization of 12 major cell populations across CA1 and the Dentate Gyrus — and uncovered a previously uncharacterized structural principle:

### 🔵 The Glial-Neuronal Fence

At the boundary between **DG excitatory neurons** (`037 DG Glut`) and the surrounding **astrocyte population**, there exists a molecularly distinct layer of **Interface Astrocytes** — the 25% of astrocytes spatially closest to the DG neurogenic zone.

These Interface Astrocytes are not simply passive bystanders. They express a distinct genetic program, including:

| Gene | Role | Direction |
|------|------|-----------|
| **C1ql2** | Complement-related synaptic organizer; boundary-defining | ⬆ Upregulated |
| **Tnc** (Tenascin-C) | Extracellular matrix repulsion molecule | ⬆ Upregulated |
| **Lama3** | Laminin alpha-3; basement membrane anchoring | ⬆ Upregulated |
| **Slc7a10** | Neutral amino acid transporter; metabolic gating | ⬆ Upregulated |
| **Thbs4** | Thrombospondin-4; synaptogenesis inhibitor | ⬆ Upregulated |

This constitutes what we term the **"Glial-Neuronal Fence"** — an active molecular barrier that likely regulates which new neurons successfully integrate into the hippocampal circuit.

### 🔴 The Glial Hub

Simultaneously, astrocytes and microglia exhibit **strong mutual attraction** (neighborhood enrichment Z-score > +2), suggesting co-localized surveillance of the neurogenic niche. This astrocyte-microglia cluster forms the "hub" of the niche, providing neurotrophic support.

### The Two-Structure Model

```
  ┌──────────────────────────────────────────────────┐
  │           HIPPOCAMPAL NEUROGENIC NICHE           │
  │                                                  │
  │  [GLIAL HUB]               [FENCE]               │
  │  Astro ←→ Microglia   →→→  Interface Astro       │
  │  (Attraction, Z>+2)         (C1ql2, TNC, Lama3)  │
  │                              ↓↓↓                  │
  │                          DG Neurons               │
  │                        (037 DG Glut)              │
  │                                                   │
  └───────────────────────────────────────────────────┘
```

> **Disease Hypothesis:** Alzheimer's disease and temporal lobe epilepsy may initiate specifically when the Fence is disrupted — when Interface Astrocytes lose their boundary-defining gene signature and begin expressing reactive astrogliosis markers instead. This provides the first quantitative baseline to test that hypothesis.

---

## Dataset

| Property | Value |
|----------|-------|
| **Dataset** | Zhuang-ABCA-1 (Allen Brain Cell Atlas) |
| **Species** | *Mus musculus* |
| **Technology** | MERFISH (Multiplexed Error-Robust FISH) |
| **Gene Panel** | 1,122 genes |
| **Total Cells** | 4,000,000+ brain-wide |
| **This Study** | **65,956 hippocampal cells** (CA1 + DG) |
| **Brain Regions** | Cornu Ammonis 1 (CA1), Dentate Gyrus (DG) |
| **Cell Types** | 12 subclasses across 4 major classes |
| **Data Access** | [Allen Brain Cell Atlas Portal](https://portal.brain-map.org/atlases-and-data/bkp/abc-atlas) |

---

## Analysis Pipeline

```
Phase 1: Setup (abc_atlas_access, scanpy, squidpy)
    ↓
Phase 2: Data Acquisition
         → cell_metadata, cluster_to_cluster_annotation,
           ccf_coordinates, parcellation terms, gene table
    ↓
Phase 3: Master Metadata Table
         → 4-way join: cells + types + 3D coords + anatomy
    ↓
Phase 4: Hippocampal Filtering
         → CA1 | DG → 65,956 cells
    ↓
Phase 5: AnnData Construction
         → 65,956 cells × 1,122 genes + spatial coords
    ↓
Phase 7a: Neighborhood Enrichment
          → sq.gr.spatial_neighbors → sq.gr.nhood_enrichment
          → Z-score matrix: attraction/repulsion between all 12 types
    ↓
Phase 7b: CA1 vs DG Regional Comparison
          → Separate neighborhood analyses per subregion
    ↓
Phase 7c/d: Spatial Expression Maps
            → Pvalb, Gfap mapped onto 2D tissue coordinates
    ↓
Advancement 4: Fence Gene Discovery (KDTree DGE)
               → Interface Astrocytes (25th percentile proximity)
               → Wilcoxon DGE vs Distant Astrocytes
               → GO Enrichment (gProfiler) → ECM, boundary, adhesion
    ↓
Advancement 4b: Microglial Fence Genes
                → Same KDTree approach for Microglia-DG boundary
                → Upregulated + Downregulated gene sets
    ↓
Advancement 6: Astrocyte-Microglia Attraction Genes
               → Interface Astrocytes nearest to Microglia
               → GO: inflammatory regulation, phagocytosis signaling
    ↓
Advancement 7: Ligand-Receptor Analysis
               → sq.gr.ligrec (100 permutations)
               → C1ql2, TNC, Lama3 pathways: Astro → DG Neurons
    ↓
Robustness Checks: k=5, k=15 (default), k=45 neighbor variation
```

---

## Reproduction

### Requirements

```bash
pip install scanpy squidpy gprofiler-official holoviews bokeh seaborn
pip install "abc_atlas_access[notebooks] @ git+https://github.com/alleninstitute/abc_atlas_access.git"
```

### Environment

```
Python       3.10+
scanpy       1.9+
squidpy      1.3+
anndata      0.9+
scipy        1.10+
seaborn      0.12+
gprofiler    1.0+
```

### Run

The full analysis is in a single Google Colab notebook:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1hpBLpGIY5xauIu6pjCYd9v0AhBcZpsZ3)

Or run the standalone script:
```bash
python zhuang_abca_1_analysis.py
```

> ⚠️ **Note:** The full gene expression file (~3GB) is downloaded automatically from the Allen Brain Cell Atlas on first run. Subsequent runs use a cached version. Ensure ~5GB free disk space.

---

## Repository Structure

```
hippocampal-glial-atlas/
├── zhuang_abca_1_analysis.py     # Main analysis script (all phases)
├── app/
│   ├── app.py                    # Streamlit atlas web app
│   ├── requirements.txt          # App dependencies
│   └── data/
│       └── healthy_blueprint_data.h5ad  # Processed AnnData (download separately)
├── figures/
│   ├── spatial_map_pvalb_gfap.png        # Gene expression spatial map
│   ├── neighborhood_heatmap.png          # Cell type attraction/repulsion
│   ├── volcano_interface_astrocytes.png  # Fence gene DGE
│   ├── go_enrichment_fence_genes.png     # GO biological processes
│   └── chord_diagram.html                # Interactive network graph
├── results/
│   ├── fence_genes_astrocyte_dg.csv       # Full DGE results table
│   ├── fence_genes_microglia_dg.csv       # Microglial DGE results
│   └── zscore_matrix.csv                  # Full neighborhood Z-score matrix
└── README.md
```

---

## Key Results Summary

| Analysis | Finding | Z-score / p-value |
|----------|---------|-------------------|
| Astrocyte ↔ Microglia neighborhood | **Strong mutual attraction** (Glial Hub) | Z > +2.0 |
| Astrocyte ↔ DG Neuron neighborhood | **Active repulsion** (The Fence) | Z < −2.0 |
| Interface vs Distant Astrocytes (DGE) | **C1ql2, TNC, Lama3, Thbs4** upregulated | adj. p < 0.05 |
| Interface Microglia (DGE) | Immune surveillance genes upregulated | adj. p < 0.05 |
| GO Enrichment (Fence Genes) | Extracellular matrix, cell adhesion, boundary formation | p < 0.01 |
| Robustness (k=5, 15, 45) | Fence pattern stable across neighbor parameters | — |
| Ligand-Receptor | C1ql2, TNC pathways: Astro → DG signaling | p < 0.05 |

---

## Citation

If you use this analysis, data, or code, please cite:

```bibtex
@software{solagbade2025hippocampal,
  author    = {Solagbade, Oluwapelumi S. and Yusuf, J. A.},
  title     = {Computational Pipeline for High-Resolution Spatial Mapping
               of the Glial Neurogenic Niche in the Mouse Hippocampus
               (MERFISH Analysis)},
  year      = {2025},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.17778234},
  url       = {https://doi.org/10.5281/zenodo.17778234},
  note      = {Presented at AfInMic 2025 International Student Microscopy
               Symposium, 3rd Place Oral Presentation}
}
```

**Original dataset:**
> Zhuang, X. et al. (2023). A molecularly defined and spatially resolved
> cell atlas of the whole mouse brain. *Nature*, 624, 343–354.
> https://doi.org/10.1038/s41586-023-06808-9

---

## About

This work was conducted by **Oluwapelumi Solagbade** (MBChB Candidate, Obafemi Awolowo University, Nigeria) as an independent computational study in the Eagle Research Laboratory, Ladoke Akintola University of Technology (LAUTECH), Nigeria.

This project demonstrates that high-impact spatial transcriptomics research can be conducted in resource-limited settings using entirely open-access data and infrastructure — a deliberate proof-of-concept for African computational neuroscience.

**Contact:** solagbadepelumi@gmail.com | [ORCID](https://orcid.org/0009-0005-7748-356X) | [GitHub](https://github.com/Vulcan-Spark)

---

*Built in Nigeria. With open data. For the world.*
