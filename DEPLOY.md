# DEPLOYMENT GUIDE
# Hippocampal Glial Niche Atlas
# ===============================

## STEP 1: Create the GitHub Repository

1. Go to github.com → click "+" → "New repository"
2. Name it: hippocampal-glial-atlas
3. Set to Public
4. Click "Create repository"

---

## STEP 2: Upload These Files to GitHub

Upload in this exact folder structure:

hippocampal-glial-atlas/
├── README.md                          ← Upload this
├── app.py                             ← Upload this
├── requirements.txt                   ← Upload this
├── zhuang_abca_1_analysis.py          ← Upload your analysis script
└── data/
    └── .gitkeep                       ← Create this empty file

To create the data/ folder on GitHub:
  - Click "Add file" → "Create new file"
  - Type: data/.gitkeep
  - Commit it

---

## STEP 3: Upload Your h5ad File (Optional but Recommended)

Your healthy_blueprint_data.h5ad file is likely too large for GitHub (>25MB limit).
Upload it to Zenodo instead (you already have an account):

1. Go to zenodo.org → New Upload
2. Upload healthy_blueprint_data.h5ad
3. Set as "Open Access"
4. Add it to your existing DOI record OR create a new one
5. Add the download link in your README

The app runs in "Demo Mode" without the h5ad file, showing
simulated data that reflects your real findings.

---

## STEP 4: Deploy on Streamlit Cloud (Free)

1. Go to: share.streamlit.io
2. Sign in with GitHub
3. Click "New app"
4. Settings:
   - Repository: Vulcan-Spark/hippocampal-glial-atlas
   - Branch: main
   - Main file path: app.py
5. Click "Deploy!"

Your app will be live at:
https://vulcan-spark-hippocampal-glial-atlas.streamlit.app

(Or a similar URL that Streamlit assigns)

---

## STEP 5: Link Everything Together

In your website (index.html), add the atlas link to your research card:

<a class="rc-link" href="https://vulcan-spark-hippocampal-glial-atlas.streamlit.app" target="_blank">↗ Live Atlas</a>

In your Zenodo record description, add:
"Interactive atlas available at: [Streamlit URL]"

---

## TROUBLESHOOTING

Problem: App fails to load on Streamlit Cloud
Fix: Check requirements.txt — remove anndata if not loading h5ad

Problem: scipy import error
Fix: Add this to requirements.txt: scipy>=1.10.0

Problem: App is slow
Fix: The demo mode uses n=5000 cells for fast rendering.
     With real h5ad, add st.cache_data to all loading functions (already done).

---

## LOADING REAL DATA (When Ready)

Place your file at: data/healthy_blueprint_data.h5ad

The app will automatically detect it and load real data.
All 65,956 cells will be available for exploration.

If the file is too large for Streamlit Cloud memory (>2GB):
1. Pre-compute a subset AnnData with only CA1/DG cells (already done — that IS your file)
2. Or use Streamlit's file_uploader to let users upload it themselves
