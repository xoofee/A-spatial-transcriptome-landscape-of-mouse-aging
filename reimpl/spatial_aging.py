"""
ipython.run_line_magic('matplotlib', 'tk') if (ipython := get_ipython()) else None


Stereo-seq Spot Configuration: (web search by cursor)
Spot diameter: ~220 nanometers
Center-to-center distance: 500 nanometers between spots
This means there's about 280nm spacing between adjacent spots (500nm - 220nm = 280nm gap)

For 2-month-old young adult mice:
Cell Length: ~100-150 micrometers (μm)
Cell Width: ~15-20 micrometers (μm)

"""

#%% Spatial Transcriptomics Analysis Pipeline
# Implementation of "A spatial transcriptome landscape of mouse aging" (Cell 2024)
# Following the roadmap: Phase 1 - Data Preprocessing & Quality Control

import pandas as pd
import numpy as np
import scipy.sparse as sp
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import gzip
import warnings
warnings.filterwarnings('ignore')
from anndata import AnnData

def load_gem_file(file_path: str) -> AnnData:

    # hallmark of aging\STT0000039\Analysis\STSA0000367\STTS0000720\Y-2M-Heart-1_1_filtered_bin50.gem\Y-2M-Heart-1_1_filtered_bin50.gem

    """ 
    # male C57BL/6J mice 

        geneID	x	y	MIDCounts
        1810058I24Rik	121	37	1
        2310034G01Rik	121	37	1
        4930513N10Rik	121	37	1
        4932414N04Rik	121	37	1
        AC124561.1	121	37	1
        AC154200.1	121	37	3
        AC157554.3	121	37	4
    """

    # file_path = r"D:\work\bio\spatial\hallmark of aging\STT0000039\Analysis\STSA0000367\STTS0000720\Y-2M-Heart-1_1_filtered_bin50.gem.gz"

    with gzip.open(file_path, 'rt') as f:
        gem = pd.read_csv(f, sep='\t')

    # adata = AnnData(X=df,  obsm={"spatial": df[['x', 'y']].values})

    agg = gem.groupby(["x", "y", "geneID"])["MIDCounts"].sum().reset_index()

    # Pivot
    expr_matrix = agg.pivot_table(
        index=["x", "y"], columns="geneID", values="MIDCounts", fill_value=0
    )

    # Create AnnData
    adata = sc.AnnData(X=expr_matrix.values)
    adata.obs = pd.DataFrame(index=[f"{x}_{y}" for x, y in expr_matrix.index])
    adata.var = pd.DataFrame(index=expr_matrix.columns)
    adata.obsm["spatial"] = np.array(expr_matrix.index.tolist())

    return adata


# Now usable in Squidpy
import squidpy as sq

# adata = load_gem_file(r"D:\work\bio\spatial\hallmark of aging\STT0000039\Analysis\STSA0000367\STTS0000720\Y-2M-Heart-1_1_filtered_bin50.gem.gz")
# adata = load_gem_file(r"D:\work\bio\spatial\hallmark of aging\STT0000039\Analysis\STSA0000368\STTS0000721\Y-2M-Heart-2_1_filtered_bin50.gem.gz")
adata = load_gem_file(r"D:\work\bio\spatial\hallmark of aging\STT0000039\Analysis\STSA0000400\STTS0000753\Y-2M-Liver-2_2_filtered_bin50.gem.gz")


#%%
adata.obs["total_counts"] = adata.X.sum(axis=1)
adata.obs['n_genes'] = (adata.X > 0).sum(1)

sq.pl.spatial_scatter(adata, shape=None, size=1, color="total_counts")


#%% aggregation test


n = 5  # number of rows to aggregate
rows, cols = adata.X.shape

# Make sure number of rows is divisible by n
trimmed_rows = (rows // n) * n
X_trimmed = adata.X[:trimmed_rows]  # drop extra rows if not divisible by n

# Reshape and sum
adata.X_agg = X_trimmed.reshape(-1, n, cols).sum(axis=1)

new_genes = (adata.X_agg > 0 ).sum(1)

plt.hist(adata.obs['n_genes'], bins=100)


