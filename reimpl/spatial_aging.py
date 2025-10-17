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
        geneID	x	y	MIDCounts
        1810058I24Rik	121	37	1
        2310034G01Rik	121	37	1
        4930513N10Rik	121	37	1
        4932414N04Rik	121	37	1
        AC124561.1	121	37	1
        AC154200.1	121	37	3
        AC157554.3	121	37	4
    """

    with gzip.open(file_path, 'rt') as f:
        df = pd.read_csv(f, sep='\t', header=None)

    return adata

