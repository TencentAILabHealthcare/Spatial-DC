# Visualization of spatial proteomics data
# 20241119
# Yuan

import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mat
import os
import sys
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# --------------------------------------
os.chdir("xxx")
sp_adata = sc.read_h5ad("xxx.h5ad")

# sc.pp.filter_cells(sc_adata, min_genes=10)
# sc.pp.filter_genes(sc_adata, min_cells=3)

sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)

sc.pp.neighbors(sp_adata)
sc.tl.umap(sp_adata)
sc.tl.leiden(sp_adata,resolution=0.1)

sc.pl.umap(sp_adata, color=['leiden'],size=100)
sc.pl.spatial(sp_adata, spot_size=1, color='leiden')

