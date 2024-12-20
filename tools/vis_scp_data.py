# visualization of single-cell proteomics data
# Yuan
# 20241119

import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr,spearmanr
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import matplotlib as mat
import os
import sys
from scipy import stats
import warnings
from sklearn.preprocessing import MinMaxScaler,StandardScaler
warnings.filterwarnings("ignore")

plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ---------------------------------
os.chdir("")

sc_adata = sc.read_h5ad("xxx.h5ad")
# sc_adata = sc_adata[sc_adata.obs.celltype!= "cycling myeloid"].copy()

sc.pp.filter_cells(sc_adata, min_genes=10)
sc.pp.filter_genes(sc_adata, min_cells=3)

# with annotated cell types
ct_order = ['xxx','xxxx']
sc_adata.obs["celltype"] = sc_adata.obs["celltype"].astype("category")
sc_adata.obs["celltype"] = sc_adata.obs["celltype"].cat.set_categories(ct_order,ordered=False)

my_palette = ["#66C2A5", "#E41A1C"]
sc.pl.umap(sc_adata, color=["celltype"], size=10, palette=my_palette) 

# with un-supervised clustering
sc.pp.normalize_total(sc_adata)
sc.pp.log1p(sc_adata)

sc.pp.neighbors(sc_adata)
sc.tl.umap(sc_adata)

sc.tl.leiden(sc_adata,resolution=0.2)
sc.pl.umap(sc_adata, color=['leiden'],size=10)
