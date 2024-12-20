# Methods for set varying colors for categories
# Yuan
# 20241119

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

import matplotlib.colors as mcolors

plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

os.chdir("")
sc_adata = sc.read_h5ad("xxx.h5ad")

# with categorical colors
# If provided, values of adata.uns["{var}_colors"] will be set.
ct_order = ['xxx','xxxx']
sc_adata.obs["celltype"] = sc_adata.obs["celltype"].astype("category")
sc_adata.obs["celltype"] = sc_adata.obs["celltype"].cat.set_categories(ct_order,ordered=False)
my_palette = ["#66C2A5", "#84CEB7"] 
sc.pl.umap(sc_adata, color=["celltype"], size=10,palette=my_palette)

# with continous colors
# Can be a name or a Colormap instance (e.g. "magma”, "viridis" or mpl.cm.cividis)
colors = ['#001f3f', '#545b62', '#7fc97f', '#ffb74d', 'red']
my_cmap = mcolors.LinearSegmentedColormap.from_list('mycmap', colors)

sc.pl.spatial(sp_adata,color=["CD19","CD3","CD45","CD324"],spot_size=1, color_map=my_cmap,vmin=-2,vmax=2)
