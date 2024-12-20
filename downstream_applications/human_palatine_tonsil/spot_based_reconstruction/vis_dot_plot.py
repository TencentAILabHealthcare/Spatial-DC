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
# plt.rcParams["font.family"] = "Arial"



os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/02_proteomic_profiles")
scp_adata = sc.read_h5ad("01_data/SpatialDC/SpatialDC_reconstruct_norm.h5ad")
sc.pp.filter_cells(scp_adata, min_counts=10)
sc.pp.filter_genes(scp_adata, min_cells=3)
scp_adata = scp_adata[scp_adata.obs["cellperc"] > 0.01]

sc.pp.normalize_total(scp_adata)
sc.pp.log1p(scp_adata)
sc.pp.scale(scp_adata)


ct_order = ['NBC','Activated NBC','GCBC','MBC', 'PC', 
              'CD4 T', 'Naive CD4 T','CD8 T','Naive CD8 T', 'cycling T', 'DN', 'NK', 
              'DC','FDC','PDC',  'Mono/Macro',  'Granulocytes', 
       'ILC',  'epithelial']

scp_adata.obs["celltype"] = scp_adata.obs["celltype"].cat.set_categories(ct_order,ordered=False)

import matplotlib.colors as mcolors

colors = ['#001f3f', '#545b62', '#7fc97f', '#ffb74d', 'red']
my_cmap = mcolors.LinearSegmentedColormap.from_list('mycmap', colors)

markers_list = ["CD19",#B"
    
    "CD3","CD4","CD8", # T cells
    "CD56","CD16",#NK
    "CD11c","CD123",#cDC,pDC
    "CD11b",
    "CD324",
]
ax=sc.pl.dotplot(scp_adata, markers_list,return_fig=True, groupby='celltype', swap_axes=True,standard_scale="var", cmap=my_cmap)
# ax.show()
ax.savefig("figures/selected_markers.pdf")


