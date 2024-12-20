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
import os

plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ------------------------------------------
from sklearn.cluster import KMeans

# Panel-2, cluster of spatial proteomics
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/01_cell_perc")

sp_adata = sc.read_h5ad("01_data/HumanTonsil_Spatial_2492.h5ad")

sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)

kmeans = KMeans(n_clusters=3, random_state=42).fit(np.array(sp_adata.to_df()))
sp_adata.obs["true_labels"] = kmeans.labels_.astype(str)

sc.pp.neighbors(sp_adata)
sc.tl.umap(sp_adata)

sc.tl.leiden(sp_adata,resolution=0.2)
sc.pl.umap(sp_adata, color=["true_labels",'leiden'],size=10, save="HumanTonsil_spatial_domain_clustering_umap.pdf")
sc.pl.spatial(sp_adata, spot_size=1, color=["true_labels"],title=["kmeans"],palette = ['#ff7f0e','#1f77b4', '#2ca02c'],save="HumanTonsil_spatial_domain_clustering_spatial_map_kmeans.pdf") #final decide using this
sc.pl.spatial(sp_adata, spot_size=1, color=["leiden"],title=["leiden"],palette = ['#ff7f0e','#2ca02c','#1f77b4'],save="HumanTonsil_spatial_domain_clustering_spatial_map_leiden.pdf")
