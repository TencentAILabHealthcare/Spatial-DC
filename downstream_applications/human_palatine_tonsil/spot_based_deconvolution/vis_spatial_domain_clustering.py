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

from sklearn.metrics import adjusted_rand_score
from scipy.spatial.distance import jensenshannon

from sklearn.metrics import silhouette_score
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
# Panel- 
# 基于细胞比例进行聚类，看看与基于count 聚类的区别
# panel-4 umap of predicted cell-type perc
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/01_cell_perc")

ct_order = ['Activated NBC','NBC','GCBC','MBC', 'PC', # 5 PC, Plasma cells; B
              'cycling T', 'DN','ILC','NK',#NK   #4 DN, double-negative T cells with a profile of proinflamatory activation
              'CD4 T', 'Naive CD4 T', #2
              'CD8 T','Naive CD8 T',  #2                            
              'DC','PDC', 'FDC',  #3 PDC, plasmactyoid DC; DC
              'Mono/Macro','Granulocytes',  #2 MO
              #GN
              'epithelial' # 1 Epi and ILC
       ]

sp_adata = sc.read_h5ad("01_data/HumanTonsil_Spatial_2492.h5ad")
sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)
# sc.pp.neighbors(sp_adata)
# sc.tl.umap(sp_adata)

# sc.tl.leiden(sp_adata,resolution=0.2)
kmeans = KMeans(n_clusters=3, random_state=42).fit(np.array(sp_adata.to_df()))
sp_adata.obs["true_labels"] = kmeans.labels_.astype(str)
categories_mapping = {'0':'Follicle', '1':'Inter-Follicle','2':'Crypt'}
sp_adata.obs[f"true_labels"] = sp_adata.obs[f"true_labels"].replace(categories_mapping)
sp_adata.uns[f"true_labels_colors"] = ['#6B705C','#B5838D','#FFE8D6']

methods=["SpatialDC_initial","SpatialDC_reconstruct","Tangram","Cell2location","Stereoscope","DestVI","CARD","spotlight","RCTD","spatialdwls","Seurat"]

for m in methods:
    if(m == "SpatialDC_initial"):
        ct_perc_pred = pd.read_csv("03_benchmark_methods/exp_ref46/SpatialDC/initial.csv",index_col=0)
    elif(m == "SpatialDC_reconstruct"):
        ct_perc_pred = pd.read_csv("03_benchmark_methods/exp_ref46/SpatialDC/reconstruct.csv",index_col=0)
    else:
        ct_perc_pred = pd.read_csv(f"03_benchmark_methods/exp_ref46/{m}/{m}.csv",index_col=0)

    ct_perc_adata = ad.AnnData(X=csr_matrix(np.array(ct_perc_pred[ct_order]).clip(0)))    
    ct_perc_adata.obsm["spatial"] = sp_adata.obsm["spatial"]    
    sc.pp.filter_genes(ct_perc_adata, min_cells=1)
    kmeans = KMeans(n_clusters=3, random_state=42).fit(np.array(ct_perc_pred[ct_order]).clip(0))

    sp_adata.obs[f"{m}_labels"] = kmeans.labels_.astype(str)

    if(m in ["SpatialDC_initial"]):

        categories_mapping = {'0':'Crypt','1':'Follicle','2':'Inter-Follicle'}
        sp_adata.obs[f"{m}_labels"] = sp_adata.obs[f"{m}_labels"].replace(categories_mapping)
        sp_adata.uns[f"{m}_labels_colors"] = ['#6B705C','#B5838D','#FFE8D6']
    elif(m in ["SpatialDC_reconstruct","Stereoscope","DestVI"]):
        
        categories_mapping = {'0':'Follicle','1':'Inter-Follicle','2':'Crypt'}
        sp_adata.obs[f"{m}_labels"] = sp_adata.obs[f"{m}_labels"].replace(categories_mapping)
        sp_adata.uns[f"{m}_labels_colors"] = ['#6B705C','#B5838D','#FFE8D6']
    elif(m in ["Tangram"]):
        
        categories_mapping = {'0':'Crypt','1':'Inter-Follicle','2':'Follicle'}
        sp_adata.obs[f"{m}_labels"] = sp_adata.obs[f"{m}_labels"].replace(categories_mapping)
        sp_adata.uns[f"{m}_labels_colors"] = ['#6B705C','#B5838D','#FFE8D6']
    elif(m in ["Cell2location","CARD","spatialdwls","spotlight","RCTD"]):
        
        categories_mapping = {'0':'Inter-Follicle','1':'Follicle','2':'Crypt'}
        sp_adata.obs[f"{m}_labels"] = sp_adata.obs[f"{m}_labels"].replace(categories_mapping)
        sp_adata.uns[f"{m}_labels_colors"] = ['#6B705C','#B5838D','#FFE8D6']
    else:
        categories_mapping = {'0':'Follicle','1':'Crypt','2':'Inter-Follicle'}
        sp_adata.obs[f"{m}_labels"] = sp_adata.obs[f"{m}_labels"].replace(categories_mapping)
        sp_adata.uns[f"{m}_labels_colors"] = ['#6B705C','#B5838D','#FFE8D6']

# plot
sc.pl.spatial(sp_adata, color=sp_adata.obs.columns[6:],size=1,spot_size=1, save="comparision_kmean_cluster_20240925.pdf")    
  
# sc.pl.spatial(sp_adata, color=sp_adata.obs.columns[6:],size=1,spot_size=1,palette=['#ff7f0e','#1f77b4', '#2ca02c'],save="comparision_kmean_cluster.pdf")    
