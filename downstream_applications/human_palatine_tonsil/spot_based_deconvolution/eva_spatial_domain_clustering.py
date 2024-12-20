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

methods=["SpatialDC_initial","SpatialDC_reconstruct","Tangram","Cell2location","Stereoscope","DestVI","CARD","spotlight","RCTD","spatialdwls","Seurat"]
save_ARI = []  
save_Purity = []  

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

    ari = adjusted_rand_score(sp_adata.obs['true_labels'].values, sp_adata.obs[f"{m}_labels"])    
   
    
    # 计算Purity值
    y_pred = sp_adata.obs[f"{m}_labels"].values.astype(int)
    y_true = sp_adata.obs['true_labels'].values.astype(int)
    N = len(y_pred)
    purity = 0
    for k in np.unique(y_pred):
        idx = (y_pred == k)
        labels = y_true[idx]
        count = np.bincount(labels)
        purity += np.max(count)
    purity /= N

    print(f'{m}_{ari}_{purity}')
    save_ARI.append(ari)
    save_Purity.append(purity)

# plot
sc.pl.spatial(sp_adata, color=sp_adata.obs.columns[6:],size=1,spot_size=1,palette=['#ff7f0e','#1f77b4', '#2ca02c'],save="comparision_kmean_cluster.pdf")    

# save 
eva_df = pd.DataFrame({  
    'ARI': save_ARI,  
    'Purity': save_Purity
}, index=methods)  

eva_df.to_csv("figures/eva_kmean_cluster.csv")

# plot
import matplotlib.pyplot as plt
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/01_cell_perc")
x=["SpatialDC","Tangram","Cell2location","Stereoscope","DestVI","Seurat","spotlight","spatialdwls","CARD"]
# 数值数据
y = [0.649, 0.103, 0.209, 0.572, 0.428,
    0.122,0.02,0.112,0.063]

fig = plt.figure(figsize=(6, 3))
# 绘制bar plot
plt.bar(x, y,width=0.7)
for i in range(len(x)):
    plt.text(x=i, y=y[i]+0.05, s=y[i], ha='center')
# 添加标题和标签
plt.title('ARI')
plt.xticks(rotation=45)
plt.ylim(0,0.9)
plt.xticks([])
# 显示图形
# plt.show()
plt.savefig("figures/benchmark/cellperc_corr_ARI.pdf")
