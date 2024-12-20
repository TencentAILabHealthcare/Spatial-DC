# -------------- 信号通讯分析-------------------
# COMMOT
import sys
import matplotlib.pyplot as plt
import seaborn as sns

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

# Commot
import gc
import ot
import pickle
from scipy import sparse
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt

import commot as ct

os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/01_cell_perc")
scp_adata = sc.read_h5ad("03_benchmark_methods/exp_ref46/SpatialDC/SpatialDC_reconstruct_not_norm.h5ad")

# set zero of cellperc<0
set_zero = scp_adata.obs["cellperc"] < 0.01
temp_arr = np.array(scp_adata.to_df())
temp_arr[set_zero,:]=0
scp_adata.X = temp_arr

sender_cell = "FDC"
sender_pid = "CD257" # TNFSF13B
receiver_cell = "GCBC"
receiver_pid = "CD267" #TNFRSF13B

keep_index = scp_adata.var["name2"].isin([sender_pid,receiver_pid])
keep_ct = scp_adata.obs["celltype"].isin([sender_cell])
target_adata = scp_adata[keep_ct,keep_index].copy()

temp_index = [f"spot_{i}" for i in range(target_adata.shape[0])]
target_adata.obs.index = temp_index

temp_df = target_adata.to_df()

keep_index = scp_adata.var["name2"].isin([receiver_pid])
keep_ct = scp_adata.obs["celltype"].isin([receiver_cell])
temp_df[receiver_pid]  = scp_adata[keep_ct,keep_index].to_df().values
target_adata.X = np.array(temp_df)

# prepare adata
target_adata.obs[["Y"]] = -target_adata.obs[["Y"]]
target_adata.obsm["spatial"] = np.array(target_adata.obs[["X","Y"]])
sc.pp.log1p(target_adata)

# scaler = MinMaxScaler(feature_range=(0, 1))
# inputs = scaler.fit_transform(target_adata.to_df())
# target_adata.X = inputs

# 准备PPI
target_ppi_df = pd.DataFrame({0:sender_pid,1:receiver_pid,2:'sig_ppi',3:'Secreted Signaling'},index=[0,1,2,3])
ct.tl.spatial_communication(target_adata,
    database_name='selected_db', df_ligrec=target_ppi_df, dis_thr=500, heteromeric=True, pathway_sum=True) 

# 展示信号传递的方向
ct.tl.communication_direction(target_adata, database_name='selected_db', pathway_name='sig_ppi', k=5)

ct.pl.plot_cell_communication(target_adata, database_name='selected_db', pathway_name='sig_ppi',normalize_v=True,normalize_v_quantile=0.995,summary="sender",
grid_density=0.4,scale=0.5,ndsize=0,grid_thresh=10,grid_scale=10,filename=f"{sender_cell}_{receiver_cell}_sender_0.5.pdf",plot_method="cell")

ct.pl.plot_cell_communication(target_adata, database_name='selected_db', pathway_name='sig_ppi',normalize_v=True,normalize_v_quantile=0.995,clustering=None,summary="receiver",
grid_density=0.4,scale=0.5,ndsize=0,grid_thresh=10,grid_scale=10,filename=f"{sender_cell}_{receiver_cell}_receiver_0.5.pdf",plot_method="cell")


# # 展示信号传递的强度
# adata = target_adata.copy()
# pts = adata.obsm['spatial']
# s = adata.obsm['commot-selected_db-sum-sender'][f's-{sender_pid}-{receiver_pid}']
# r = adata.obsm['commot-selected_db-sum-receiver'][f'r-{sender_pid}-{receiver_pid}']
# fig, ax = plt.subplots(1,2, figsize=(5,2))
# ax[0].scatter(pts[:,0], pts[:,1], c=s, s=1, cmap='coolwarm')
# ax[0].set_title(f'Sender_{sender_pid}-{receiver_pid}')
# ax[1].scatter(pts[:,0], pts[:,1], c=r, s=1, cmap='coolwarm')
# ax[1].set_title(f'Receiver_{sender_pid}-{receiver_pid}')

target_adata.obs[["Y"]] = -target_adata.obs[["Y"]]
target_adata.obsm["spatial"] = np.array(target_adata.obs[["X","Y"]])
target_adata.obs[f's-{sender_pid}-{receiver_pid}'] = target_adata.obsm['commot-selected_db-sum-sender'][f's-{sender_pid}-{receiver_pid}']
target_adata.obs[f'r-{sender_pid}-{receiver_pid}'] = target_adata.obsm['commot-selected_db-sum-receiver'][f'r-{sender_pid}-{receiver_pid}']
sc.pl.spatial(target_adata,spot_size=1,color=[f's-{sender_pid}-{receiver_pid}',f'r-{sender_pid}-{receiver_pid}'],cmap="coolwarm",save=f"{sender_cell}_{receiver_cell}_strength.pdf")

# plt.savefig(f"figures/Sender_Receiver_{sender_pid}-{receiver_pid}.pdf")
