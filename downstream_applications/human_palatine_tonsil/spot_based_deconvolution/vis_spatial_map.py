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

# Panel-6, The distribution of total cells
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/01_cell_perc")

methods = ["SpatialDC","Tangram","Cell2location","Stereoscope","DestVI","Seurat","spotlight","spatialdwls","CARD","RCTD","MuSiC"]

for m in methods:
    sp_adata = sc.read_h5ad("01_data/HumanTonsil_Spatial_2492.h5ad")

    if (m == "SpatialDC"):
        pred_df = pd.read_csv(f"03_benchmark_methods/exp_ref46/{m}/reconstruct.csv",index_col=0)
    else:
        pred_df = pd.read_csv(f"03_benchmark_methods/exp_ref46/{m}/{m}.csv",index_col=0)
    
    scaler = MinMaxScaler(feature_range=(0, 1))
    
    
    pred_df["B"] = pred_df[['Activated NBC','NBC','GCBC','MBC', 'PC']].sum(axis=1).values        
    pred_df["T"] = pred_df[['CD4 T', 'CD8 T','Naive CD4 T','Naive CD8 T', 'cycling T', 'DN']].sum(axis=1).values
    pred_df["Leukocytes"] = pred_df[['B','T','NK','DC','PDC', 'FDC', 'Mono/Macro',  'Granulocytes','ILC']].sum(axis=1).values
    pred_df["Epi"] = pred_df[['epithelial']].values
    
    ct_order = ['Activated NBC','NBC','GCBC','MBC', 'PC', # Plasma cells
                 'CD4 T', 'CD8 T','Naive CD4 T','Naive CD8 T', 'cycling T', 'DN', 'NK',  #DN, double-negative T cells with a profile of proinflamatory activation
                 'DC','PDC', 'FDC', 'Mono/Macro',  'Granulocytes', # plasmactyoid DC
          'ILC',  'epithelial',"B","T","Leukocytes"]
    
    pred_df = pred_df[ct_order]
    ct_order_rename = ['Activated NBCs','NBCs','GCBCs','MBCs', 'PCs', # Plasma cells
                 'T4', 'T8','Naive T4','Naive T8', 'Cycling T cells', 'DN', 'NK',  #DN, double-negative T cells with a profile of proinflamatory activation
                 'DCs','PDCs', 'FDCs', 'MO',  'GN', # plasmactyoid DC
          'ILC',  'Epi',"B cells","T cells","Leukocytes"]
    
    pred_df = pd.DataFrame(scaler.fit_transform(pred_df), index=pred_df.index, columns=pred_df.columns)
    
    sp_adata.obs = pred_df
    
    sc.pl.spatial(sp_adata, color=ct_order, spot_size=1, color_map="magma", save=f"HumanTonsil_spatial_map_{m}.pdf", title=ct_order_rename)
