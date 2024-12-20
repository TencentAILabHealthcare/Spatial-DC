# Methods for scale during pre-process stage
# Yuan
# 20241119

import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import os
import sys
from scipy import stats
from sklearn.preprocessing import scale, MinMaxScaler,StandardScaler
import warnings
warnings.filterwarnings("ignore")

os.chdir("xxx")

# -----------------
# MinMaxScaler
# scaled by feature for data.frame
scaler = MinMaxScaler(feature_range=(0, 1))

sp_adata = sc.read_h5ad("xxx.h5ad")
sp_adata.X = scaler.fit_transform(sp_adata.to_df()) 

pred_df = pd.read_csv("xxx.csv",index_col=0)
pred_df = pd.DataFrame(scaler.fit_transform(pred_df), index=pred_df.index, columns=pred_df.columns)

# ------------------
# StandardScale
# scaled by feature for adata
sp_adata = sc.read_h5ad("xxx.h5ad")
sc.pp.scale(sp_adata)

# StandardScale
# scaled by feature for matrix
inputs = scale(inputs,axis=0) # 0, by feature; 1, by sample
