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

os.chdir("")

# vis of cell-type distribution and construction of empty adata
# downloaded from [https://nanostring.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/], "SMI_Giotto_Object.RData" from "All+SMI+Giotto+object.tar.gz" file was utilzed to obtain sp_meta_df and sp_loc_df, which record the cell location of primary cells
sp_meta_df = pd.read_csv("sp_meta.csv", index_col=0)
sp_loc_df = pd.read_csv("spatial_loca_raw.csv", index_col=0)

sp_meta_df[["x","y"]] = sp_loc_df[["sdimx","sdimy"]].values * 100
sp_meta_df["celltype"] = sp_meta_df["cell_type"]

empty_adata = sc.AnnData(X=csr_matrix(np.ones([sp_meta_df.shape[0],10])), obs=sp_meta_df)

# begin filter
empty_adata = empty_adata[empty_adata.obs["patient"] == "Lung5"]
empty_adata = empty_adata[empty_adata.obs["Run_Tissue_name"] == "Lung5_Rep1"]
empty_adata = empty_adata[empty_adata.obs["celltype"].isin(["macrophage","tumor 5","endothelial"])]

# sort and order
empty_adata.obs["celltype"][empty_adata.obs["celltype"] == "tumor 5"] = "C10"#"Group1" #C10
empty_adata.obs["celltype"][empty_adata.obs["celltype"] == "macrophage"] = "RAW"#"Group2" #RAW
empty_adata.obs["celltype"][empty_adata.obs["celltype"] == "endothelial"] = "SVEC"#"Group3" # SVEC
empty_adata.obs["celltype"] = empty_adata.obs["celltype"].astype("category")

empty_adata.obs.sort_values('celltype', ascending=True, inplace=True)
empty_adata.obsm["spatial"] = np.array(empty_adata.obs[['x', 'y']])

empty_adata.write_h5ad("NSCLC_ct3_scp_empty_Lung5_R1.h5ad")
