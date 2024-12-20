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

# for NSCLC reference single-cell proteomics data [scp2021]
data_df = pd.read_excel("scp2021.xlsx",index_col=0) # LFQ Quantified
meta_df = pd.read_excel("scp2021_meta.xlsx", index_col=0)

data_df.fillna(0, inplace=True)
reference_adata = sc.AnnData(X=csr_matrix(np.array(data_df.T)), obs=meta_df)
reference_adata.obs.index = reference_adata.obs.SampleID
reference_adata.var_names = data_df.index

reference_adata.obs["celltype"] = reference_adata.obs.Group
reference_adata.obs.celltype = reference_adata.obs.celltype.astype('category')

sc.pp.filter_genes(reference_adata, min_cells=1) # 108 * 1437

# for NSCLC projection single-cell proteomics data [scp2019]
data_df = pd.read_csv("scp2019_1225.csv",sep="\t",index_col=0)
projection_adata.X = csr_matrix(data_df.T)

# retain commonly identified proteins
intersect_id = np.intersect1d(projection_adata.var_names,reference_adata.var_names)

# save file
projection_adata[:,intersect_id].write_h5ad("scp2019_1003_Projection.h5ad")
reference_adata[:,intersect_id].write_h5ad("scp2021_1003_Reference.h5ad")
