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

import anndata as ad

os.chdir("")

def empty_project_scp(empty_adata, scp_adata, celltype_order="celltype"):    
    sc_data = pd.DataFrame(scp_adata.X.toarray(), index=scp_adata.obs[celltype_order], columns=scp_adata.var.index)
    sc_data["celltype_simu"] = sc_data.index
    sc_data.index = range(len(sc_data))

    celltype_groups = sc_data.groupby("celltype_simu").groups    
    for key, value in celltype_groups.items():
        celltype_groups[key] = np.array(value)
    allcellname = celltype_groups.keys()
    sc_data.drop(columns=["celltype_simu"], inplace=True)
    sc_data = ad.AnnData(sc_data)
    sc_data = sc_data.X
    sc_data = np.ascontiguousarray(sc_data, dtype=np.float32)
    
    sample_count = pd.DataFrame(empty_adata.obs[celltype_order].value_counts()) # define the cell counts of each cell-type

    # sample the defined cells and integrate
    projection_value = []
    for j, cellname in enumerate(allcellname): # 遍历cell types
        select_index = np.random.choice(celltype_groups[cellname], size=int(sample_count.loc[cellname]["count"]), replace=True) 
        if (len(projection_value)  == 0):
            projection_value = sc_data[select_index]
        else:
            projection_value = np.concatenate((projection_value, sc_data[select_index]))
    print(projection_value.shape)
    projection_scp_adata = sc.AnnData(X=csr_matrix(projection_value), obs=empty_adata.obs)
    projection_scp_adata.obsm["spatial"] = empty_adata.obsm["spatial"]
    projection_scp_adata.var_names = scp_adata.var_names    

    return(projection_scp_adata)

empty_adata_for_projection = sc.read_h5ad("NSCLC_ct3_scp_empty_Lung5_R1.h5ad") # spatial tissue structure
scp_adata_for_projection = sc.read_h5ad("scp2019_1003_Reference.h5ad") # projection data
projection_scp_adata = empty_project_scp(empty_adata=empty_adata_for_projection, scp_adata=scp_adata_for_projection)

projection_scp_adata.write_h5ad("NSCLC_ct3_scp_Lung5_R1_scp2019.h5ad") # after integration, obtain synthetic NSCLC single-cell and spatially resolved proteomics data
