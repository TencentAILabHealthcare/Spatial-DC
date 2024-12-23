# Load necessary packages
import os
import sys
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from scipy.sparse import csr_matrix
from scipy.stats import pearsonr,spearmanr
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import matplotlib as mat
from scipy import stats
from sklearn.preprocessing import MinMaxScaler,StandardScaler
import warnings
warnings.filterwarnings("ignore")

import time
start_time = time.time()

# set the working directory
os.chdir("./")
sys.path.append("Spatial_DC")

from Spatial_DC import SpatialDC
SpatialDC.get_SpatialDC_version() # 1.0.0

dataset_dir = "datasets"
model_dir = "trained_model"
output_dir = "output"

datasets = ["NSCLC", "human_palatine_tonsil", "mouse_brain_coronal", "mouse_PDAC"]

# Begin train
celltype_key = "celltype"
reference_data_type = ""
dataset_type = ""

for dataset in datasets:
    print(f"Dataset: {dataset}")
    # Set the file path and datasets type
    if dataset == "NSCLC":
        sc_file_path = f"{dataset_dir}/{dataset}/synthetic_noise_levels/reference_proteomics_noise/reference_noise0.h5ad"        
        sp_file_path = f"{dataset_dir}/{dataset}/synthetic_noise_levels/spatial_proteomics_spotsize_noise/spatial_spotsize200_noise0.h5ad"
        dataset_type = "synthetic"
    else:
        sc_file_path = f"{dataset_dir}/{dataset}/intersected_reference_proteomics.h5ad"        
        sp_file_path = f"{dataset_dir}/{dataset}/intersected_spatial_proteomics.h5ad"
        dataset_type = "real"

    # Set the reference data type
    if dataset in ["human_palatine_tonsil", "NSCLC"]:
        reference_data_type = "single_cell"
        
    elif dataset in ["mouse_brain_coronal", "mouse_PDAC"]:
        reference_data_type = "single_cell_type"
    else:
        raise ValueError('Invalid dataset name')
        

    sc_adata = sc.read_h5ad(sc_file_path)
    sp_adata = sc.read_h5ad(sp_file_path)

    intersect = np.intersect1d(sc_adata.var_names, sp_adata.var_names)        
    sp_adata = sp_adata[:, intersect].copy()
    sc_adata = sc_adata[:, intersect].copy()     

    sc.pp.normalize_total(sc_adata)
    sc.pp.normalize_total(sp_adata)
    model_dir = f"output/{dataset}/model/"
    model_path = f"{model_dir}/trained_model.pt"

    # Construct the SpatialDC object    
    spatial_dc = SpatialDC(sc_adata=sc_adata, sp_adata=sp_adata, celltype_key=celltype_key, reference_data_type=reference_data_type, dataset_type=dataset_type) 
    spatial_dc.load_distribution_model(load_model_path = model_path)    

    spatial_dc.transfer_distribution_model() # inference stage

    # Set the parameters and train the self-supvervised model
    spatial_dc.setup_self_supervised_model(k_graph=30, epochs=200)
    spatial_dc.train_self_supervised_model()
    spatial_dc.inference_self_supervised_model()

    # Output-1, cell-type composition
    pred_ct = spatial_dc.sp_adata.uns["cellperc_refined"] # Obtain the predicted cell-type composition

    # Output-2, reconstructed single-cell data
    rec_scp_adata = spatial_dc.reconstruction(norm=True)    
    
end_time = time.time()
print("Total time consumption: seconds")
print(end_time - start_time)
