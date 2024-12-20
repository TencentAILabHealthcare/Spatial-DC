import sys
from Spatial_DC_beta import SpatialDC
import scanpy as sc
import  numpy as np
import os

import time
start_time = time.time()


def train_SpatialDC(sc_adata, sp_adata, celltype_key,output_file_path):
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)           

    if not os.path.exists(model_path):        
        spatial_dc = SpatialDC(sc_adata=sc_adata, sp_adata=sp_adata, print_info=True) 
        spatial_dc.setup_distribution_model(spot_num=10000, epochs=200, batch_size=128, lr=0.001) # Then, setup distribution_model, users can change default params such as learning rate or epochs
        spatial_dc.train_distribution_model() # Begin train this model, which need some times
        spatial_dc.save_distribution_model(save_model_path = model_path)

# -----------------------------------------------------------------------
# Train the distriubtion model for each reference datasets
os.chdir("")
datasets = ["MousePDAC", "NSCLC", "HumanTonsil", "MouseBrain"]

method = "SpatialDC_V1"
celltype_key = "celltype"

for dataset in datasets:
    print(f"============Dataset: {dataset}==============")
    if dataset == "NSCLC":
        sc_file_path = f"01_data/{dataset}/scp2021_1003_Reference.h5ad"        
        sp_file_path = f"01_data/{dataset}/Simu_seed0_cells10_noise0.h5ad"
    elif dataset == "HumanTonsil":
        sc_file_path = f"01_data/{dataset}/HumanTonsil_reference_ct19_46_intersected.h5ad"        
        sp_file_path = f"01_data/{dataset}/HumanTonsil_Spatial_2492_intersected.h5ad"
    elif dataset == "MouseBrain":
        sc_file_path = f"01_data/{dataset}/MouseBrain2022_ct4_4351.h5ad"        
        sp_file_path = f"01_data/{dataset}/MouseBrain2022_spot208_4351.h5ad"
    elif dataset == "MousePDAC":
        sc_file_path = f"01_data/{dataset}/MousePDAC2023_ct10_2837.h5ad"        
        sp_file_path = f"01_data/{dataset}/MousePDAC2023_spot108_2837.h5ad"

    output_dir = f"03_benchmark_methods/{dataset}/{method}/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)   

    model_dir = f"{output_dir}/{method}_model/"
    model_path = f"{output_dir}/{method}_model/model_epochs200.pt"

    output_file_path = f"{output_dir}/{method}"            

    sc_adata = sc.read_h5ad(sc_file_path)
    sp_adata = sc.read_h5ad(sp_file_path)

    intersect = np.intersect1d(sc_adata.var_names, sp_adata.var_names)        
    sp_adata = sp_adata[:, intersect].copy()
    sc_adata = sc_adata[:, intersect].copy()     

    sc.pp.normalize_total(sc_adata)
    sc.pp.normalize_total(sp_adata)
    run_SpatialDC(sc_adata=sc_adata, sp_adata=sp_adata, celltype_key=celltype_key,output_file_path=output_file_path)

end_time = time.time()
print("Total time consumption: seconds")
print(end_time - start_time)
