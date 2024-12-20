
import sys
from Spatial_DC import SpatialDC

import scanpy as sc
import  numpy as np
import os

import time
start_time = time.time()

def run_SpatialDC(sc_adata, sp_adata, celltype_key,output_file_path):
    if not os.path.exists(model_path):        
        exit() # should set correct path of model
   
    spatial_dc = SpatialDC(sc_adata=sc_adata, sp_adata=sp_adata) 
    spatial_dc.load_distribution_model(load_model_path = model_path) # Load trained DIS model

    # predict intermediate results with transfer model
    pred_sp_adata = spatial_dc.transfer_distribution_model()
    pred_sp_adata.uns["cellperc_initial"].to_csv(output_file_path + "_initial.csv")
    pred_sp_adata.write_h5ad(output_file_path + "_initial.h5ad")

    purified_sp_adata = spatial_dc.initial_purify_spots(norm=True, spatial_info=False)
    purified_sp_adata.write_h5ad(output_file_path + "_initial_norm.h5ad")

    purified_sp_adata = spatial_dc.initial_purify_spots(norm=False, spatial_info=False)
    purified_sp_adata.write_h5ad(output_file_path + "_initial_not_norm.h5ad")

    # refined intermediate results with self-supervisd model
    spatial_dc.setup_reconstruction_model(k_graph=30, epochs=200, w_cls=25, w_rec=25) # With trained DIS model, users can transfer it with register spatial proteomics data.    
    reconstruct_sp_adata = spatial_dc.reconstruct()
    reconstruct_sp_adata.write_h5ad(output_file_path + "_reconstruct.h5ad")

    purified_sp_adata = spatial_dc.purify_spots(norm=True, spatial_info=True)
    purified_sp_adata.write_h5ad(output_file_path + "_reconstruct_norm.h5ad")

    purified_sp_adata = spatial_dc.purify_spots(norm=False, spatial_info=True)
    purified_sp_adata.write_h5ad(output_file_path + "_reconstruct_not_norm.h5ad")

    reconstruct_sp_adata.uns["cellperc_reconstruct"].to_csv(output_file_path + "_reconstruct.csv")


# ------------------------------------------------------------
# deconvolution of synthetic, public and in-house spatial proteomics data
os.chdir("")
datasets = ["NSCLC", "HumanTonsil", "MouseBrain", "MousePDAC"]

method = "SpatialDC"
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
