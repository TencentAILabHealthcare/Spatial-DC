import pandas as pd
import scanpy as sc
import numpy as np
import tangram as tg
import torch
import os
device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

import time

# -----------------------------------------------------------------------
start_time = time.time()

def run_Tangram(sc_adata, sp_adata, celltype_key,cond, output_file_path):

    sc.pp.log1p(sc_adata)
    sc.tl.rank_genes_groups(sc_adata, groupby=celltype_key, use_raw=False)

    markers_df = pd.DataFrame(sc_adata.uns["rank_genes_groups"]["names"]).iloc[0:2000, :]

    genes_sc = np.unique(markers_df.melt().value.values)
    genes_st = sp_adata.var_names.values
    genes = list(set(genes_sc).intersection(set(genes_st)))
    tg.pp_adatas(sc_adata, sp_adata, genes=genes)

    ad_map = tg.map_cells_to_space(
                    sc_adata,
                    sp_adata,
                    mode='clusters',
                    cluster_label=celltype_key,device=device)

    tg.project_cell_annotations(ad_map, sp_adata, annotation=celltype_key)

    celltype_density = sp_adata.obsm['tangram_ct_pred']
    celltype_density = (celltype_density.T/celltype_density.sum(axis=1)).T

    celltype_density.to_csv(f"{output_file_path}.csv")    

os.chdir("")

# cell-type data
reference_data_dir = "01_data/reference/reference_variable_proteins"
celltype_key = "celltype"
sp_data_dir = "01_data/simulations/synthetic_cellnum_noise_v4"

# output dir
method = "Tangram"
output_dir = f"03_output/exp_conditions_v4/{method}/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)                

# Run    
for cond in ["top100","top200","top400","top600","top800","tail600","tail800","tail100","tail200","tail400"]: 
# for cond in ["top600","top800","top1003","tail600","tail800"]: 
    for seed in [0]:
        for cells in [10,15,20]:
            for noise in [0]: # without noise
                        
                print(f"seed:{seed}_cells:{cells}_noise{int(noise*100)}")
                            
                sc_file_path = f"{reference_data_dir}/scp2021_1003_Reference_{cond}.h5ad"
                sp_file_path = f"{sp_data_dir}/Simu_seed{seed}_cells{str(cells)}_noise{int(noise*100)}.h5ad"
                                
                sc_adata = sc.read_h5ad(sc_file_path)
                sp_adata = sc.read_h5ad(sp_file_path)

                intersect = np.intersect1d(sc_adata.var_names, sp_adata.var_names)        
                sp_adata = sp_adata[:, intersect].copy()
                sc_adata = sc_adata[:, intersect].copy()

                sc.pp.normalize_total(sc_adata)
                sc.pp.normalize_total(sp_adata)

                output_file_path = f"{output_dir}/{method}_reference_{cond}_seed{seed}_cells{str(cells)}_noise{int(noise*100)}"
 
                run_Tangram(sc_adata=sc_adata, sp_adata=sp_adata, celltype_key=celltype_key,cond=cond,output_file_path=output_file_path)

end_time = time.time()
print("Total time consumption: seconds")
print(end_time - start_time)
