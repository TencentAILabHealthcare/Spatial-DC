import sys
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
import pandas as pd
import scanpy as sc
import numpy as np
import warnings
warnings.filterwarnings("ignore")

import cell2location
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

import time

# ---------------------------------------------
start_time = time.time()

def train_Cell2location(sc_adata,sp_adata, celltype_key,sc_epochs):

    ## train sc-model
    cell2location.models.RegressionModel.setup_anndata(adata=sc_adata,labels_key=celltype_key)
    mod = cell2location.models.RegressionModel(sc_adata)

    # Use all data for training (validation not implemented yet, train_size=1)
    mod.train(max_epochs=sc_epochs, use_gpu=True)
    mod.save(model_dir)

    mod.history["elbo_train"].plot()
    plt.title("Cell2location_sc_train")
    plt.savefig(model_path + '.png')
    plt.close()
        

def run_Cell2location(sc_adata,sp_adata, output_file_path, celltype_key,sc_epochs,sp_epochs):

    sc_adata.X = round(sc_adata.X)
    sp_adata.X = round(sp_adata.X)

    intersect = np.intersect1d(sc_adata.var_names, sp_adata.var_names)        
    sc_adata = sc_adata[:, intersect].copy()
    sp_adata = sp_adata[:, intersect].copy()

    if not os.path.exists(model_dir):
        # os.makedirs(model_dir)    # not need to create model_dir
        train_Cell2location(sc_adata, sp_adata, celltype_key,sc_epochs=sc_epochs)


    cell2location.models.RegressionModel.setup_anndata(adata=sc_adata,labels_key=celltype_key)
    mod = cell2location.models.RegressionModel(sc_adata)
    mod.load(model_dir,adata=sc_adata)

    sc_adata = mod.export_posterior(
        sc_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
        # sc_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2, 'use_gpu': True}
    )
    
    # train sp
    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in sc_adata.varm.keys():
        inf_aver = sc_adata.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in sc_adata.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = sc_adata.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in sc_adata.uns['mod']['factor_names']]].copy()
    inf_aver.columns = sc_adata.uns['mod']['factor_names']

    # train sp    
    intersect = np.intersect1d(sp_adata.var_names, inf_aver.index)
    sp_adata = sp_adata[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=sp_adata)

    # create and train the model
    sp_mod = cell2location.models.Cell2location(
        sp_adata, cell_state_df=inf_aver#,
    )

    sp_mod.train(max_epochs=sp_epochs,
            # train using full data (batch_size=None)
            batch_size=None,
            # use all data points in training because
            # we need to estimate cell abundance at all locations
            train_size=1,
            use_gpu=True)



    sp_adata = sp_mod.export_posterior(
        sp_adata, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
        # sp_adata, sample_kwargs={'num_samples': 1000, 'batch_size': 2, 'use_gpu': True}
    )
    
    cell2loc_results = sp_adata.obsm['q05_cell_abundance_w_sf']
    cell2loc_results = (cell2loc_results.T/cell2loc_results.sum(axis=1)).T

    # change columns 
    celltype2=[]
    for i in cell2loc_results.columns:
        temp = i.replace("q05cell_abundance_w_sf_","")
        celltype2.append(temp)    
    cell2loc_results.columns = celltype2

    cell2loc_results.to_csv(output_file_path+".csv")
    
os.chdir("")

# cell-type data
reference_data_dir = "01_data/reference/reference_variable_proteins"
celltype_key = "celltype"
sp_data_dir = "01_data/simulations/synthetic_cellnum_noise_v4"

# output dir
method = "Cell2location"
output_dir = f"03_output/exp_conditions_v4/{method}/"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)   

for cond in ["top100","top200","top400","top600","top800","tail600","tail800","tail100","tail200","tail400"]: 

    model_dir = f"{output_dir}/{method}model_raw_raw_Reference_{cond}/"
    model_path = f"{model_dir}/model.pt"            

    for seed in [0]:
        for cells in [10,15,20]:
            for noise in [0]:   
                                            
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

                run_Cell2location(sc_adata=sc_adata, sp_adata=sp_adata,
                        output_file_path=output_file_path,sp_epochs=5000,sc_epochs=10000,
                        celltype_key=celltype_key)            

end_time = time.time()
print("Total time consumption: seconds")
print(end_time - start_time)
