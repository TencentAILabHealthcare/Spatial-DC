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

import scvi
from scvi.model import CondSCVI, DestVI
os.environ["CUDA_VISIBLE_DEVICES"] = "0"

import time

# ----------------------------------------------------------
start_time = time.time()

def train_DestVI(sc_adata,sp_adata, celltype_key,sc_epochs):

    CondSCVI.setup_anndata(sc_adata, labels_key=celltype_key)
    sc_model = CondSCVI(sc_adata, weight_obs=False)
    sc_model.train(max_epochs=sc_epochs,use_gpu=True)
    sc_model.save(model_dir)

    plt.plot(sc_model.history["elbo_train"], label="train")
    plt.title("DestVI_sc_train")
    plt.legend()
    plt.savefig(model_path + '.png')
    plt.close()

def run_DestVI(sc_adata, sp_adata, celltype_key, output_file_path, sc_epochs, sp_epochs):  
    
    ## train DestVI
    sc_adata.X = round(sc_adata.X)
    sp_adata.X = round(sp_adata.X)
    # intersect = np.intersect1d(sc_adata.var_names, sp_adata.var_names)        
    # sc_adata = sc_adata[:, intersect].copy()
    # sp_adata = sp_adata[:, intersect].copy() # need intersect

    if not os.path.exists(model_dir):
        # os.makedirs(model_dir)    # not need to create model_dir
        train_DestVI(sc_adata, sp_adata, celltype_key,sc_epochs=sc_epochs)

    CondSCVI.setup_anndata(sc_adata, labels_key=celltype_key)
    sc_model = CondSCVI.load(model_dir, adata=sc_adata)
    
    DestVI.setup_anndata(sp_adata)

    # train sp-model
    spatial_model = DestVI.from_rna_model(sp_adata, sc_model)
    spatial_model.train(max_epochs=sp_epochs,use_gpu=True)

    spatial_model.get_proportions().to_csv(output_file_path+".csv")
 
os.chdir("")

# cell-type data
reference_data_dir = "01_data/reference/reference_variable_proteins"
celltype_key = "celltype"
sp_data_dir = "01_data/simulations/synthetic_cellnum_noise_v4"

# output dir
method = "DestVI"
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

                run_DestVI(sc_adata=sc_adata, sp_adata=sp_adata,
                        output_file_path=output_file_path,sc_epochs=2000,sp_epochs=2000,
                        celltype_key=celltype_key)            

end_time = time.time()
print("Total time consumption: seconds")
print(end_time - start_time)
