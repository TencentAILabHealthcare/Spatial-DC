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
import os

plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# Panel-3, Pie plot for All Benchmark methods
# using Python Script
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/01_cell_perc")
project_dir = "03_benchmark_methods/exp_ref46"
sp_adata = sc.read_h5ad("01_data/HumanTonsil_Spatial_2492.h5ad")

ct_order = ['Activated NBC','NBC','GCBC','MBC', 'PC', # 5 PC, Plasma cells; B
              'cycling T', 'DN','ILC','NK',#NK   #4 DN, double-negative T cells with a profile of proinflamatory activation
              'CD4 T', 'Naive CD4 T', #2
              'CD8 T','Naive CD8 T',  #2                            
              'DC','PDC', 'FDC',  #3 PDC, plasmactyoid DC; DC
              'Mono/Macro','Granulocytes',  #2 MO
              'epithelial' # 1 Epi and ILC
       ]


my_palette = [
    "#66C2A5", "#84CEB7", "#A3DAC9", "#C2E6DB", "#E1F3ED", # B
    "#FC8D62", "#FCA481", 
    "#377EB8",#ILC
     "#4DAF4A",  #NK 
    "#8DA0CB", "#A4B3D5", # T4
    "#984EA3", "#AC71B5", # T8
    "#A6D854", "#B7DF76", "#C9E799",#DC
    "#FFD92F", "#E5C494", #MO and GN
    "#E41A1C", #Epi
    ]


method_list = ["SpatialDC"]#["Tangram","SpatialDC"]#"Cell2location","Stereoscope","DestVI","SpatialDC_Pseudo",
spatial = pd.DataFrame(sp_adata.obsm['spatial'], index=sp_adata.obs_names,columns=["x","y"])


for method in method_list:

    if(method =="SpatialDC_initial"):
        temp_pred_path = f"{project_dir}/SpatialDC/SpatialDC_initial.csv"
    elif(method =="SpatialDC"):
        temp_pred_path = f"{project_dir}/SpatialDC/SpatialDC.csv"

    else:
        temp_pred_path = f"{project_dir}/{method}/{method}.csv"
        
    pred_df = pd.read_csv(temp_pred_path,index_col=0)
    pred_df = pred_df[ct_order].clip(0)
    
    fig, axs = plt.subplots(50, 50, figsize=(10, 10))
    plt.subplots_adjust(hspace=-0.2, wspace=-0.2)  # 调整这里的值来缩小间距  

    for i in range(50):  
        for j in range(50): 
            axs[i, j].set_xticks([])
            axs[i, j].set_yticks([])
            axs[i, j].axis('off')
            
    for spot_num in range(spatial.shape[0]):
        i = spatial["y"][spot_num] - 1
        j = spatial["x"][spot_num] - 1
        
        
        axs[i, j].pie(np.array(pred_df)[spot_num].clip(0), colors=my_palette,radius=0.01, wedgeprops=dict(width=1, edgecolor='none'))
        axs[i, j].set_aspect('equal', 'box')
        axs[i, j].set_xticks([])
        axs[i, j].set_yticks([])
    
    
    plt.savefig(f"figures/HumanTonsil_panel_3_pie_{method}.pdf")    

    plt.show()
