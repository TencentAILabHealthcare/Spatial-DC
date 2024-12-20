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


# Panel-1, umap of single-cell proteomics 
# 2024-07-30
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/01_cell_perc")

sc_adata = sc.read_h5ad("../00_raw/20220215_tonsil_atlas_cite_seurat_obj.h5ad")
sc_adata.obsm["X_umap"] = np.array(sc_adata.obs[["UMAP_1_level_1","UMAP_2_level_1"]])

sc_adata = sc_adata[sc_adata.obs.celltype!= "cycling myeloid"].copy()
sc_adata = sc_adata[sc_adata.obs.celltype!= "cycling FDC"].copy()
sc_adata = sc_adata[sc_adata.obs.celltype!= "Mast"].copy()
sc_adata = sc_adata[sc_adata.obs.celltype!= "preB/T"].copy()

sc_adata = sc_adata[sc_adata.obs.subproject == "BCLLATLAS_46"].copy()

sc.pp.filter_cells(sc_adata, min_genes=10)
sc.pp.filter_genes(sc_adata, min_cells=3)

ct_order = ['Activated NBC','NBC','GCBC','MBC', 'PC', # 5 PC, Plasma cells; B
              'cycling T', 'DN','ILC','NK',#NK   #4 DN, double-negative T cells with a profile of proinflamatory activation
              'CD4 T', 'Naive CD4 T', #2
              'CD8 T','Naive CD8 T',  #2                            
              'DC','PDC', 'FDC',  #3 PDC, plasmactyoid DC; DC
              'Mono/Macro','Granulocytes',  #2 MO
              #GN
              'epithelial' # 1 Epi and ILC
       ]
sc_adata.obs["celltype"] = sc_adata.obs["celltype"].astype("category")
sc_adata.obs["celltype"] = sc_adata.obs["celltype"].cat.set_categories(ct_order,ordered=False)
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
sc_adata.shape
sc.pl.umap(sc_adata, color=["celltype"], size=10,palette=my_palette, save="HumanTonsil_panel_1_scp_ref46.pdf")
