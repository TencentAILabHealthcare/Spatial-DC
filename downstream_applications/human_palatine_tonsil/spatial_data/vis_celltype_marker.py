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

# ------------------------------------------
import matplotlib.colors as mcolors
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/01_cell_perc")
sp_adata = sc.read_h5ad("01_data/HumanTonsil_Spatial_2492.h5ad")

sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)
sc.pp.scale(sp_adata)

colors = ['#001f3f', '#545b62', '#7fc97f', '#ffb74d', 'red']
my_cmap = mcolors.LinearSegmentedColormap.from_list('mycmap', colors)

protein_list = [
                "CD45",
                "CD19",#pan-B
                "IgD","CD23",#naive-B,activated naive-B
                "CD38",# "CD77",#GC-B
                "CD95",#"CD27","CD80",#memory-B
                "CD32",#"CD184"#Plasma-B  
                "CD3",#pan-T
                "CD4",#T4
                "CD45RA",#Naive T4
                "CD8",#
                "CD56",#NK
                "HLA_DRA","CD11c","CD123","CD21",#Itgax #DC,Myeloid DC, PDC, FDC
                "CD14","CD16",#MO/MAC
                "CD324",#Epi
               ]

sc.pl.spatial(sp_adata,color=protein_list,spot_size=1, color_map=my_cmap,vmin=-2,vmax=2, save="HumanTonsil_panel_5_known_marker_total.pdf")

