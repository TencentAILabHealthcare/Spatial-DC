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


# 绘制不同细胞之间的共定位特征
import cell2location
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/01_cell_perc")

methods = ["SpatialDC"]

sp_adata = sc.read_h5ad("01_data/HumanTonsil_Spatial_2492.h5ad")
pred_df = pd.read_csv(f"03_benchmark_methods/exp_ref46/SpatialDC/reconstruct.csv",index_col=0)
scaler = MinMaxScaler(feature_range=(0, 1))

pred_df = pd.DataFrame(scaler.fit_transform(pred_df), index=pred_df.index, columns=pred_df.columns)
sp_adata.obs = pred_df

sender_cell = ["DC","CD4 T", "FDC"]
to_cell = ["Naive CD4 T", "NBC", "GCBC"]

for i in range(len(sender_cell)):

    plt = cell2location.plt.plot_spatial(sp_adata,show_img=False,labels=[sender_cell[i],to_cell[i]],
                                        color=[sender_cell[i], to_cell[i]],circle_diameter=3.0,reorder_cmap=[1,3],max_color_quantile=1)
    plt.savefig(f"figures/CCC/{sender_cell[i]}_{to_cell[i]}_distribution.pdf")
    plt.show()
