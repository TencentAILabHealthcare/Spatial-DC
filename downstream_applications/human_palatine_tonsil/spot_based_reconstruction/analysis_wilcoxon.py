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

plt.rcParams['font.size'] = 10
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
# plt.rcParams["font.family"] = "Arial"

# wilcoxon
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/02_proteomic_profiles")
scp_adata = sc.read_h5ad("01_data/SpatialDC/SpatialDC_reconstruct_norm.h5ad")
sc.pp.filter_cells(scp_adata, min_counts=10)
sc.pp.filter_genes(scp_adata, min_cells=3)
scp_adata = scp_adata[scp_adata.obs["cellperc"] > 0.01]

ct_order = ['NBC','Activated NBC','GCBC','MBC', 'PC', 
              'CD4 T', 'CD8 T','Naive CD4 T','Naive CD8 T', 'cycling T', 'DN', 'NK', 
              'DC','PDC', 'FDC', 'Mono/Macro',  'Granulocytes', 
       'ILC',  'epithelial']

scp_adata.obs["celltype"] = scp_adata.obs["celltype"].cat.set_categories(ct_order,ordered=False)

sc.pp.normalize_total(scp_adata)
scp_adata.raw = scp_adata.copy()
sc.pp.log1p(scp_adata)

sc.tl.rank_genes_groups(scp_adata, 'celltype', method='wilcoxon',use_raw=True)

result = scp_adata.uns['rank_genes_groups']

result_dict = result
df_list = []
for celltype in result_dict["names"].dtype.names:
    df = pd.DataFrame(result_dict["scores"][celltype], columns=["score"])
    df["pvals_adj"] = result_dict["pvals_adj"][celltype]
    df["logfoldchanges"] = result_dict["logfoldchanges"][celltype]
    df["celltype"] = celltype
    df["feature"] = result_dict["names"][celltype]
    df_list.append(df)
merged_df = pd.concat(df_list)

merged_df.to_csv("analysis/dep.csv")
