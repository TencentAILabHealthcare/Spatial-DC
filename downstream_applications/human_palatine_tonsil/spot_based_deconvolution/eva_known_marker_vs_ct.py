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


# 比较已知marker与预测出来的细胞比例的PCC value
# Panel-4, Benchmark metrics with known markers
# Panel-4, Benchmark metrics with known markers
os.chdir("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/01_cell_perc")

methods = ["SpatialDC","Tangram","Cell2location","Stereoscope","DestVI","Seurat","spotlight","spatialdwls","CARD"]
sp_adata = sc.read_h5ad("01_data/HumanTonsil_Spatial_2492.h5ad")
sc.pp.normalize_total(sp_adata)
sc.pp.log1p(sp_adata)
# scaler = MinMaxScaler(feature_range=(0, 1))
# sp_adata.X = scaler.fit_transform(sp_adata.to_df())
sc.pp.scale(sp_adata)

B_value = sp_adata.to_df()["CD19"]
T_value = sp_adata.to_df()["CD3"]
Other_value = sp_adata.to_df()["CD45"]
Epi_value = sp_adata.to_df()["CD324"]

# 用于存储所有方法的Pearson相关系数的列表  
pcc_Bs = []  
pcc_Ts = []  
pcc_Leukocytes = []  
pcc_Epis = [] 

spcc_Bs = []  
spcc_Ts = []  
spcc_Leukocytes = []  
spcc_Epis = [] 

for m in methods:
    sp_adata = sc.read_h5ad("01_data/HumanTonsil_Spatial_2492.h5ad")
    ct_order = ['Activated NBC','NBC','GCBC','MBC', 'PC', # Plasma cells
                 'CD4 T', 'CD8 T','Naive CD4 T','Naive CD8 T', 'cycling T', 'DN', 'NK',  #DN, double-negative T cells with a profile of proinflamatory activation
                 'DC','PDC', 'FDC', 'Mono/Macro',  'Granulocytes', # plasmactyoid DC
          'ILC',  'epithelial']
    
    pred_df = pd.read_csv(f"03_benchmark_methods/exp_ref46/{m}/{m}.csv",index_col=0)
    scaler = MinMaxScaler(feature_range=(0, 1))
    
    pred_df["B"] = pred_df[['Activated NBC','NBC','GCBC','MBC', 'PC']].sum(axis=1).values        
    pred_df["T"] = pred_df[['CD4 T', 'CD8 T','Naive CD4 T','Naive CD8 T', 'cycling T', 'DN']].sum(axis=1).values
    pred_df["Leukocytes"] = pred_df[['B','T','NK','DC','PDC', 'FDC', 'Mono/Macro',  'Granulocytes','ILC']].sum(axis=1).values
    pred_df["Epi"] = pred_df[['epithelial']].values
#     pred_df = pd.DataFrame(scaler.fit_transform(pred_df), index=pred_df.index, columns=pred_df.columns)
    
    temp_adata = sp_adata.copy()
    temp_adata.obs = pred_df

    

    pcc1 = pearsonr(pred_df["B"], B_value)[0]
    pcc2 = pearsonr(pred_df["T"], T_value)[0]
    pcc3 = pearsonr(pred_df["Leukocytes"], Other_value)[0]
    pcc4 = pearsonr(pred_df["Epi"], Epi_value)[0]
    
    spcc1 = spearmanr(pred_df["B"], B_value)[0]
    spcc2 = spearmanr(pred_df["T"], T_value)[0]
    spcc3 = spearmanr(pred_df["Leukocytes"], Other_value)[0]
    spcc4 = spearmanr(pred_df["Epi"], Epi_value)[0]
           
    # 计算平均值  
    mean = np.mean((pcc1,pcc2,pcc3,pcc4))
    mean_2 = np.mean((spcc1,spcc2,spcc3,spcc4))
    
    pcc_Bs.append(pcc1)  
    pcc_Ts.append(pcc2)  
    pcc_Leukocytes.append(pcc3)  
    pcc_Epis.append(pcc4)   
    
    spcc_Bs.append(spcc1)  
    spcc_Ts.append(spcc2)  
    spcc_Leukocytes.append(spcc3)  
    spcc_Epis.append(spcc4)       

#     print(f"{m}:B_{spcc1:.3f}; T_{spcc2:.3f}; Leukocytes_{spcc3:.3f}; Epi_{spcc4:.3f}; SPCC_Mean_{mean_2:.3f}\n")

# PCC
# 将结果转换为DataFrame以便绘图  
pcc_df = pd.DataFrame({  
    'B_PCC': pcc_Bs,  
    'T_PCC': pcc_Ts,  
    'Leukocytes_PCC': pcc_Leukocytes,  
    'Epi_PCC': pcc_Epis  
}, index=methods)  

# 绘制箱型图  
plt.figure(figsize=(5, 3))  
ax=pcc_df.T.boxplot()  
plt.title('Pearson Correlation Coefficients by Method')  
plt.ylabel('Pearson Correlation Coefficient')  
  
# 旋转x轴标签  
plt.xticks(rotation=45, ha='right')  # ha='right' 可以帮助避免标签和箱线图重叠  
# 去除网格线  
ax.grid(False)  # 这将去除网格线  
ax.spines['bottom'].set_color('grey')  # 设置下边框颜色为灰色（或其他您喜欢的颜色）  
ax.spines['left'].set_color('grey')  # 设置左边框颜色为灰色  

# 使用tight_layout自动调整子图参数, 使之填充整个图像区域  
plt.tight_layout()  
plt.savefig(f"figures/HumanTonsil_eva_markers_pcc_boxplot.pdf")
plt.show()

# SPCC
# 将结果转换为DataFrame以便绘图  
spcc_df = pd.DataFrame({  
    'B_PCC': spcc_Bs,  
    'T_PCC': spcc_Ts,  
    'Leukocytes_PCC': spcc_Leukocytes,  
    'Epi_PCC': spcc_Epis  
}, index=methods)  

# 绘制箱型图  
plt.figure(figsize=(5, 3))  
ax=spcc_df.T.boxplot()  
plt.title('Spearman Correlation Coefficients by Method')  
plt.ylabel('Spearman Correlation Coefficient')  
  
# 旋转x轴标签  
plt.xticks(rotation=45, ha='right')  # ha='right' 可以帮助避免标签和箱线图重叠  
# 去除网格线  
ax.grid(False)  # 这将去除网格线  
ax.spines['bottom'].set_color('grey')  # 设置下边框颜色为灰色（或其他您喜欢的颜色）  
ax.spines['left'].set_color('grey')  # 设置左边框颜色为灰色  

# 使用tight_layout自动调整子图参数, 使之填充整个图像区域  
plt.tight_layout()  
plt.savefig(f"figures/HumanTonsil_eva_markers_spcc_boxplot.pdf")
plt.show()

spcc_df.to_csv("figures/HumanTonsil_eva_markers_spcc_boxplot.csv")
pcc_df.to_csv("figures/HumanTonsil_eva_markers_pcc_boxplot.csv")
