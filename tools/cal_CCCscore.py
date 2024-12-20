import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
from scipy.stats import pearsonr,spearmanr
from sklearn.metrics import mean_squared_error

def CCCscore(x, y):
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    rho = np.corrcoef(x, y)[0][1]
    sd_x = np.std(x)
    sd_y = np.std(y)
    numerator = 2 * rho * sd_x * sd_y
    denominator = sd_x ** 2 + sd_y ** 2 + (x_mean - y_mean) ** 2
    ccc = numerator / denominator
    return ccc

# -----------------------------
# Pseudo code
results = []
celltype_orders = ['xxx','xxx']
pred_df = pd.read_csv(temp_pred_path, index_col=0)
gd_df = pd.read_csv(temp_gd_path, index_col=0)

for j, ct in enumerate(celltype_orders):
    p = pred_df[ct].values
    g = gd_df[ct].values

    try:                            
        pcc = pearsonr(p, g)[0]                            
    except ValueError:
        pcc = 0
        
    try:
        rmse = np.sqrt(mean_squared_error(p,g))
    except ValueError:                            
        rmse = 1  

    try:
        ccc = CCCscore(p, g)
    except ValueError:                            
        ccc = 0

    try:
        spcc = spearmanr(p, g)[0]
    except ValueError:                            
        spcc = 0

    result = {'Seed':seed,'Method':method,'Cell Type': ct,'CCC': ccc,'PCC': pcc, 'SPCC': spcc, 'RMSE': rmse}
    results.append(result)

results_df = pd.DataFrame(results)
    
    
