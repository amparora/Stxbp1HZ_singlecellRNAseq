""" 
author: Amparo Roig Adam

Description: Prepare data for easy reading in R to Seurat and running conos packcage
Inputs: Anndata file of filtered data

Returns: Loom file with necessary extra columns for loom to seurat transformation

"""

import scanpy as sc
import anndata
from datetime import date

#read data
adataStxbp1 = anndata.read('../../../../results/20210519_final_filter' + '.h5ad')
print(adataStxbp1)

#shuffle data to avoid plotting bias
adataStxbp1 = sc.pp.subsample(adataStxbp1, fraction = 1, copy = True)

#Save for conos/Seurat/R loading
    #important to have this cell name and gene names as fields here:
adataStxbp1.var['Gene']= adataStxbp1.var_names
adataStxbp1.obs['CellID']= adataStxbp1.obs_names

adataStxbp1.write_loom('results/' + date.today().strftime("%Y%m%d") + '_conos_filt_norm.loom', write_obsm_varm = True)
