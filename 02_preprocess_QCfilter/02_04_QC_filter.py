""" 
author: Amparo Roig Adam
created: 2021/03/23
Description: Filter low quality cells according to QC metrics for Stxbp1 HZ and WT mice S1

Inputs:
    -infile     input file to run the QC on. Should be an Anndata file.
    -outfile    outout file name. Def = NULL won't save the results.
    -outdir     output directory where plots will be saved.

Returns:
    File with calculated QC metrics and pass/fail flag for QC filters for each cell
    File with calculated QC metrics and only cells passing QC filter
"""

from datetime import date
import anndata
import scanpy as sc

#THRESHOLDS
th_min_n_genes = 1000
#th_max_n_genes = 12000
th_max_mito = 10 #maximum fraction of mitochondrial genes
th_min_umi = 500

#read input file
adataStxbp1 = anndata.read('../../../../results/20210928_stxbp1_prefilter_relabelled_QC.h5ad')
print(adataStxbp1)

#Flag cells to filter
adataStxbp1.obs['QC_pass'] = 'NA'

adataStxbp1.obs.loc[(adataStxbp1.obs.n_genes > th_min_n_genes) & (adataStxbp1.obs.total_counts > th_min_umi) & (adataStxbp1.obs.percent_mito < th_max_mito), 'QC_pass'] = 'pass'
adataStxbp1.obs.loc[(adataStxbp1.obs.n_genes <= th_min_n_genes)| (adataStxbp1.obs.total_counts <= th_min_umi) | (adataStxbp1.obs.percent_mito >= th_max_mito), 'QC_pass'] = 'fail'

sc.pl.violin(adataStxbp1, ['n_genes', 'total_counts', 'percent_mito'], groupby = 'QC_pass',jitter=0.4, size = 0., multi_panel=True,save = False)

adataStxbp1.obs['QC_pass'].value_counts()

#save QC flagged object
adataStxbp1.write(filename= 'results/'+date.today().strftime("%Y%m%d")+'_QCflag_stxbp1.h5ad')

#save filtered object
adataStxbp1_filtered = adataStxbp1[adataStxbp1.obs['QC_pass'] == 'pass',:]

    #filter low expressed genes
sc.pp.filter_genes(adataStxbp1_filtered, min_cells=10)

print('n_genes median')
adataStxbp1_filtered.obs['n_genes'].median()

adataStxbp1_filtered.write('results/' + date.today().strftime("%Y%m%d") + '_filtered_stxbp1_sc.h5ad', compression='gzip')
