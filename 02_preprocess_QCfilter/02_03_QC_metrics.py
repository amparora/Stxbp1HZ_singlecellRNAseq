""" scRNA-seq QC calculations and plotting for Stxb1 project

author: Amparo Roig Adam
created: 2021/03/23
Description: scRNA-seq QC calculations and plotting for Stxbp1 HZ and WT mice S1

Inputs:
    -infile     input file to run the QC on. Should be an Anndata file.
    -outfile    outout file name. Def = NULL won't save the results.
    -outdir     output directory where plots will be saved.

Returns:
    QC plots
    File with calculated QC metrics 
"""

from datetime import date
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import scvelo
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn import preprocessing
import seaborn as sns

################################################################################

##Arguments - Input and outputs
infile = '../../../data/seq_data/20210504_stxbp1_original.h5ad'
outfile = date.today().strftime("%Y%m%d")+'_stxbp1_prefilter'
outdir = './results/'+date.today().strftime("%Y%m%d")+'_QC/'
pltname = ''

###########

if pltname:
    saveplots = True
else:
    saveplots = False
    pltname = ''

showplots = False

#Run/session info
print(date.today())
sc.logging.print_header()

#OPTIONS
sc.set_figure_params(format = 'pdf', transparent = True, vector_friendly = True, dpi_save = 200)
sc.settings.figdir = outdir
scvelo.set_figure_params(vector_friendly=True)

##Save plots
#saveplots = True
def save_plots(pltname = None, saveplots = saveplots):
    if saveplots:
        if pltname is None:
            pltname = '1'
        saveplots = '_' + pltname + '.pdf'
        return saveplots
    else:
        return False

##Plotting parameters
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=11)
plt.rcParams['pdf.fonttype'] = 42

################################################################################
#Read file and compute QC metrics
#adataStxbp1 = anndata.read('../../data/seq_data/20210504_stxbp1_original.h5ad')
adataStxbp1 = anndata.read(infile)
print(adataStxbp1)
#adataStxbp1.obs = adataStxbp1.obs[['SampleID', 'Sex', 'AgeWeek', 'Genotype', 'TissueExtraction', 'LiveCellsPercent']] #this line allows recalculation of QC metrics when input file already had (for compatibility of pre and postfiltering of cells)

adataStxbp1.var['mt'] = adataStxbp1.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adataStxbp1, qc_vars=['mt'], percent_top=[100], log1p=False, inplace=True)
adataStxbp1.obs = adataStxbp1.obs.rename(columns = {'n_genes_by_counts':'n_genes'})
adataStxbp1.obs = adataStxbp1.obs.rename(columns = {'pct_counts_mt':'percent_mito'})

#print(adataStxbp1.obs['n_genes'].median())

#THRESHOLDS
th_min_n_genes = 1000
#th_max_n_genes = 12000
th_max_mito = 10 #maximum fraction of mitochondrial genes
th_min_umi = 500

##############################################################################################################################
#PLOTS

#overall
sc.pl.violin(adataStxbp1, ['n_genes', 'total_counts', 'percent_mito'],size = 0.2, jitter = 0.4, multi_panel=True, save = save_plots(pltname+'_all'), show = showplots)

sc.pl.violin(adataStxbp1, ['n_genes', 'total_counts', 'percent_mito'],groupby = 'TissueExtraction',size = 0.2, jitter = 0.4, multi_panel=True, save = save_plots(pltname+'_expday'), show = showplots)
sc.pl.violin(adataStxbp1, ['n_genes', 'total_counts', 'percent_mito'], groupby = 'Sex', size = 0.2, jitter = 0.4, multi_panel=True, save = save_plots(pltname+'_sex'), show = showplots)
sc.pl.violin(adataStxbp1, ['n_genes', 'total_counts', 'percent_mito'], groupby = 'Genotype',size = 0.2, jitter = 0.4, multi_panel=True, save = save_plots(pltname+'_genotype'), show = showplots)


sc.pl.scatter(adataStxbp1, x='total_counts', y='percent_mito', save = save_plots(pltname+'_mitoumi'), show = showplots)
sc.pl.scatter(adataStxbp1, x='total_counts', y='n_genes', color='percent_mito', legend_loc = 'best', save = save_plots(pltname+'_NgenesUMI'), show = showplots)

sc.pl.highest_expr_genes(adataStxbp1, n_top=20, save = save_plots(pltname), show = showplots)

#QC metrics per sample
for feature in ['n_genes', 'total_counts', 'percent_mito']:

     figure, axes = plt.subplots()
     axes.set_xticklabels(adataStxbp1.obs.SampleID.values.unique(),rotation=45, ha='right', size = 8)

     sc.pl.violin(adataStxbp1,feature, groupby = 'SampleID', ax = axes,jitter=0, multi_panel=True,save = save_plots(pltname +'_'+ feature), show = showplots)

#Histograms
#MITO GENES & N GENES
def multi_hist(adata, ftrs, thrs):
    nftr = -1
    for ftr in ftrs:
        nftr +=1

        figure, axes = plt.subplots(5,3, sharex = True, figsize = (10,12))

        if ftr == 'percent_mito':
            title = 'Fraction of mitochondrial genes'
        elif ftr == 'n_genes':
            title = 'Number of genes'
        figure.suptitle(title, fontsize = 16, y = 0.95)

        i = -1
        row = -1
        for sample in adata.obs.SampleID.values.unique():
             i += 1
             col = i%3
             if col == 0:
                 row += 1

             axes[row,col].set_title(sample)
             axes[row, col].hist(adata.obs[ftr][adata.obs.SampleID == sample],bins = 50)
             axes[row, col].axvline(thrs[nftr], color = 'orange')

        if saveplots:
             plt.savefig(str(sc.settings.figdir )+'/hists_'+pltname+'_'+ftr+'.pdf')

        #plt.show()

multi_hist(adata = adataStxbp1, ftrs = ['percent_mito', 'n_genes'], thrs = [th_max_mito, th_min_n_genes])

#Scatter plots
def all_scatter(adata, x, y, thrx, thry, save ,color = None):

    figure, axes = plt.subplots()
    axes.plot([0,max(adata.obs[x])], [thry,thry])
    axes.plot([thrx ,thrx ],[0,max(adata.obs[y])])
    #axes.set_xlim(-100,max(adataStxbp1.obs.n_genes)+100)
    #axes.set_ylim(-1)
    if saveplots:
        save = str(sc.settings.figdir )+ '/scatter'+ save

    scvelo.pl.scatter(adata, x=x, y=y ,ax = axes, save = False, color = color, legend_loc = 'right margin', show = showplots, color_map = 'viridis')

    figure.savefig(save)

# nUMI vs % mitochondrial genes + threshold lines
all_scatter(adata = adataStxbp1, x = 'n_genes', y = 'percent_mito', thrx = th_min_n_genes, thry = th_max_mito, save = save_plots(pltname+'_mitogenes'))

#nUMI vs nGenes. Colour code % mitochondrial genes + threshold lines
all_scatter(adata = adataStxbp1, x = 'total_counts', y = 'n_genes', thrx = th_min_umi, thry = th_min_n_genes, color = 'percent_mito', save = save_plots(pltname+'_NgenesUMI'))

#same plot for individual samples
def multi_scatter(adata, x, y, thrx, thry, color = None):
    figure, axes = plt.subplots(5,3, figsize = (20,20))
    figure.subplots_adjust(hspace=0.4, wspace=0.2)

    i = -1
    row = -1
    for sample in adata.obs.SampleID.values.unique():
        i += 1
        col = i%3
        if col == 0:
            row += 1

        axes[row, col].set_title(sample)

        axes[row, col].plot([0,max(adata.obs[x])], [thry,thry])
        axes[row, col].axvline(thrx, color = 'orange')
        axes[row, col].axhline(thry, color = 'blue')

        scvelo.pl.scatter(adata[adata.obs.SampleID == sample], x=x, y=y, ax = axes[row, col], legend_loc ='none' ,color = color, title = sample, save= False, show = False)

    if saveplots:
        plt.savefig(str(sc.settings.figdir )+'/scatter_'+pltname+'_samples_mito.pdf')

multi_scatter(adata = adataStxbp1, x = 'total_counts', y = 'n_genes', thrx = th_min_umi, thry = th_min_n_genes, color = 'percent_mito')

##############################################################################################################################

#run PCA QC metrics to confirm/detect outlier samples 

#collect data
qcmetrics = ['LiveCellsPercent', 'n_genes', 'total_counts', 'pct_counts_in_top_100_genes', 'total_counts_mt','percent_mito']
data = adataStxbp1.obs[qcmetrics]

    # First center and scale the data  
scaled_data = preprocessing.scale(data)

pca = PCA(n_components = 2) # create a PCA object
pca.fit(scaled_data) # do the math
pca_data = pca.transform(scaled_data) # get PCA coordinates for scaled_data


     #The following code constructs the Scree plot
per_var = np.round(pca.explained_variance_ratio_* 100, decimals=1)
print('Explained variance ratio', per_var)
labels = ['PC' + str(x) for x in range(1, len(per_var)+1)]
#
#     plt.bar(x=range(1,len(per_var)+1), height=per_var, tick_label=labels)
#     plt.ylabel('Percentage of Explained Variance')
#     plt.xlabel('Principal Component')
#     plt.title('Scree Plot')
#     plt.show()
#
pca_df = pd.DataFrame(pca_data, index = data.index.values ,columns=labels)
pca_df['SampleID'], pca_df['Sex'], pca_df['Genotype'], pca_df['TissueExtraction'] = adataStxbp1.obs.SampleID, adataStxbp1.obs.Sex, adataStxbp1.obs.Genotype, adataStxbp1.obs.TissueExtraction


#plt.scatter(pca_df.PC1, pca_df.PC2, label = pca_df.SampleID)
if saveplots:
    plt = sns.lmplot(data = pca_df, scatter = False, x = 'PC1', y = 'PC2', hue = 'SampleID', col = 'Sex', row = 'Genotype', order = 2)
    plt.savefig(str(sc.settings.figdir )+'/QC_PCAfitline_samples_'+pltname+'.pdf')

    plt = sns.lmplot(data = pca_df, scatter = True, x = 'PC1', y = 'PC2', hue = 'SampleID', col = 'Sex', row = 'Genotype', order = 2)
    plt.savefig(str(sc.settings.figdir )+'/QC_PCAdatapoints_samples_'+pltname+'.pdf')

    plt = sns.lmplot(data = pca_df, scatter = True, x = 'PC1', y = 'PC2', hue = 'Sex', col = 'Genotype', order = 2)
    plt.savefig(str(sc.settings.figdir )+'/QC_PCAdatapoins_sex_'+pltname+'.pdf')

    plt = sns.lmplot(data = pca_df, scatter = False, x = 'PC1', y = 'PC2', hue = 'Sex', col = 'Genotype', order = 2)
    plt.savefig(str(sc.settings.figdir )+'/QC_PCAfitline_sex_'+pltname+'.pdf')


## get the name of the top 10 measurements (QC metrics) that contribute most to pc1 & pc2.
for pc in range(0,2):
    loading_scores = pd.Series(pca.components_[pc], index=qcmetrics)
    sorted_loading_scores = loading_scores.abs().sort_values(ascending=False)
    top_10_genes = sorted_loading_scores[0:10].index.values

    print(f'PC{pc+1} loadings:\n',loading_scores[top_10_genes])


##save final object with QC_PCA for each cell
adataStxbp1.obs['QC_PCA1'] = pca_df.PC1
adataStxbp1.obs['QC_PCA2'] = pca_df.PC2

adataStxbp1.write(filename= outdir+date.today().strftime("%Y%m%d")+outfile+'.h5ad')
