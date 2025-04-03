""" 
author: Amparo Roig Adam

Description: Doublet removal uses the raw UMI cout matrix and calculates a doublet score for each cell

Inputs:
    anndata file with filtered cells according to QC metrics

Returns:
    Anntada file with doublet score and T/F label on doublet annotation in metadata
    Anndata file with doublets filtered out
    Scrublet plots

"""

import scanpy as sc
import anndata
import matplotlib.pyplot as plt
from datetime import date
import scrublet as scr

#OPTIONS
#Save plots
sc.set_figure_params(format = 'pdf', transparent = True, vector_friendly = True, dpi_save = 200)
sc.settings.figdir = './results/'+date.today().strftime("%Y%m%d")+'_scrublet'

saveplots = True

def save_plots(pltname = None, saveplots = saveplots):
    if saveplots:
        if pltname is None:
            pltname = ''
        saveplots = 'results/' + date.today().strftime("%Y%m%d")+ '_scrublet/' + pltname + '.pdf'
        return saveplots
    else:
        return False

#Plotting parameters
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=12)
plt.rcParams['pdf.fonttype'] = 42

#############################
#Load prefiltered data
#adataStxbp1 = anndata.read('results/'+date.today().strftime("%Y%m%d")+'_filtered_stxbp1_sc.h5ad')
adataStxbp1 = anndata.read('../../../../results/20210504_filtered_stxbp1_sc.h5ad')
adataStxbp1

#list of sample IDs
samples = adataStxbp1.obs.SampleID.values.unique()

print(*samples)

#manually set thresholds for doublet score - adjusted after first run with automatic thresholds
manual_thr = [0.39, 0.21, 0.15, 0.22, 0.20, 0.21, 0.20, 0.22, 0.21, 0.19, 0.21, 0.25, 0.22, 0.20]


# Keep doublet_score in final object
adataStxbp1.obs['scrublet_score'] = 0
adataStxbp1.obs['predicted_doublet'] = 'NA'
adataStxbp1.obs['doublet_threshold'] = 0

## Run SCRUBLET for each sample individually
i = -1
for sample_name in samples:
    i += 1
    print(sample_name)
    adataStxbp1_sample = adataStxbp1[adataStxbp1.obs.SampleID == sample_name]

    scrub = scr.Scrublet(adataStxbp1_sample.X, expected_doublet_rate=0.039)

    adataStxbp1_sample.obs['doublet_scores'], adataStxbp1_sample.obs['predicted_doublets'] = scrub.scrub_doublets()

    print('AUTOMATIC THRESHOLD')
    #scrub.plot_histogram()

    if saveplots:
        plt.savefig(save_plots(sample_name +'_scrublet_autothr_score_hist'))

    #sc.pl.umap(adataStxbp1_sample, color = 'doublet_scores', color_map = 'Wistia')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding('UMAP', order_points=True)

    if saveplots:
        plt.savefig(save_plots(sample_name + 'scrublet_autothr_score_umap'))

    plt.show()
    print('MANUALLY ADJUSTED THRESHOLD')
    print('Manually set threshold at doublet score = ' + str (manual_thr[i]))

    scrub.call_doublets(threshold= manual_thr[i])

    scrub.plot_histogram()
    #plt.show()
    if saveplots:
        plt.savefig(save_plots(sample_name + 'scrublet_manthr_score_hist'))
    #sc.pl.umap(adataStxbp1_sample, color = 'doublet_scores', color_map = 'Wistia')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    scrub.plot_embedding('UMAP', order_points=True)

    if saveplots:
        plt.savefig(save_plots(sample_name +'scrublet_manthr_score_umap'))

    plt.show()

    print('Total cells in sample: ' + str(adataStxbp1_sample.shape[0]))
    print('Cells kept after filtering doublets: ' + str(adataStxbp1_sample[adataStxbp1_sample.obs.doublet_scores < manual_thr[i]].shape[0]))

    #Save scrublet results into anndata object
    adataStxbp1.obs.loc[adataStxbp1_sample.obs_names,'scrublet_score'] = adataStxbp1_sample.obs['doublet_scores']
    adataStxbp1.obs.loc[adataStxbp1_sample.obs_names,'predicted_doublet'] = adataStxbp1_sample.obs.doublet_scores > manual_thr[i]
    adataStxbp1.obs.loc[adataStxbp1_sample.obs_names,'doublet_threshold'] = manual_thr[i]

#Check object was updated correctly
adataStxbp1.obs['scrublet_score'].value_counts()
adataStxbp1.obs['predicted_doublet'].value_counts()
adataStxbp1.obs['doublet_threshold'].value_counts()

#Filter out doublets
adataStxbp1_doublet_filter = adataStxbp1[adataStxbp1.obs.predicted_doublet == False]

#Filter out bad quality samples - determined by QC metrics and outlier in PCA 
adataStxbp1_doublet_filter = adataStxbp1_doublet_filter[adataStxbp1_doublet_filter.obs.SampleID != samples[0]]
adataStxbp1_doublet_filter = adataStxbp1_doublet_filter[adataStxbp1_doublet_filter.obs.SampleID != samples[1]]

#Final dataset
print(adataStxbp1_doublet_filter)
adataStxbp1_doublet_filter.obs.SampleID.values.unique()

#Plot cells
figure, axes = plt.subplots()
xlim, ylim = 150000, 10500

axes.set_xlim(-5500,xlim)
axes.set_ylim(-500, ylim)

sc.pl.scatter(adataStxbp1_doublet_filter, x='total_counts', y='n_genes', ax = axes,color = 'percent_mito',save= save_plots(pltname = 'prefilter' + 'ngenesUMI'))

#Filtered dataset - write to file
adataStxbp1_doublet_filter.write('results/' + date.today().strftime("%Y%m%d") + '_final_filter.h5ad', compression='gzip')
#Unfiltered dataset - metadata indicating doublet pass or not
adataStxbp1.write('results/' + date.today().strftime("%Y%m%d") + '_scrublet_unfiltered.h5ad', compression='gzip')


