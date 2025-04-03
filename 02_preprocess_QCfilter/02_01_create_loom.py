# author: Amparo Roig Adam
# created: 2021/02/09
# description: Generate loom file from cellranger aggr output

import loompy

path_data = "/proj/snic2020-6-62/STXBP1/data/"

sample = 'aggregated_data/aggregated_stxbp1'

indir = path_data + sample
outdir = path_data + "loom_files/"

loompy.create_from_cellranger(indir,outdir,"refdata-cellranger-mm10-3.0.0")
