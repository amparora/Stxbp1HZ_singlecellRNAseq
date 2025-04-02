#! /bin/bash -l
#SBATCH -A snic2018-8-143
#SBATCH -p node -n 20
#SBATCH -t 40:00:00
#SBATCH -J cr5_aggregate_stxbp1
#SBATCH -o 20200211_aggregate_stxbp1.out

module load bioinfo-tools
module load bcl2fastq/2.20.0

path_cr=/proj/snic2020-6-62/STXBP1/tools/cellranger-5.0.1/

${path_cr}cellranger aggr --id=aggregated_stxbp1 --csv=/proj/snic2020-6-62/STXBP1/data_cr5/aggregated_data/agg_stxbp1_libraries.csv --normalize=none --nosecondary

