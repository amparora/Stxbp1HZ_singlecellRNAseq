#! /bin/bash -l
#SBATCH -A snic2018-8-143
#SBATCH -p node -n 16
#SBATCH -t 24:00:00
#SBATCH -J demultiplexing_stxbp1
#SBATCH -o 20210209_cr5_demultiplexing.out

module load bioinfo-tools
module load bcl2fastq/2.20.0

path_data=/proj/snic2020-6-62/STXBP1/data/raw_data/delivery04233/INBOX/P19312/

path_cr=/proj/snic2020-6-62/STXBP1/tools/cellranger-5.0.1/

${path_cr}cellranger mkfastq --run=${path_data}210119_A00187_0416_AHNHLWDSXY/ --csv=samplesID.csv --localmem=128 --localcores=20
