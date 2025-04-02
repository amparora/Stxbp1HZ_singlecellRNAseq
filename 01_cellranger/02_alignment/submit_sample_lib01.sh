#! /bin/bash -l
#SBATCH -A snic2018-8-143
#SBATCH -p node -n 16
#SBATCH -t 40:00:00

module load bioinfo-tools
module load bcl2fastq/2.20.0

path_ref=/proj/snic2020-6-62/SCZ_mice_sc/tools/

path_cr=/proj/snic2020-6-62/STXBP1/tools/cellranger-5.0.1/

path_fastq=/proj/snic2020-6-62/STXBP1/data_cr5/demultiplexed_data/HNHLWDSXY/

iden=$1
${path_cr}cellranger count --id="Counts_${iden}" --transcriptome=${path_ref}refdata-cellranger-mm10-3.0.0/ --fastqs ${path_fastq}outs/fastq_path/ --sample $1 --expect-cells 5000 --localmem=120 --localcores=20

