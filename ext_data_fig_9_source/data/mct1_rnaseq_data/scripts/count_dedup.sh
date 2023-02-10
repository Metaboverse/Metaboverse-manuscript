#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o /uufs/chpc.utah.edu/common/home/rutter-group1/jordan/sce_mct1_rnaseq/slurmjob-%j
#SBATCH --partition=notchpeak

source /uufs/chpc.utah.edu/common/home/u0690617/miniconda3/etc/profile.d/conda.sh
source activate xpresspipe

INPUT=/scratch/general/lustre/u0690617/1314327/output/alignments_coordinates
REF=/uufs/chpc.utah.edu/common/home/rutter-group1/jordan/references/yeast_ensembl_v100_se49
PROJ=/uufs/chpc.utah.edu/common/home/rutter-group1/jordan/sce_mct1_rnaseq/deduped_counts

xpresspipe count -i $INPUT -o $PROJ -g $REF/transcripts.gtf -e sce_mct1_deduped_counts -c htseq --bam_suffix _dedupRemoved.bam





