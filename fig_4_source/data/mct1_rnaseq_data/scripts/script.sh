#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH -o /uufs/chpc.utah.edu/common/home/rutter-group1/jordan/sce_mct1_rnaseq/slurmjob-%j
#SBATCH --partition=notchpeak

source /uufs/chpc.utah.edu/common/home/u0690617/miniconda3/etc/profile.d/conda.sh
source activate xpresspipe

SCRDIR=/scratch/general/lustre/$USER/$SLURM_JOBID
mkdir -p $SCRDIR

INPUT=/uufs/chpc.utah.edu/common/home/rutter-group1/jordan/sce_mct1_rnaseq/files
REF=/uufs/chpc.utah.edu/common/home/rutter-group1/jordan/references/yeast_ensembl_v100_se49

mkdir $SCRDIR/input
mkdir $SCRDIR/output
cp $INPUT/*.txt.gz $SCRDIR/input/.

cd $SCRDIR/.

#xpresspipe curateReference -o $REF -f $REF/genome_fastas -g $REF/transcripts.gtf --sjdbOverhang 49 --genome_size 12100000

xpresspipe seRNAseq -i $SCRDIR/input -o $SCRDIR/output -r $REF --gtf $REF/transcripts.gtf -e sce_mct1 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA  --sjdbOverhang 49 --genome_size 12100000
