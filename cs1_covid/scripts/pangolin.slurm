#!/bin/bash
#SBATCH --job-name=pangolin_usher
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --constraint='douglas&256gb'
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --time=100-00:00:00
#SBATCH -e pango_%j.err
#SBATCH -o pango_%j.txt

export PATH=/home/tkchafin/miniconda3/bin:/home/tkchafin/miniconda3/condabin:$PATH

module load boost 
source activate usher

export PATH=/home/tkchafin/miniconda3/envs/pango/bin:$PATH
#source activate pango

FASTA=all_seqs.align_targets.gapfix.fasta
OUTPUT=pango_usher_test.csv
RUNMODE="usher"

cp $SLURM_SUBMIT_DIR/$FASTA /scratch/$SLURM_JOB_ID/$FASTA

cd /scratch/$SLURM_JOB_ID

~/miniconda3/envs/pango/bin/pangolin $FASTA --outfile $OUTPUT -t 16 --tempdir /scratch/$SLURM_JOB_ID --max-ambig 0.5 --analysis-mode $RUNMODE

#rm $FASTA
#mv $OUTPUT $SLURM_SUBMIT_DIR/
mv * $SLURM_SUBMIT_DIR/

