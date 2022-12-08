#!/bin/bash
#SBATCH --job-name=degap
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --constraint='douglas&256gb'
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
#SBATCH --time=100-00:00:00
#SBATCH -e degap_%j.err
#SBATHC -o degap_%j.txt

export PATH=/home/tkchafin/miniconda3/bin:/home/tkchafin/miniconda3/condabin:$PATH

FASTA="all_seqs.align_targets.phylip"

cp $SLURM_SUBMIT_DIR/$FASTA /scratch/$SLURM_JOB_ID/$FASTA

cd /scratch/$SLURM_JOB_ID

python3 ~/local/src/scripts/terminalGapRemover.py -p $FASTA 

mv ./* $SLURM_SUBMIT_DIR/
