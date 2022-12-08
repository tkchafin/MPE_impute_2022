#!/bin/bash
#SBATCH --job-name=subs
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --constraint='douglas&256gb'
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --time=100-00:00:00
#SBATCH -e subset_%j.err
#SBATHC -o subset_%j.txt

export PATH=/home/tkchafin/miniconda3/bin:/home/tkchafin/miniconda3/condabin:$PATH

source activate pgsui

#cp $SLURM_SUBMIT_DIR/$FASTA /scratch/$SLURM_JOB_ID/$FASTA

cd $SLURM_SUBMIT_DIR

python3 ~/local/src/scripts/alignment_subsetter.py -p all_seqs.align_targets.phylip.uniqueseq.phy -o alignment -s 10000 | tee log_subset2.txt
