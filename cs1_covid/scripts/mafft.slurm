#!/bin/bash
#SBATCH --job-name=bigboi_mafft
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --constraint='abi&3072gb'
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --exclusive
#SBATCH --time=100-00:00:00
#SBATCH -e mafft_%j.err
#SBATHC -o mafft_%j.txt

export PATH=/home/tkchafin/miniconda3/bin:/home/tkchafin/miniconda3/condabin:$PATH

TREE="pruned_tree.nwk"
FASTA="all_seqs.fasta"
REF="cov_ref_Wuhan-Hu-1_NC_045512.fasta"

cp $SLURM_SUBMIT_DIR/$TREE /scratch/$SLURM_JOB_ID/$TREE
cp $SLURM_SUBMIT_DIR/$FASTA /scratch/$SLURM_JOB_ID/$FASTA
cp $SLURM_SUBMIT_DIR/$REF /scratch/$SLURM_JOB_ID/$REF

cd /scratch/$SLURM_JOB_ID

mafft --6merpair --thread -1 --keeplength --addfragments $FASTA $REF > all_seqs.align.fasta

mv ./* $SLURM_SUBMIT_DIR/
