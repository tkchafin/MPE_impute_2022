#!/bin/bash
#SBATCH --job-name=turt_iqtree
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --constraint='douglas&256gb'
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --time=100-00:00:00
#SBATCH -e iqtree_%j.err
#SBATHC -o iqtree_%j.txt

export PATH=/home/mussmann/local/bin:/home/tkchafin/miniconda3/condabin:$PATH

NEX="BOX_concat_N214.nex"

cp $SLURM_SUBMIT_DIR/$NEX /scratch/$SLURM_JOB_ID/$NEX

cd /scratch/$SLURM_JOB_ID

iqtree -s $NEX -fast -m MFP+ -nt 16 | tee $SLURM_SUBMIT_DIR/log.txt

#iqtree2 -fast -mem 2800G -nt 48 -s $FASTA -te $TREE -m "GTR{0.14012,0.66091,0.12269,0.09998,2.61404}+F{0.299,0.184,0.196,0.321}+I{0.616}+R5{0.210,0.334,0.219,0.342,0.284,0.344,0.202,1.456,0.085,5.456}" --no-opt-gamma-inv -blmin 0.0000000001 -wsr | tee $SLURM_SUBMIT_DIR/log2.txt

mv ./* $SLURM_SUBMIT_DIR/
