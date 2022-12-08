#!/bin/bash
#SBATCH --job-name=treecl
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --constraint='douglas&256gb'
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --time=100-00:00:00
#SBATCH -e treecl_%j.err
#SBATCH -o treecl_%j.txt

export PATH=/home/mussmann/local/bin:/home/tkchafin/local/src/treeCl/bin:$PATH

cp -r $SLURM_SUBMIT_DIR/* /scratch/$SLURM_JOB_ID/$NEX

cd /scratch/$SLURM_JOB_ID

python3 ./treeCl_pipeline.py

mv *.tsv $SLURM_SUBMIT_DIR/
