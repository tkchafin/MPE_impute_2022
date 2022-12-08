#!/bin/bash
#SBATCH --job-name=iqtreeS10
#SBATCH --partition condo
#SBATCH --qos condo
#SBATCH --constraint='abi&3072gb'
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --exclusive
#SBATCH --time=100-00:00:00
#SBATCH -e iqtree_%j.err
#SBATHC -o iqtree_%j.txt

export PATH=/home/tkchafin/miniconda3/bin:/home/tkchafin/miniconda3/condabin:$PATH

cp $SLURM_SUBMIT_DIR/alignment_*.phylip /scratch/$SLURM_JOB_ID/$FASTA

module load parallel

cd /scratch/$SLURM_JOB_ID

doit(){
  iqtree -s $1 -safe -fast -m GTR+I+R3 -nt 4 -wsr
}
export -f doit
 
ls *.phylip | parallel doit {}

#iqtree2 -fast -mem 2800G -nt 48 -s $FASTA -te $TREE -m "GTR{0.14012,0.66091,0.12269,0.09998,2.61404}+F{0.299,0.184,0.196,0.321}+I{0.616}+R5{0.210,0.334,0.219,0.342,0.284,0.344,0.202,1.456,0.085,5.456}" --no-opt-gamma-inv -blmin 0.0000000001 -wsr | tee $SLURM_SUBMIT_DIR/log2.txt

mv ./* $SLURM_SUBMIT_DIR/
