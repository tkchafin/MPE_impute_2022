#!/bin/bash 
#SBATCH --job-name=iqtreeB
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

module load parallel

FASTA="fasta"

cp -r $SLURM_SUBMIT_DIR/$FASTA /scratch/$SLURM_JOB_ID

cd /scratch/$SLURM_JOB_ID


CORES=48 #number of parallel IQ-TREE runs to do 

run_iqtree(){

  ALN=$1
  base=`echo $ALN | sed "s/.fas.gz//g"`
  gunzip $ALN
  if [ ! -f $base".nwk" ]; then
    echo "$base missing treefile, running IQ-TREE..."
    iqtree -s $base".fas" -nt 1 -safe -m MFP -redo -pre $base
    gzip $base".fas"
    rm $base".iqtree"
    rm $base".model.gz"
    rm $base".log"
    rm $base".bionj"
    rm $base".mldist"
    rm $base".uniqueseq.phy"
    rm $base".ckp.gz"
    mv $base".treefile" $base".nwk"
  fi
}

export -f run_iqtree

cd $FASTA

ls *.fas.gz | parallel -j $CORES run_iqtree {}

cp * $SLURM_SUBMIT_DIR/$FASTA

