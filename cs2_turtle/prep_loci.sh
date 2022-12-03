#!/bin/bash

# decompress
tar -xvzf BOX.loci.tar.gz 

# filter -- requires 25% presence and 3 PIS per alignment 
python3 ./filterLoci.py -i BOX.loci -s 105 -p 3 -f -o box_loc | tee filter_log.txt

# move into fasta directory
for f in *.fasta; do 
 new=`echo $f | sed "s/.fasta/.fas/g"`
 mv $f $new 
 gzip $new
done
mkdir fasta
mv *.fas.gz fasta/

# re-compress
tar -cvzf BOX.loci.tar.gz BOX.loci
