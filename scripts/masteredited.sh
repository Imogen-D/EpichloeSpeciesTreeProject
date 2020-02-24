#!/bin/bash

mkdir ./aligned
mkdir ./unaligned

file_name=Epiall  #no extensions required

python3 ./orthologs.py ortho_long.tsv seq_map.tsv 20

#arguments for above; ortho_subset, seq_files, minimum species

for file in ./unaligned/*
do
./auto_align.sh $file
done

mkdir ${file_name}_trees
for file in ./aligned/*
do
filebase=$(basename "$file" _aligned.fna)
echo $filebase
raxmlHPC -s $file -m GTRGAMMA -n ${filebase}_${file_name} -p 10
mv *.${filebase}_${file_name} ${file_name}_trees
done



