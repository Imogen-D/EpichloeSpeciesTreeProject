#!/bin/bash

name=$(basename $1 .fna)
echo "Aligning "${name}
mafft --auto --adjustdirection $1 > aligned/${name}_aligned.fna


