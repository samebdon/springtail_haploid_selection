#!/bin/bash

#BSUB -o logs/mosdepth.out.%J
#BSUB -e logs/mosdepth.err.%J
#BSUB -q basement
#BSUB -n 4
#BSUB -M 4096
#BSUB -R "select[mem>4096] rusage[mem=4096]"

#mamba activate coverage
cat data/results/bam_paths.txt | cut -f6- -d'/' | cut -f-1 -d'.' | parallel -j1 'mosdepth -n -t 4 --fast-mode --by 500000 data/results/mosdepth/{} {}.GCA_947179485.1.simple_header.earlGrey_masked.bam' 