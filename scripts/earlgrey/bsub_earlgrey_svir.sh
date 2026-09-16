#!/bin/bash

#BSUB -o logs/earlGrey.out.%J
#BSUB -e logs/earlGrey.err.%J
#BSUB -q hugemem
#BSUB -n 64
#BSUB -M 400000
#BSUB -R "select[mem>400000] rusage[mem=400000]"

conda activate earlgrey 
earlGrey -g data/genomes/GCA_965194885.1_qlSmiViri2.1_genomic.simple_header.fna -s sminthurusViridis -o data/results/earlGrey/sminthurus_viridis/results -t 64 -d yes
