#!/bin/bash

#BSUB -o logs/earlGrey.out.%J
#BSUB -e logs/earlGrey.err.%J
#BSUB -q long
#BSUB -n 64
#BSUB -M 204800
#BSUB -R "select[mem>204800] rusage[mem=204800]"

##module load earlgrey/4.4.0--h4ac6f70_0
conda activate earlgrey
earlGrey -g data/results/genomes/lipothrix_lubbocki/GCA_034696945.1_ASM3469694v1_genomic.fna -s lipothrixLubbocki -o data/results/earlGrey/lipothrix_lubbocki/results -t 64 -d yes
