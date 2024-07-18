#!/bin/bash

#BSUB -o logs/earlGrey.out.%J
#BSUB -e logs/earlGrey.err.%J
#BSUB -q basement
#BSUB -n 64
#BSUB -M 204800
#BSUB -R "select[mem>204800] rusage[mem=204800]"

earlGrey -g /lustre/scratch126/tol/teams/jaron/data/assemblies_Sanger/insects/Bradysia_impatiens/assembly/curated/idBraImpa2.1/idBraImpa2.1.primary.curated.fa -s bradysiaImpatiens -o data/results/earlGrey/bradysia_impatiens/results -t 64