#!/bin/bash

#BSUB -o logs/earlGrey.out.%J
#BSUB -e logs/earlGrey.err.%J
#BSUB -q basement
#BSUB -n 64
#BSUB -M 204800
#BSUB -R "select[mem>204800] rusage[mem=204800]"

earlGrey -g data/results/bradysia_earlgrey/idBraImpa2.1.primary.curated.fa -s bradysiaImpatiens -o data/results/bradysia_earlgrey/ -t 64