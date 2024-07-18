#!/bin/bash

#BSUB -o logs/earlGrey.out.%J
#BSUB -e logs/earlGrey.err.%J
#BSUB -q long
#BSUB -n 64
#BSUB -M 102400
#BSUB -R "select[mem>102400] rusage[mem=102400]"

earlGrey -g data/results/bradysia_earlgrey/idBraImpa2.1.primary.curated.fa -s bradysiaImpatiens -o data/results/bradysia_earlgrey/ -t 64