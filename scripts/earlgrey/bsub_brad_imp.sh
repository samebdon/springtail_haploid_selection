#!/bin/bash

#BSUB -o logs/earlGrey.out.%J
#BSUB -e logs/earlGrey.err.%J
#BSUB -q long
#BSUB -n 64
#BSUB -M 500000
#BSUB -R "select[mem>500000] rusage[mem=500000]"

earlGrey -g data/results/bradysia_earlgrey/idBraImpa2.1.primary.curated.fa -s bradysiaImpatiens -o data/results/bradysia_earlgrey/ -t 64
