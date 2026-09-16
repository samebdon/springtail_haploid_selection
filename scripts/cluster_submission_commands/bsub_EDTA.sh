#!/bin/bash

#BSUB -o logs/earlGrey.out.%J
#BSUB -e logs/earlGrey.err.%J
#BSUB -q week
#BSUB -n 64
#BSUB -M 102400
#BSUB -R "select[mem>102400] rusage[mem=102400]"

conda activate EDTA
module load earlgrey/3.0-c1
EDTA.pl --genome data/genomes/GCA_947179485.1_qeAllFusc8.1_genomic.simple_header.fna --anno 1 --threads 64
