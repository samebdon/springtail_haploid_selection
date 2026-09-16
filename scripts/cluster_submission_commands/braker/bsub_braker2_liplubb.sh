#!/bin/bash

#BSUB -o logs/braker2.liplubb.out.%J
#BSUB -e logs/braker2.liplubb.err.%J
#BSUB -q long
#BSUB -n 32
#BSUB -M40960
#BSUB -R "select[mem>40960] rusage[mem=40960]"

module load braker3/3.0.8-c1
braker.pl --genome=./data/results/earlGrey/lipothrix_lubbocki/results/lipothrixLubbocki_EarlGrey/lipothrixLubbocki_summaryFiles/lipothrixLubbocki.softmasked.fasta --softmasking --workingdir=./data/results/braker2/lipothrix_lubbocki --threads 32 --species=lipothrix_lubbocki --gff3 --prot_seq=./data/results/braker2/dbs/Arthropoda.allacmaFusca.fa --useexisting
