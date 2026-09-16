#!/bin/bash

#BSUB -o logs/braker2.smivir.out.%J
#BSUB -e logs/braker2.smivir.err.%J
#BSUB -q long
#BSUB -n 32
#BSUB -M40960
#BSUB -R "select[mem>40960] rusage[mem=40960]"

module load braker3/3.0.8-c1
braker.pl --genome=./data/genomes/GCA_965194885.1_qlSmiViri2.1_genomic.simple_header.earlgrey_masked.fna --softmasking --workingdir=./data/results/braker2/sminthurus_viridis --threads 32 --species=sminthurus_viridis --gff3 --prot_seq=./data/Arthropoda.allacmaFusca.fa --useexisting
