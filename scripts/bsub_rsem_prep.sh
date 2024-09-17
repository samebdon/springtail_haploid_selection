#!/bin/bash

#BSUB -o logs/rsem.out.%J
#BSUB -e logs/rsem.err.%J
#BSUB -q normal
#BSUB -n 4
#BSUB -M 4096
#BSUB -R "select[mem>4096] rusage[mem=4096]"

conda activate dosage_comp
rsem-prepare-reference -gtf data/results/braker3/allacma_fusca/braker.gtf --bowtie data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.fasta data/results/diff_expr/allacma_fusca/rsem/afusca
