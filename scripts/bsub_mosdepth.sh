#!/bin/bash

#BSUB -o logs/mosdepth.out.%J
#BSUB -e logs/mosdepth.err.%J
#BSUB -q normal
#BSUB -n 4
#BSUB -M 4096
#BSUB -R "select[mem>4096] rusage[mem=4096]"

#mamba activate coverage
# need sorted and indexed bams
# with sorting and indexing this is maybe already enough to be a little nextflow workflow
#ls *earlGrey_masked.bam | cut -f-1 -d'.' | parallel -j1 'samtools sort -@ 8 -o {}.GCA_947179485.1.simple_header.earlGrey_masked.sorted.bam {}.GCA_947179485.1.simple_header.earlGrey_masked.bam'
cat data/results/bam_paths.txt | cut -f6- -d'/' | cut -f-1 -d'.' | parallel -j1 'mosdepth -n -t 4 --fast-mode --by 500000 data/results/mosdepth/{} data/results/var_call/allacma_fusca/bwamem/{}.GCA_947179485.1.simple_header.earlGrey_masked.sorted.bam' 