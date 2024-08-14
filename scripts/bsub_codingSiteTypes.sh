#!/bin/bash

#BSUB -o logs/codingSiteTypes.out.%J
#BSUB -e logs/codingSiteTypes.err.%J
#BSUB -q normal
#BSUB -n 4
#BSUB -M 4096
#BSUB -R "select[mem>4096] rusage[mem=4096]"

/software/team360/smartin/genomics_general/codingSiteTypes.py -a data/results/braker3/allacma_fusca/agat/braker.agat.fix.models.phases.overlap.pseudo.longest_isoform.gff3 -f gff3 -o data/workdir/codingSiteTypessiteTypes.txt -v data/results/var_call/allacma_fusca_rm_repeats/allacma_fusca.hard_filtered.sorted.vcf.gz -r data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.earlGrey_masked.fasta