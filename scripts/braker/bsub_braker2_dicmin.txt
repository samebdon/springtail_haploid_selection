#!/bin/bash

#BSUB -o logs/braker2.dicmin.out.%J
#BSUB -e logs/braker2.dicmin.err.%J
#BSUB -q long
#BSUB -n 32
#BSUB -M40960
#BSUB -R "select[mem>40960] rusage[mem=40960]"

module load ISG/singularity/3.9.0
export BRAKER_SIF=/software/team360/src/braker3.sif
singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl --genome=data/results/genomes/dicyrtomina_minuta/qeDicMinu4_1.earlGrey_masked.fa --softmasking --workingdir=./data/results/braker2/dicyrtomina_minuta --GENEMARK_PATH=${ETP}/gmes --threads 32 --species=dicyrtomina_minuta --gff3 --prot_seq=./data/results/braker2/dbs/Arthropoda.allacmaFusca.fa --useexisting
