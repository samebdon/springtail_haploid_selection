#!/bin/bash

#BSUB -o logs/braker3.out.%J
#BSUB -e logs/braker3.err.%J
#BSUB -q long
#BSUB -n 32
#BSUB -M40960
#BSUB -R "select[mem>40960] rusage[mem=40960]"

module load ISG/singularity/3.11.4
export BRAKER_SIF=/software/team360/src/braker3.sif
singularity exec -B ${PWD}:${PWD} ${BRAKER_SIF} braker.pl \
--genome=./data/lycoriella_braker3_test/core.fa \
--bam=./data/lycoriella_braker3_test/All_rnaseq.bam \
--softmasking \
--workingdir=./data/lycoriella_braker3_test/braker3/lycoriella \
--threads 32 \
--species=lycoriella_core_b3 \
--gff3 \
--prot_seq=./data/lycoriella_braker3_test/Arthropoda+contarina.simple_header.fa \
--useexisting

# genemark key in home directory?
# ~/.gm_key