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
--genome=/lustre/scratch126/tol/teams/jaron/users/fede/Braker/Ling/core.fa \
--bam=/lustre/scratch126/tol/teams/jaron/users/fede/RNAseq/Ling/hisat2/All_rnaseq.bam \
--softmasking --workingdir=./data/results/braker3/lycoriella \
--GENEMARK_PATH=${ETP}/gmes \
--threads 32 \
--species=lycoriella_core_b3 \
--gff3 \
--prot_seq=/lustre/scratch126/tol/teams/jaron/users/fede/Braker/Arthropoda+contarina.simple_header.fa \
--useexisting

# genemark key in home directory?
# ~/.gm_key