#!/bin/bash

#BSUB -o logs/annotate.out.%J
#BSUB -e logs/annotate.err.%J
#BSUB -q normal
#BSUB -n 4
#BSUB -M 8192
#BSUB -R "select[mem>8192] rusage[mem=8192]"


#module load braker3/3.0.8--hdfd78af_0
#module load earlgrey/3.0-c1

module load nextflow/23.04.1-5866
module load bedtools/2.31.1--hf5e1c6e_1 
module load ISG/singularity/3.11.4
export BRAKER_SIF=/software/team360/src/braker3.sif
nextflow run nextflow/annotate/annotate_main.nf -c nextflow/conf/nextflow.config -w data/workdir/annotate -resume # -with-report -with-trace -with-timeline -with-dag
