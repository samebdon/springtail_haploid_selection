#!/bin/bash

#BSUB -o logs/annotate.out.%J
#BSUB -e logs/annotate.err.%J
#BSUB -q normal
#BSUB -n 4
#BSUB -M 8192
#BSUB -R "select[mem>8192] rusage[mem=8192]"

module load braker3/3.0.8-c1
conda activate earlgrey
nextflow run nextflow/annotate/annotate_main.nf -c nextflow/conf/nextflow.config -w data/workdir/annotate -resume # -with-report -with-trace -with-timeline -with-dag
