#!/bin/bash

#BSUB -o logs/ortholog_pop_gen.out.%J
#BSUB -e logs/ortholog_pop_gen.err.%J
#BSUB -q basement
#BSUB -n 4
#BSUB -M 4096
#BSUB -R "select[mem>4096] rusage[mem=4096]"

module load nextflow/23.04.1-5866
nextflow run nextflow/ortholog_pop_gen/ortholog_pop_gen.nf -c nextflow/conf/nextflow.config -with-conda /software/team360/minisam/envs/ortholog_pop_gen -w data/workdir/ortholog_pop_gen -resume #-with-report -with-trace -with-timeline -with-dag
