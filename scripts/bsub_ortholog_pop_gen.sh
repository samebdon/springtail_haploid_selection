#!/bin/bash

#BSUB -o logs/ortholog_pop_gen.out.%J
#BSUB -e logs/ortholog_pop_gen.err.%J
#BSUB -q normal
#BSUB -n 4
#BSUB -M 8192
#BSUB -R "select[mem>8192] rusage[mem=8192]"

module load nextflow/23.04.1-5866
nextflow run nextflow/ortholog_pop_gen/ortholog_pop_gen_main.nf -c nextflow/conf/nextflow.config -with-conda /software/team360/minisam/envs/ortholog_pop_gen -w data/workdir/ortholog_pop_gen -resume -with-report -with-trace -with-timeline -with-dag
