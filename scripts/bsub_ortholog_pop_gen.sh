#!/bin/bash

#BSUB -o logs/ortholog_pop_gen.out.%J
#BSUB -e logs/ortholog_pop_gen.err.%J
#BSUB -q long
#BSUB -n 4
#BSUB -M 8192
#BSUB -R "select[mem>8192] rusage[mem=8192]"

module load nextflow/23.04.1-5866
module load bcftools/1.20--h8b25389_0 
nextflow run $NF_PATH/ortholog_pop_gen/ortholog_pop_gen_main.nf -params-file $NF_PATH/ortholog_pop_gen/afusca_smivir_params.json -c $NF_PATH/conf/nextflow.config -with-conda /software/treeoflife/conda/users/envs/team360/se13/ortholog_pop_gen -w data/workdir/ortholog_pop_gen -resume # -with-report -with-trace -with-timeline -with-dag
