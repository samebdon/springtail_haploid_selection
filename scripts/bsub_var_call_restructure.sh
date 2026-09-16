#!/bin/bash

#BSUB -o logs/var_call.out.%J
#BSUB -e logs/var_call.err.%J
#BSUB -q normal
#BSUB -n 4
#BSUB -M 4096
#BSUB -R "select[mem>4096] rusage[mem=4096]"

export MOSDEPTH_Q0=NO_COVERAGE
export MOSDEPTH_Q1=LOW_COVERAGE
export MOSDEPTH_Q2=CALLABLE
export MOSDEPTH_Q3=HIGH_COVERAGE

module load nextflow/23.04.1-5866
module load bcftools/1.20--h8b25389_0
module load picard/3.1.1--hdfd78af_0

nextflow run $NF_PATH/var_call -params-file $NF_PATH/var_call/params/afusca_params.json -with-conda /software/treeoflife/conda/users/envs/team360/se13/var_call -w data/workdir/var_call -resume
