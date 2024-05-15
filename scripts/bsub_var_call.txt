#!/bin/bash

#BSUB -o logs/var_call.out.%J
#BSUB -e logs/var_call.err.%J
#BSUB -q basement
#BSUB -n 4
#BSUB -M 4096
#BSUB -R "select[mem>4096] rusage[mem=4096]"

export MOSDEPTH_Q0=NO_COVERAGE
export MOSDEPTH_Q1=LOW_COVERAGE
export MOSDEPTH_Q2=CALLABLE
export MOSDEPTH_Q3=HIGH_COVERAGE

module load nextflow/23.04.1-5866
nextflow run nextflow/var_call/var_call_main.nf -c nextflow/conf/nextflow.config -with-conda /software/team360/miniconda/envs/lg_het -w data/workdir/var_call -resume #-with-report -with-trace -with-timeline -with-dag
