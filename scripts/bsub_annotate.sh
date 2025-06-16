#!/bin/bash

#BSUB -o logs/annotate.out.%J
#BSUB -e logs/annotate.err.%J
#BSUB -q normal
#BSUB -n 4
#BSUB -M 8192
#BSUB -R "select[mem>8192] rusage[mem=8192]"

module load nextflow/23.04.1-5866
module load braker3/3.0.8-c1
#module load earlgrey/4.4.0--h4ac6f70_0 DOESNT WORK
conda activate earlgrey
nextflow run $NF_PATH/annotate/annotate_main.nf -c $NF_PATH/conf/nextflow.config -w data/workdir/annotate -params-file $NF_PATH/annotate/s_viridis_draft_params.json -resume # -with-report -with-trace -with-timeline -with-dag
