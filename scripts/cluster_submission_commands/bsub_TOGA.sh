#!/bin/bash

#BSUB -o logs/TOGA.out.%J
#BSUB -e logs/TOGA.err.%J
#BSUB -q basement
#BSUB -n 4
#BSUB -M 4096
#BSUB -R "select[mem>4096] rusage[mem=4096]"

module load nextflow/23.04.1-5866

toga.py data/workdir/make_chains/afusca_svir/GCA_947179485.Svir.final.chain.gz data/results/braker3/allacma_fusca/braker.agat.simpler.TOGA.bed data/results/genomes/allacma_fusca/GCA_947179485.1.simpler_header.earlGrey_masked.2bit data/results/genomes/sminthurus_viridis/Svir.primary_haploid_assembly.earlGrey_masked.2bit --kt --pd data/results/toga/afusca_svir_clem --nc /software/team360/src/TOGA/nextflow_config_files/ --cb 10,100 --cjn 500 --nextflow_dir data/workdir/toga -i data/results/braker3/allacma_fusca/afusca.braker3.isoforms.tsv
toga.py data/workdir/make_chains/afusca_smiaqu/GCA_947179485.SmiAqu.final.chain.gz data/results/braker3/allacma_fusca/braker.agat.simpler.TOGA.bed data/results/genomes/allacma_fusca/GCA_947179485.1.simpler_header.earlGrey_masked.2bit data/results/genomes/sminthurides_aquaticus/SmiAqu.scaffolds.final.earlGrey_masked.2bit --kt --pd data/results/toga/afusca_smiaqu --nc /software/team360/src/TOGA/nextflow_config_files/ --cb 10,100 --cjn 500 --nextflow_dir data/workdir/toga -i data/results/braker3/allacma_fusca/afusca.braker3.isoforms.tsv
