#!/bin/bash

#BSUB -o logs/rsem_calc.out.%J
#BSUB -e logs/rsem_calc.err.%J
#BSUB -q long
#BSUB -n 32
#BSUB -M 4096
#BSUB -R "select[mem>4096] rusage[mem=4096]"

conda activate dosage_comp

#rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_10.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_10.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_F_10
#rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_1.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_1.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_F_1
#rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_2.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_2.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_F_2
#rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_3.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_3.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_F_3
#rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_4.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_4.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_F_4
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_5.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_5.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_F_5
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_6.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_6.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_F_6
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_7.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_7.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_F_7
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_8.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_8.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_F_8
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_9.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_F_9.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_F_9
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_10.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_10.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_M_10
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_1.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_1.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_M_1
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_2.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_2.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_M_2
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_3.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_3.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_M_3
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_4.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_4.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_M_4
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_5.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_5.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_M_5
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_6.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_6.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_M_6
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_7.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_7.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_M_7
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_8.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_8.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_M_8
rsem-calculate-expression -p 32 --paired-end data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_9.1.fastp.fastq.gz data/results/rnaseq_aln/allacma_fusca/fastp/AF_M_9.2.fastp.fastq.gz data/results/diff_expr/allacma_fusca/rsem/afusca data/results/diff_expr/allacma_fusca/rsem/AF_M_9
