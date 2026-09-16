cat data/results/diff_expr/allacma_fusca/rsem/AF_F_1.genes.results | awk '{print $1'} > data/results/diff_expr/allacma_fusca/rsem/rsem_genes.tsv

parallel -j1 "awk '{print \$6'} data/results/diff_expr/allacma_fusca/rsem/{}.genes.results | sed -e 's/TPM/{}/g' > data/results/diff_expr/allacma_fusca/rsem/TPMs/{}_TPM.tsv" :::: scripts/rnaseq_prefixes.txt

parallel -j1 "awk '{print \$7'} data/results/diff_expr/allacma_fusca/rsem/{}.genes.results | sed -e 's/FPKM/{}/g' > data/results/diff_expr/allacma_fusca/rsem/FPKMs/{}_FPKM.tsv" :::: scripts/rnaseq_prefixes.txt

parallel -j1 "awk '{print \$5'} data/results/diff_expr/allacma_fusca/rsem/{}.genes.results | sed -e 's/expected_count/{}/g' > data/results/diff_expr/allacma_fusca/rsem/counts/{}_est_count.tsv" :::: scripts/rnaseq_prefixes.txt

paste -d '\t' data/results/diff_expr/allacma_fusca/rsem/rsem_genes.tsv data/results/diff_expr/allacma_fusca/rsem/TPMs/* > data/results/diff_expr/allacma_fusca/rsem/rsem_TPM_agg.tsv

paste -d '\t' data/results/diff_expr/allacma_fusca/rsem/rsem_genes.tsv data/results/diff_expr/allacma_fusca/rsem/FPKMs/* > data/results/diff_expr/allacma_fusca/rsem/rsem_FPKM_agg.tsv
