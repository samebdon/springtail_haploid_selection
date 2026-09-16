agat_sp_filter_feature_from_keep_list.pl -gff data/results/braker3/allacma_fusca/agat/braker.agat.fix.models.phases.overlap.pseudo.gff3 --keep_list data/results/diff_expr/allacma_fusca/deseq2/v_3/male_biased.txt --output data/results/braker3/allacma_fusca/male_biased.gff3

agat_sp_filter_feature_from_keep_list.pl -gff data/results/braker3/allacma_fusca/agat/braker.agat.fix.models.phases.overlap.pseudo.gff3 --keep_list data/results/diff_expr/allacma_fusca/deseq2/v_3/female_biased.txt --output data/results/braker3/allacma_fusca/female_biased.gff3

agat_sp_filter_feature_from_keep_list.pl -gff data/results/braker3/allacma_fusca/agat/braker.agat.fix.models.phases.overlap.pseudo.gff3 --keep_list data/results/diff_expr/allacma_fusca/deseq2/v_3/unbiased.txt --output data/results/braker3/allacma_fusca/unbiased.gff3

agat_sp_filter_feature_from_keep_list.pl -gff data/results/braker3/allacma_fusca/male_biased.gff3 --keep_list data/results/ortholog_pop_gen/allacma_fusca.vs.sminthurus_viridis/allacma_fusca.sp1.SCO_genes.txt --output data/results/braker3/allacma_fusca/male_biased.sco.gff3

agat_sp_filter_feature_from_keep_list.pl -gff data/results/braker3/allacma_fusca/female_biased.gff3 --keep_list data/results/ortholog_pop_gen/allacma_fusca.vs.sminthurus_viridis/allacma_fusca.sp1.SCO_genes.txt --output data/results/braker3/allacma_fusca/female_biased.sco.gff3

agat_sp_filter_feature_from_keep_list.pl -gff data/results/braker3/allacma_fusca/unbiased.gff3 --keep_list data/results/ortholog_pop_gen/allacma_fusca.vs.sminthurus_viridis/allacma_fusca.sp1.SCO_genes.txt --output data/results/braker3/allacma_fusca/unbiased.sco.gff3
