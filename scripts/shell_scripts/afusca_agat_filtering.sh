agat_convert_sp_gxf2gxf.pl --gff data/results/braker3/allacma_fusca/braker.gtf -o data/results/braker3/allacma_fusca/braker.agat.gff3
agat_sp_filter_incomplete_gene_coding_models.pl --gff data/results/braker3/allacma_fusca/braker.agat.gff3 --fasta data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.earlGrey_masked.fasta -o data/results/braker3/allacma_fusca/braker.agat.fix.models.gff3
agat_sp_fix_cds_phases.pl --gff data/results/braker3/allacma_fusca/braker.agat.fix.models.gff3 -f data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.earlGrey_masked.fasta -o data/results/braker3/allacma_fusca/braker.agat.fix.models.phases.gff3
agat_sp_fix_overlaping_genes.pl -f data/results/braker3/allacma_fusca/braker.agat.fix.models.phases.gff3 -o data/results/braker3/allacma_fusca/braker.agat.fix.models.phases.overlap.gff3
agat_sp_flag_premature_stop_codons.pl --gff data/results/braker3/allacma_fusca/braker.agat.fix.models.phases.overlap.gff3 --fasta data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.earlGrey_masked.fasta  --out data/results/braker3/allacma_fusca/braker.agat.fix.models.phases.overlap.pseudo.gff3
agat_sp_statistics.pl --gff data/results/braker3/allacma_fusca/braker.agat.fix.models.phases.overlap.pseudo.gff3

agat_sp_keep_longest_isoform.pl -gff data/results/braker3/allacma_fusca/braker.agat.fix.models.phases.overlap.pseudo.gff3 -o data/results/braker3/allacma_fusca/braker.agat.fix.models.phases.overlap.pseudo.longest_isoform.gff3
agat_sp_statistics.pl --gff data/results/braker3/allacma_fusca/braker.agat.fix.models.phases.overlap.pseudo.longest_isoform.gff3
