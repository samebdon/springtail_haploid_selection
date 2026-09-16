perl /software/team360/gFACs-master/gFACs.pl \
	-f braker_2.1.2_gtf \
	--statistics \
	--rem-all-incompletes \
	--sort-by-chromosome \
	--fasta data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.earlGrey_masked.fasta \
	--rem-genes-without-start-and-stop-codon \
	--splice-table \
	--nt-content \
	-O data/results/braker3/gFACS/ \
	data/results/braker3/allacma_fusca/braker.gtf