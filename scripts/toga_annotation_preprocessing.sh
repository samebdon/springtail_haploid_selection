faToTwoBit GCA_947179485.1.simpler_header.earlGrey_masked.fasta GCA_947179485.1.simpler_header.earlGrey_masked.2bit

agat_convert_sp_gff2bed.pl --gff data/results/braker3/allacma_fusca/braker.agat.gff3 -o data/results/braker3/allacma_fusca/braker.agat.bed

sed 's/\.1//1' data/results/braker3/allacma_fusca/braker.agat.bed > data/results/braker3/allacma_fusca/braker.agat.simpler.bed

awk 'BEGIN {FS=OFS="\t"}{$5=sprintf("%3.0f",$5)}1' data/results/braker3/allacma_fusca/braker.agat.simpler.bed > data/results/braker3/allacma_fusca/braker.agat.simpler.TOGA.bed
