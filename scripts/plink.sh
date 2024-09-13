module load bcftools/1.20--h8b25389_0
module load bedtools/2.31.1--hf5e1c6e_1
conda activate plink

bedtools sort -i data/results/earlGrey/allacma_fusca/results/allacmaFusca_EarlGrey/allacmaFusca_summaryFiles/allacmaFusca.filteredRepeats.bed -g data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.earlGrey_masked.fasta.fai > data/results/roh/afusca.repeats.sorted.bed
bedtools complement -i data/results/roh/afusca.repeats.sorted.bed -g data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.earlGrey_masked.fasta.fai > data/results/roh/afusca_norepeats.bed
bcftools filter -O z --include "N_MISSING=0" data/results/var_call/allacma_fusca_rm_repeats/allacma_fusca.hard_filtered.sorted.vcf.gz > data/results/roh/afusca.nomissing.vcf.gz
bcftools index -c data/results/roh/afusca.nomissing.vcf.gz -o data/results/roh/afusca.nomissing.vcf.gz.csi
bcftools view -R data/results/roh/afusca_norepeats.bed -o data/results/roh/afusca.nomissing.norepeats.vcf.gz data/results/roh/afusca.nomissing.vcf.gz
bcftools index -c data/results/roh/afusca.nomissing.norepeats.vcf.gz -o data/results/roh/afusca.nomissing.norepeats.vcf.gz.csi
plink --vcf data/results/roh/afusca.nomissing.norepeats.vcf.gz --homozyg-kb 100 --homozyg-window-het 10 --homozyg-window-snp 1000 --homozyg-window-threshold 0.001 --allow-extra-chr
bcftools view data/results/roh/afusca.nomissing.norepeats.vcf.gz --regions OX359245.1,OX359246.1,OX359247.1,OX359248.1 > data/results/roh/afusca.nomissing.norepeats.autosomes.vcf.gz
plink --vcf data/results/roh/afusca.nomissing.norepeats.autosomes.vcf.gz --homozyg-kb 100 --homozyg-window-het 10 --homozyg-window-snp 1000 --homozyg-window-threshold 0.001 --allow-extra-chr

plink --vcf data/results/roh/afusca.nomissing.norepeats.vcf.gz --homozyg-window-het 10 --homozyg-window-snp 1000 --homozyg-window-threshold 0.001 --allow-extra-chr
plink --vcf data/results/roh/afusca.nomissing.norepeats.autosomes.vcf.gz --homozyg-window-het 10 --homozyg-window-snp 1000 --homozyg-window-threshold 0.001 --allow-extra-chr


grep -v 'OX359249' afusca_norepeats.bed > tmp.bed
grep -v 'OX359250' tmp.bed > afusca_norepeats.autosomes.bed
rm tmp.bed
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' afusca_norepeats.bed
awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}' afusca_norepeats.autosomes.bed
