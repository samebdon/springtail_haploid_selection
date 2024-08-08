remove repeats and missing sites

bedtools complement -i ../SBGE_dataset_scripting/iphiclides/annotation/iphiclides_podalirius.IP_504.v2_0.filteredRepeats.bed -g diem_files/concat/gimble_fix/IP_genomefile.tsv > iphiclides_norepeats.bed

bcftools filter -O z --include "N_MISSING=0" iphiclides.IP_IF.20220207.vcf.gz > iphiclides.IP_IF.20220207.no_missing.vcf.gz && bcftools index -c iphiclides.IP_IF.20220207.no_missing.vcf.gz -o iphiclides.IP_IF.20220207.no_missing.vcf.gz.csi

bcftools view -R iphiclides_norepeats.bed -o iphiclides.IP_IF.20220207.no_missing.no_repeats.vcf.gz iphiclides.IP_IF.20220207.no_missing.vcf.gz && bcftools index -c iphiclides.IP_IF.20220207.no_missing.no_repeats.vcf.gz -o iphiclides.IP_IF.20220207.no_missing.no_repeats.vcf.gz.csi

plink --vcf iphiclides.IP_IF.20220207.no_missing.no_repeats.vcf.gz --homozyg-kb 100 --homozyg-window-het 10 --homozyg-window-snp 1000 --homozyg-window-threshold 0.001 --allow-extra-chr






bedtools complement -i iphiclides_podalirius.IP_504.v2_0.filteredRepeats.bed -g IP_genomefile.tsv > iphiclides_norepeats.bed

bcftools filter -O z --include "N_MISSING=0" iphiclides.IP_IF.20220207.vcf.gz > iphiclides.IP_IF.20220207.no_missing.vcf.gz

bcftools index -c iphiclides.IP_IF.20220207.no_missing.vcf.gz -o iphiclides.IP_IF.20220207.no_missing.vcf.gz.csi

bcftools view -R iphiclides_norepeats.bed -o iphiclides.IP_IF.20220207.no_missing.no_repeats.vcf.gz iphiclides.IP_IF.20220207.no_missing.vcf.gz

bcftools index -c iphiclides.IP_IF.20220207.no_missing.no_repeats.vcf.gz -o iphiclides.IP_IF.20220207.no_missing.no_repeats.vcf.gz.csi

plink --vcf iphiclides.IP_IF.20220207.no_missing.no_repeats.vcf.gz --homozyg-kb 100 --homozyg-window-het 10 --homozyg-window-snp 1000 --homozyg-window-threshold 0.001 --allow-extra-chr