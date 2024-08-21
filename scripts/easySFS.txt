module load bcftools/1.20--h8b25389_0
module load bedtools/2.31.1--hf5e1c6e_1
conda activate easySFS

bcftools query -l allacma_fusca.hard_filtered.sorted.vcf.gz | awk '{split($0,a,"."); print $1,a[2]}' > pop_file.txt

## add pops after pop file

/software/team360/easySFS/easySFS.py -i allacma_fusca.hard_filtered.sorted.vcf.gz -p pop_file.txt -a -f --preview

/software/team360/easySFS/easySFS.py -i allacma_fusca.hard_filtered.sorted.vcf.gz -p pop_file.txt -a -f --proj 24

bedtools intersect -a allacma_fusca.hard_filtered.sorted.vcf.gz -b allacma_fusca.longest_isoforms.0D.bed -header > afusca_0d.vcf

bgzip afusca_0d.vcf 
tabix afusca_0d.vcf.gz

bedtools intersect -a allacma_fusca.hard_filtered.sorted.vcf.gz -b allacma_fusca.longest_isoforms.4D.bed -header > afusca_4d.vcf

bgzip afusca_4d.vcf 
tabix afusca_4d.vcf.gz

bcftools view afusca_0d.vcf.gz --regions OX359245.1,OX359246.1,OX359247.1,OX359248.1 > afusca_0d.A.vcf

bcftools view afusca_4d.vcf.gz --regions OX359245.1,OX359246.1,OX359247.1,OX359248.1 > afusca_4d.A.vcf

bcftools view afusca_0d.vcf.gz --regions OX359249.1,OX359250.1 > afusca_0d.X.vcf

bcftools view afusca_4d.vcf.gz --regions OX359249.1,OX359250.1 > afusca_4d.X.vcf

/software/team360/easySFS/easySFS.py -i afusca_0d.A.vcf -p pop_file.txt -a -f --preview
# 24 highest
/software/team360/easySFS/easySFS.py -i afusca_4d.A.vcf -p pop_file.txt -a -f --preview
# 24 highest
/software/team360/easySFS/easySFS.py -i afusca_0d.X.vcf -p pop_file.txt -a -f --preview
# 24 highest
/software/team360/easySFS/easySFS.py -i afusca_4d.X.vcf -p pop_file.txt -a -f --preview
# 24 highest

/software/team360/easySFS/easySFS.py -i afusca_0d.A.vcf -p pop_file.txt -a -f --proj 24 -o degen_splits/0D_A/
/software/team360/easySFS/easySFS.py -i afusca_4d.A.vcf -p pop_file.txt -a -f --proj 24 -o degen_splits/4D_A/
/software/team360/easySFS/easySFS.py -i afusca_0d.X.vcf -p pop_file.txt -a -f --proj 24 -o degen_splits/0D_X/
/software/team360/easySFS/easySFS.py -i afusca_4d.X.vcf -p pop_file.txt -a -f --proj 24 -o degen_splits/4D_X/

## TO DO: 
## create intergenic bed file
## Subset whole file to intergenic
## Split to A and X
## Run on Whole, A, and X intergenic