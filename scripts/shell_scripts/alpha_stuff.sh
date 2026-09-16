module load bcftools/1.20--h8b25389_0

parallel  'bcftools view -R {}_intergenic.bed -o {}.intergenic.vcf allacma_fusca.hard_filtered.sorted.vcf.gz' :::: chroms.txt

conda activate easySFS

pops.txt
BH3-2 pop1
WW1-2 pop1
WW1-4 pop1
WW2-1 pop1
WW2-5 pop1
WW2-6 pop1
WW3-1 pop1
WW5-1 pop1
WW5-3 pop1
WW5-4 pop1
WW5-5 pop1
WW5-6 pop1

parallel './easySFS.py -i {}.intergenic.vcf -p pops.txt -a --proj 24 -o {}_sfs' :::: chroms.txt

parallel -j1 "echo {}; cat {}_sfs/dadi/pop1.sfs | cut -f-13 -d' ' | head -n2 | tail -n1" :::: chroms.txt


DFE alpha bootstrapping
> Confidence intervals for parameter estimates were obtained by 
> bootstrapping by gene 200 times, unless otherwise stated.
How did they do the bootstrapping by gene? 
Do i need to sample X sized datasets from the whole genome and compare alpha inference 
to the autosome to test for alpha significance, or do i need to bootstrap genes from x and A separately
to get confidence intervals for A and X (I think this is what they did?)
This is different to fastdfe bootstrapping by resampling SFSs from the DFE 
I can still write a script that bootstraps by gene in a script for the cluster then
writes the SFS. Dont forget monomorphic sites. How to get monomorphic sites with fastdfe approach?
Fast dfe doesnt do the demographic inference

Can I calculate alpha from the DFE estimated with fastDFE
Maybe better with dfe alpha given i have folded sfs

Resample with easySFS or fastDFE?

2 tests:
Is alpha inferred on the X higher than alpha inferred on the autosomes?
Does each bin in the NeS histogram differ between X and A?


Could try get SFS samples with fastDFE or easySFS
If easySFS write a loop based on my current pipeline
If fastDFE, as follows (test full dataset gets current dfe)
Lets start with fastDFE, if it works is better
fastDFE requires ancestral state of sites to be determined based on the AA field
I dont think we know the ancestral state so this will be random
But we should be able to fold anyway and then its ok?

Input: 
longest isoform gff
genome
VCF 						

What is degeneracy annotation and stratification
It relies on VCF info tags to determine degeneracy, will that work for my vcf?

How to split for resampling

Can I even just use fastDFE to get the discretized dfe plot with confidence intervals?
Still need dfe alpha for alpha
It will do bootstrapping for me. lets do this for plot

1) split GFF into X and A GFFs
2) check sfs looks right
3) get the discretized DFE plot
4) bootstrap
5) figure out how to bootstrap SFSs for alpha

#sh
mamba create -n fastDFE 'python>=3.10,<3.13' pip && mamba activate fastDFE && pip install fastdfe && pip install notebook
#bsub -G team360 -Is -n 16 -M 10240 -R "select[mem>10240] rusage[mem=10240]" bash -l
#conda activate fastDFE

grep -v 'OX359249.1' braker.agat.longest_isoform.gff3 | grep -v 'OX359250.1' | grep -v 'CAMXBY' > allacma_fusca.longest_isoform.autosome.gff3
grep -v 'OX359245.1' braker.agat.longest_isoform.gff3 | grep -v 'OX359246.1' | | grep -v 'OX359247.1' | grep -v 'OX359248.1' | grep -v 'CAMXBY' > allacma_fusca.longest_isoform.x.gff3

import fastdfe as fd

vcf='allacma_fusca.hard_filtered.sorted.vcf.gz'
fasta='GCA_947179485.1.simple_header.earlGrey_masked.fasta'
a_gff='allacma_fusca.longest_isoform.autosome.gff3'
x_gff='allacma_fusca.longest_isoform.x.gff3'

p = fd.Parser(
    n=8,
    vcf=vcf,
    fasta=fasta,
    gff=a_gff,
    annotations=[
        fd.DegeneracyAnnotation()
    ],
    stratifications=[fd.DegeneracyStratification()]
)

Ok fastdfe not going to work - use easySFS with dfe alpha






#Is it better to bootstrap by site or bootstrap by gene?
#With or without replacement?
#With replacement doesnt work this way since we cant take duplicate regions out of the vcf
#We would need to do this with a python script during calculate pi bed
#Maybe i need a nextflow pipeline that will do the resampling with calculate pi bed?

bsub -G team360 -Is -n 16 -M 10240 -R "select[mem>10240] rusage[mem=10240]" bash -l
conda activate easySFS
module load bcftools/1.20--h8b25389_0
module load bedtools/2.31.1--hf5e1c6e_1

## bootstrap with replacement given 1 gff (A and X separately) to get confidence intervals
# create beds of A, X for 0D and 4D sites
agat_convert_sp_gff2bed.pl -gff allacma_fusca.longest_isoform.autosome.gff3 -o tmp.bed --sub gene && cat tmp.bed | grep -v 'agat' > allacma_fusca.longest_isoform.autosome.gene.bed
agat_convert_sp_gff2bed.pl -gff allacma_fusca.longest_isoform.x.gff3 -o tmp.bed --sub gene && cat tmp.bed | grep -v 'agat' > allacma_fusca.longest_isoform.x.gene.bed
bedtools intersect -a allacma_fusca.longest_isoform.autosome.gene.bed -b allacma_fusca.longest_isoforms.0D.bed > allacma_fusca.longest_isoform.autosome.gene.0D.bed
bedtools intersect -a allacma_fusca.longest_isoform.x.gene.bed -b allacma_fusca.longest_isoforms.0D.bed > allacma_fusca.longest_isoform.x.gene.0D.bed
bedtools intersect -a allacma_fusca.longest_isoform.autosome.gene.bed -b allacma_fusca.longest_isoforms.4D.bed > allacma_fusca.longest_isoform.autosome.gene.4D.bed
bedtools intersect -a allacma_fusca.longest_isoform.x.gene.bed -b allacma_fusca.longest_isoforms.4D.bed > allacma_fusca.longest_isoform.x.gene.4D.bed

#gene list
awk '{print $4'} allacma_fusca.longest_isoform.autosome.gene.bed > autosome_genes.txt
awk '{print $4'} allacma_fusca.longest_isoform.x.gene.bed > x_genes.txt

# create a set of random gene samples (need to intersect beds with 4D or 0D sites)
NA=$(cat autosome_genes.txt | wc -l)
NX=$(cat x_genes.txt | wc -l)
for i in {1..100}
do
	#--repeat to sample with replacement
	shuf --repeat -n $NA autosome_genes.txt > gene_bootstrap_selection/autosomes/$i.genes.txt
	shuf --repeat -n $NX x_genes.txt > gene_bootstrap_selection/x/$i.genes.txt
done

# for each bed create a subsampled vcf, get the sfs for it, replace the first count with empty count
# then format it for dfe alpha

bcftools view -R subsample.bed -o bootstrap.vcf allacma_fusca.hard_filtered.sorted.vcf.gz
./easySFS.py -i bootstrap.vcf -p pops.txt -a --proj 24 -o bootstrap_n_sfs

## bootstrap taking X sized datasets from total gff without replacement to test for faster x
#combine X and A master beds
# X dataset size = 5625 genes
cat autosome_genes.txt x_genes.txt > all_genes.txt
cat allacma_fusca.longest_isoform.autosome.gene.0D.bed allacma_fusca.longest_isoform.x.gene.0D.bed > allacma_fusca.longest_isoform.all.gene.0D.bed
cat allacma_fusca.longest_isoform.autosome.gene.4D.bed allacma_fusca.longest_isoform.x.gene.4D.bed > allacma_fusca.longest_isoform.all.gene.4D.bed

for i in {1..100}
do
	shuf -n 5625 all_genes.txt > alpha_test/gene_selection/$i.genes.txt
done

# This seems like a really slow way to do it, maybe i will have to modify calculate pi bed for just doing SFS stuff
# Maybe i will have to modify calculate pi bed after all
for FILE in ./alpha_test/gene_selection/*
do
	grep -f $FILE allacma_fusca.longest_isoform.all.gene.0D.bed > alpha_test/wdir/0D.bed
	grep -f $FILE allacma_fusca.longest_isoform.all.gene.4D.bed > alpha_test/wdir/4D.bed
	bcftools view -R alpha_test/wdir/0D.bed -o alpha_test/wdir/0D.vcf allacma_fusca.hard_filtered.sorted.vcf.gz
	bcftools view -R alpha_test/wdir/4D.bed -o alpha_test/wdir/4D.vcf allacma_fusca.hard_filtered.sorted.vcf.gz
	./easySFS.py -i alpha_test/wdir/0D.vcf -p pops.txt -a --proj 24 -o alpha_test/wdir/0D_sfs
	./easySFS.py -i alpha_test/wdir/4D.vcf -p pops.txt -a --proj 24 -o alpha_test/wdir/4D_sfs
done

# Alternate option
# Get biallelic allele counts for each chrom for each degeneracy
# Make an autosomal file, x file, and total file
# Write a script which takes a file, a sample size N, and calculates N SFSs
# Can sample form autosome and x separately for total dataset size with replacement to get confidence intervals. 
# invariant site count should stay the same for this because im not resampling from all sites just variant positions
# To do test need to sample X size from everything to generate SFSs. If sampling
# Double check where i got invariant site count from? 0D and 4D bed file length minus variant site count?
# I think this becomes difficult to figure out with this kind of test and maybe the by gene test eith easysfs is better for this.
# Could open a window to run this in a tmux session while other stuff runs?
# Get confidence intervals first though then do this kind of test afterwords, thats the order of priority

# callable invariant bed from resamples more complicated - figure that out second

# Figure out if need 1 sfs file or A and X sfs file together
# Must be both right otherwise how do i have diff X and A dfe file?
# Should just run the whole dfe alpha pipeline on each file independently 

# I have dir of full SFSs for A and X and resamples
# est demography with either autosomes, summed sfs, or combined sfs file
# loop to est sel dfe with full SFSs for A and X and then each resample
# get alpha and omega alpha from each loop
# get prop mut dist for each loop
# Give alpha and omega alpha and dfe plot confidence intervals
# if this isnt good enough can resample by gene after christmas

# bootstrapping by snp feels very conservative cos the SFS doesnt vary much
# in the long run can carry on with the above bootstrapping by gene to do cross validation test
# and maybe bootstrapping by gene
# I guess it will do for now

module load dfe_alpha/2.16-c1
cd /lustre/scratch126/tol/teams/jaron/projects/springtails_haploid_selection/data/results/dfe_resampling/dfe_alpha_resample/SCOs/all

For a set of genes with allele counts from gene pop wdir 
cat <(tail -n +2 OX359245.1.0D_pi.biallelic_ac.txt) <(tail -n +2 OX359246.1.0D_pi.biallelic_ac.txt) <(tail -n +2 OX359247.1.0D_pi.biallelic_ac.txt) <(tail -n +2 OX359248.1.0D_pi.biallelic_ac.txt) > A.0D.tsv
cat <(tail -n +2 OX359245.1.4D_pi.biallelic_ac.txt) <(tail -n +2 OX359246.1.4D_pi.biallelic_ac.txt) <(tail -n +2 OX359247.1.4D_pi.biallelic_ac.txt) <(tail -n +2 OX359248.1.4D_pi.biallelic_ac.txt) > A.4D.tsv
cat <(tail -n +2 OX359249.1.0D_pi.biallelic_ac.txt) <(tail -n +2 OX359250.1.0D_pi.biallelic_ac.txt) > X.0D.tsv
cat <(tail -n +2 OX359249.1.4D_pi.biallelic_ac.txt) <(tail -n +2 OX359250.1.4D_pi.biallelic_ac.txt) > X.4D.tsv

# intersect other beds with these beds to get subsets
# module load bedtools/2.31.1--hf5e1c6e_1 
# bedtools intersect -a data/results/gene_pop_male_bias_sco/allacma_fusca.longest_isoforms.0D.bed -b data/workdir/dfe_alpha/0D_callable.invariant.sco.bed | cut -f-1 -d'.' | uniq -c
cat data/workdir/dfe_alpha/0D_callable.invariant.sco.bed | cut -f-1 -d'.' | uniq -c
0D_A = 5746147
4D_A = 4214009
cat data/workdir/dfe_alpha/4D_callable.invariant.sco.bed | cut -f-1 -d'.' | uniq -c
0D_X = 1283252
4D_X = 947167

# look at that meisel paper on syn demography problems
# settled on autosomal intergenic SFS for neutral demographic inference
cd cd /data/tol/teams/jaron/lustre/projects/springtail_haploid_selection/data/results/dfe_alpha
est_dfe -c est_dfe_config_file_neut.txt
#2=  74 N1 100 N2 5 t2 56.0321 Nw  5.35 f0 0.985579798 L -9711857.0059
#3= 264 N1 100 N2 2 t2 11.9184 N3 5 t3 5.0000 Nw  6.20 f0 0.985488850 L -9711187.4713

#2*(-9711187.4713--9711857.0059)=1339 use 3 epoch

#*.all.sfs.txt or *.SCO.sfs.txt
cp inputs/A.SCO.sfs.txt wdir/sfs.txt
est_dfe -c est_dfe_config_file_sel.txt
prop_muts_in_s_ranges -c wdir/results_dir_sel/est_dfe.out -o wdir/prop_muts.txt
awk '{print $3"     "$6"    "$9"    "$12}' wdir/prop_muts.txt > results/obs_A_dfe.txt
est_alpha_omega -c est_alpha_omega_config_file.txt
cat wdir/est_alpha_omega.out > results/obs_A_alpha.txt

cp inputs/X.SCO.sfs.txt wdir/sfs.txt
est_dfe -c est_dfe_config_file_sel.txt
prop_muts_in_s_ranges -c wdir/results_dir_sel/est_dfe.out -o wdir/prop_muts.txt
awk '{print $3"     "$6"    "$9"    "$12}' wdir/prop_muts.txt > results/obs_X_dfe.txt
est_alpha_omega -c est_alpha_omega_config_file.txt
cat wdir/est_alpha_omega.out > results/obs_X_alpha.txt

#A
# Neutral divergence 0.678827
# Selected divergence 0.031225
# Fixation prob of deleterious mutation 0.0110267481
# adaptive_divergence -0.061592 alpha -1.972515 omega_A -0.090733

# X
# Neutral divergence 0.678827
# Selected divergence 0.031225
# Fixation prob of deleterious mutation 0.0011096106
# adaptive_divergence 0.021885 alpha 0.700879 omega_A 0.032239

for FILE in ./inputs/A/*
do
	cp $FILE wdir/sfs.txt
	est_dfe -c est_dfe_config_file_sel.txt
	est_alpha_omega -c est_alpha_omega_config_file.txt
	cat wdir/est_alpha_omega.out | cut -f5- -d' ' >> results/bootstrap_A_alpha.txt
	prop_muts_in_s_ranges -c wdir/results_dir_sel/est_dfe.out -o wdir/prop_muts.txt
	awk '{print $3"     "$6"    "$9"    "$12}' wdir/prop_muts.txt >> results/bootstrap_A_dfe.txt
done

for FILE in ./inputs/X/*
do
	cp $FILE wdir/sfs.txt
	est_dfe -c est_dfe_config_file_sel.txt
	est_alpha_omega -c est_alpha_omega_config_file.txt
	cat wdir/est_alpha_omega.out | cut -f5- -d' ' >> results/bootstrap_X_alpha.txt
	prop_muts_in_s_ranges -c wdir/results_dir_sel/est_dfe.out -o wdir/prop_muts.txt
	awk '{print $3"     "$6"    "$9"    "$12}' wdir/prop_muts.txt >> results/bootstrap_X_dfe.txt
done