# I need to include 0 sites in gene pop output...
# im not sure if my 0 counts are correct... these definitely feel wrong
# How do i get a good count of invariant 0D and 4D positions using my approach?
# Can i intersect the 0D and 4D beds with the VCF and get the leftover file per chrom, and use this number?
# Yeah, invariant 0D and 4D counts on X and A. Pass count as an option 
# do i need to run X and A separately or does it assume all the SFSs are for the same thing?
# I should use SFS from SCOs because it uses dn and ds too
# they recommend to remove polymorphism for est alpha omega as described in 2012 paper
# divergence is large so not likely to be a problem but remember for later

module load bedtools/2.31.1--hf5e1c6e_1

parallel -j1 'grep {} data/results/ortholog_pop_gen/allacma_fusca.vs.dicyrtomina_minuta/allacma_fusca.cds.bed' :::: data/results/ortholog_pop_gen/allacma_fusca.vs.dicyrtomina_minuta/allacma_fusca.sp1.SCO_genes.txt > data/workdir/dfe_alpha/sco.bed
bedtools intersect -a data/workdir/dfe_alpha/sco.bed -b data/results/var_call/allacma_fusca_rm_repeats/allacma_fusca.callable.freebayes.norepeats.bed > data/workdir/dfe_alpha/sco_callable.bed 
bedtools intersect -a data/results/ortholog_pop_gen/allacma_fusca.vs.dicyrtomina_minuta/allacma_fusca.longest_isoforms.0D.bed -b data/workdir/dfe_alpha/sco_callable.bed > data/workdir/dfe_alpha/0D_callable.sco.bed
bedtools intersect -a data/results/ortholog_pop_gen/allacma_fusca.vs.dicyrtomina_minuta/allacma_fusca.longest_isoforms.4D.bed -b data/workdir/dfe_alpha/sco_callable.bed > data/workdir/dfe_alpha/4D_callable.sco.bed
bedtools subtract -a data/workdir/dfe_alpha/0D_callable.sco.bed -b data/results/var_call/allacma_fusca_rm_repeats/allacma_fusca.hard_filtered.sorted.vcf.gz > data/workdir/dfe_alpha/0D_callable.invariant.sco.bed
bedtools subtract -a data/workdir/dfe_alpha/4D_callable.sco.bed -b data/results/var_call/allacma_fusca_rm_repeats/allacma_fusca.hard_filtered.sorted.vcf.gz > data/workdir/dfe_alpha/4D_callable.invariant.sco.bed

# SFS is from SCOs
# Should I be generating the DFE from all genes not just SCOs
# even if i have to use SCOs for divergence information?
# Would just give more information for the DFE...

#sco
cat 0D_callable.invariant.sco.bed | cut -f-1 -d'.' | uniq -c
0D_A = 5746147
4D_A = 4214009
cat 4D_callable.invariant.sco.bed | cut -f-1 -d'.' | uniq -c
0D_X = 1283252
4D_X = 947167

#written config files to data/results/dfe_alpha/
#config file neutral est_dfe
data_path_1 /lustre/scratch126/tol/teams/jaron/users/sam/dfe_alpha_extra_files/data_path_1/
data_path_2 /lustre/scratch126/tol/teams/jaron/users/sam/dfe_alpha_extra_files/data_path_2/
sfs_input_file data/results/dfe_alpha/sfs.txt
est_dfe_results_dir data/results/dfe_alpha/results_dir_neut
site_class 0
fold 1
epochs 2
search_n2 1
t2_variable 1
t2 50

#config file selected est_dfe
data_path_1 /lustre/scratch126/tol/teams/jaron/users/sam/dfe_alpha_extra_files/data_path_1/
data_path_2 /lustre/scratch126/tol/teams/jaron/users/sam/dfe_alpha_extra_files/data_path_2/
sfs_input_file /lustre/scratch126/tol/teams/jaron/projects/springtails_haploid_selection/data/results/dfe_alpha/sfs.txt
est_dfe_results_dir data/results/dfe_alpha/results_dir_sel
est_dfe_demography_results_file data/results/dfe_alpha/results_dir_neut/est_dfe.out
site_class 1
fold 1
epochs 2
mean_s_variable 1
mean_s -0.1
beta_variable 1
beta 0.5

#config file est_alpha_omega
data_path_1 /lustre/scratch126/tol/teams/jaron/users/sam/dfe_alpha_extra_files/data_path_1/
divergence_file data/results/dfe_alpha/divergence_file.txt
est_alpha_omega_results_file data/results/dfe_alpha/results_dir_alpha_omega/est_alpha_omega.out
est_dfe_results_file data/results/dfe_alpha/results_dir_sel/est_dfe.out
neut_egf_file data/results/dfe_alpha/results_dir_neut/neut_egf.out
sel_egf_file data/results/dfe_alpha/results_dir_sel/sel_egf.out
do_jukes_cantor 1
remove_poly 1

#sfs.txt including A (top) and X (bottom) in one
# I think i probably need to split it separately and do different analysis on X and A
# But I could run together to see what happens?
2
24
2586 16633 7557 5467 5113 4682 4298 4220 4266 3942 3431 3898 3665 0 0 0 0 0 0 0 0 0 0 0 0 0
2482 13417 6680 4831 4334 3825 3807 3489 3800 3389 3161 3493 2731 0 0 0 0 0 0 0 0 0 0 0 0 0
23
172 1173 520 337 278 298 183 117 170 155 142 231 364 0 0 0 0 0 0 0 0 0 0 0 0 0
187 597 521 344 296 201 198 161 198 153 153 161 251 0 0 0 0 0 0 0 0 0 0 0 0 0

#divergence_file.txt
#1,total number of 0d sites, number of 4d fixed differences
#2,number of 4d sites, number of 4d fixed differences
# should be able to get these counts from ortholog_pop_gen
1 6588078 253571 
0 1345742 813059

module load dfe_alpha/2.16-c1

est_dfe -c data/results/dfe_alpha/est_dfe_config_file_neut_A.txt
est_dfe -c data/results/dfe_alpha/est_dfe_config_file_neut_X.txt
est_dfe -c data/results/dfe_alpha/est_dfe_config_file_sel_A.txt
est_dfe -c data/results/dfe_alpha/est_dfe_config_file_sel_X.txt
est_alpha_omega -c data/results/dfe_alpha/est_alpha_omega_config_file_A.txt
est_alpha_omega -c data/results/dfe_alpha/est_alpha_omega_config_file_X.txt


