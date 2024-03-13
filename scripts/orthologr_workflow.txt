bsub -q normal -G team360-grp -Is -n 32 -M 40000 -R "select[mem>40000] rusage[mem=40000]" bash -l 

mamba activate orthologr
R

### INSTALLATION
# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

# Install package dependencies
BiocManager::install(c(
        "Biostrings",
        "GenomicRanges",
        "GenomicFeatures",
        "Rsamtools",
        "rtracklayer"
))

# install CRAN dependencies
install.packages(c("doParallel", "foreach", "ape", "Rdpack", "benchmarkme", "devtools"))

# install BLAST dependency metablastr from GitHub
devtools::install_github("drostlab/metablastr")

# install DIAMOND dependency rdiamond from GitHub
devtools::install_github("drostlab/rdiamond")

# install orthologr from GitHub
devtools::install_github("drostlab/orthologr")
###

library(orthologr)
# get a dNdS table using:
# 1) reciprocal best hit for orthology inference (RBH)
# 2) Needleman-Wunsch for pairwise amino acid alignments
# 3) pal2nal for codon alignments
# 4) Comeron for dNdS estimation
# 5) single core processing 'comp_cores = 1'

dnds_result = dNdS(query_file ='./data/results/braker3/allacma_fusca/afusca.braker.codingseq', subject_file = './data/results/braker2/svir_clem/svir_clem.braker.codingseq', ortho_detection = "RBH", aa_aln_type     = "pairwise",aa_aln_tool     = "NW", codon_aln_tool  = "pal2nal", dnds_est.method = "Comeron", comp_cores      = 32)

write.table(dnds_result, file = './data/results/afusca.svir_clem.braker_cds.dnds.orthologr.csv')

dnds_result = dNdS(query_file ='./data/results/braker3/allacma_fusca/braker.codingseq', subject_file = './data/results/braker2/sminthurides_aquaticus/saquaticus.braker.codingseq', ortho_detection = "RBH", aa_aln_type     = "pairwise",aa_aln_tool     = "NW", codon_aln_tool  = "pal2nal", dnds_est.method = "Comeron", comp_cores      = 32)

write.table(dnds_result, file = './data/results/afusca.saquaticus.braker_cds.dnds.orthologr.csv')

dnds_result = dNdS(query_file ='./data/results/braker2/sminthurides_aquaticus/saquaticus.braker.codingseq', subject_file = './data/results/braker2/svir_clem/svir_clem.braker.codingseq', ortho_detection = "RBH", aa_aln_type     = "pairwise",aa_aln_tool     = "NW", codon_aln_tool  = "pal2nal", dnds_est.method = "Comeron", comp_cores      = 32)

write.table(dnds_result, file = './data/results/saquaticus.svir_clem.braker_cds.dnds.orthologr.csv')


dnds_result = dNdS(
    query_file ='./data/results/toga/afusca_svir_clem/processed/one2one.codon.reference.edited.fasta', 
    subject_file = './data/results/toga/afusca_svir_clem/processed/one2one.codon.query.edited.fasta', 
    ortho_detection = "RBH", 
    aa_aln_type     = "pairwise",
    aa_aln_tool     = "NW", 
    codon_aln_tool  = "pal2nal", 
    dnds_est.method = "Comeron", 
    comp_cores      = 32)
