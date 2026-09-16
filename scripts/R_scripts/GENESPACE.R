#agat_convert_sp_gff2bed.pl --gff braker.gtf -o braker.agat.bed
#cut -f-4 -d'	' braker.agat.bed > braker.agat.trimmed.bed

'./data/results/braker3/allacma_fusca/braker.agat.bed'
'./data/results/braker3/allacma_fusca/braker.aa'

'./data/results/braker2/dicyrtomina_minuta/braker.agat.bed'
'./data/results/braker2/dicyrtomina_minuta/braker.aa'

'./data/results/braker2/sminthurides_aquaticus/braker.agat.bed'
'./data/results/braker2/sminthurides_aquaticus/braker.aa'


# Running Genespace

bsub -G team360 -Is -n 16 -M 10240 -R "select[mem>10240] rusage[mem=10240]" bash -l

module load mcscanx/0.0.1-c1
conda activate /software/treeoflife/conda/users/envs/team360/se13/genespace

R

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("jtlovell/GENESPACE")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "rtracklayer"))

library(GENESPACE)

## wd should include /bed/ and /peptide/ directories with the bed files and aa files respectively 
## 4 column bed files, column 4 matches pep fasta headers
wd <- "./data/workdir/genespace" 
path2mcscanx <- "/software/treeoflife/shpc/0.1.26/wrapper/quay.io/sanger-tol/mcscanx/0.0.1-c1/bin/"
path2orthofinder <- "/software/treeoflife/conda/users/envs/team360/se13/genespace/bin/orthofinder"

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx,
  path2orthofinder = path2orthofinder)
# Should find valid path to orthofinder, mcscanx, and diamond

## It can pause here and give you a command to run orthofinder outside of R then join back in again :eyeroll:
out <- run_genespace(gpar, overwrite = T)