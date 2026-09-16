
#BSUB -o logs/ASE.AF_F_7.out.%J
#BSUB -e logs/ASE.AF_F_7.err.%J
#BSUB -q normal
#BSUB -n 12
#BSUB -M 20000
#BSUB -R "select[mem>20000] rusage[mem=20000]"

module load ISG/gatk/4.5.0.0
module load samtools/1.20--h50ea8bc_0

gatk SplitNCigarReads   -R /lustre/scratch126/tol/teams/jaron/projects/springtail_haploid_selection/data/genomes/GCA_947179485.1_qeAllFusc8.1_genomic.simple_header.earlgrey_masked.fna   -I /lustre/scratch126/tol/teams/jaron/projects/springtail_haploid_selection/data/results/var_call/allacma_fusca_rnaseq/AF_F_7.query_sorted.coord_sorted.bam   -O /lustre/scratch126/tol/teams/jaron/projects/springtail_haploid_selection/data/results/var_call/allacma_fusca_rnaseq/AF_F_7.split.bam

samtools sort -o /lustre/scratch126/tol/teams/jaron/projects/springtail_haploid_selection/data/results/var_call/allacma_fusca_rnaseq/AF_F_7.split.sorted.bam /lustre/scratch126/tol/teams/jaron/projects/springtail_haploid_selection/data/results/var_call/allacma_fusca_rnaseq/AF_F_7.split.bam
samtools index /lustre/scratch126/tol/teams/jaron/projects/springtail_haploid_selection/data/results/var_call/allacma_fusca_rnaseq/AF_F_7.split.sorted.bam

gatk ASEReadCounter   -R /lustre/scratch126/tol/teams/jaron/projects/springtail_haploid_selection/data/genomes/GCA_947179485.1_qeAllFusc8.1_genomic.simple_header.earlgrey_masked.fna   -I /lustre/scratch126/tol/teams/jaron/projects/springtail_haploid_selection/data/results/var_call/allacma_fusca_rnaseq/AF_F_7.split.sorted.bam   -V /lustre/scratch126/tol/teams/jaron/projects/springtail_haploid_selection/data/results/var_call/allacma_fusca_rnaseq/allacma_fusca.hard_filtered.sorted.vcf.gz   -O /lustre/scratch126/tol/teams/jaron/projects/springtail_haploid_selection/data/results/ase/AF_F_7_ase_counts.table
