params.softmasked_genome = '/lustre/scratch126/tol/teams/jaron/projects/springtails_haploid_selection/0_data/qeAllFusc8.1_refrence/genome.fasta.gz'
params.rnaseq_alignments = '/lustre/scratch126/tol/teams/jaron/projects/springtails_haploid_selection/3_hisat/bams/*.sorted.bam'
params.species = 'Allacma_fusca'
params.outdir = '/lustre/scratch126/tol/teams/jaron/projects/springtails_haploid_selection/4_braker/results'

log.info """\
         R N A S E Q   N F   P I P E L I N E    
         ===================================
         genome       : ${params.softmasked_genome}
         bams        : ${params.rnaseq_alignments}
         outdir       : ${params.outdir}
         """
         .stripIndent()

include { braker_flow } from '/lustre/scratch126/tol/teams/jaron/projects/springtails_haploid_selection/4_braker/braker_flows.nf'

workflow {
	alignment_ch = Channel .fromPath(params.rnaseq_alignments, checkIfExists:true )
	braker_flow(params.softmasked_genome, alignment_ch.collect().map { it.join(',') }, params.species)
}
