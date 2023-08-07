params.genome = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.fasta"
params.bams = "$launchDir/data/results/rnaseq_aln/allacma_fusca/hisat2/bams/*"
params.outdir = "$launchDir/data/results/lg_het/allacma_fusca/"

log.info """\
         LG HET   N F   P I P E L I N E    
         ===================================
         genome : ${params.genome}
	 bams : ${params.bams}
         outdir : ${params.outdir}
         """
         .stripIndent()

include { lg_het_flow } from './lg_het_flows.nf'

workflow {
        bam_ch = Channel.fromPath( params.bams, checkIfExists:true ).map { file -> tuple(file.simpleName, file) }
        lg_het_flow(params.genome, bam_ch)
}

