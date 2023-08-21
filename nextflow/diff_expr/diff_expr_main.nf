params.bams = "$launchDir/data/results/rnaseq_aln/allacma_fusca/hisat2/bams/*"
params.annotation = "$launchDir/data/results/braker/allacma_fusca/hint_support/tsebra.augustus.gtf"
params.outdir = "$launchDir/data/results/diff_expr/allacma_fusca/"

log.info """\
         DIFF EXPR   N F   P I P E L I N E    
         ===================================
         bams : ${params.bams}
	 annotation : ${params.annotation}
         outdir : ${params.outdir}
         """
         .stripIndent()

include { diff_expr_flow } from './diff_expr_flows.nf'

workflow {
        bam_ch = Channel.fromPath( params.bams, checkIfExists:true ).map { file -> tuple(file.simpleName, file) }
        diff_expr_flow(bam_ch, params.annotation)
}
