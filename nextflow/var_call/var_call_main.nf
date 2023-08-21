params.genome = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.fasta"
params.genome_index = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.fasta.fai"
params.genome_dict = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.dict"
params.bams = "$launchDir/data/results/rnaseq_aln/allacma_fusca/hisat2/bams/AF_M_6.query_sorted.bam"
params.outdir = "$launchDir/data/results/var_call/allacma_fusca/"
params.species = "allacma_fusca"

log.info """\
         LG HET   N F   P I P E L I N E    
         ===================================
         genome : ${params.genome}
	 bams : ${params.bams}
         outdir : ${params.outdir}
         """
         .stripIndent()

include { var_call_flow } from './var_call_flows.nf'

workflow {
        bam_ch = Channel.fromPath( params.bams, checkIfExists:true ).map { file -> tuple(file.simpleName, file) }
        var_call_flow(params.genome, params.genome_index, params.genome_dict, bam_ch, params.species)
}

