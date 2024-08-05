log.info """\
         V A R  C A L L   N F   P I P E L I N E    
         ===================================
         genome : ${params.genome}
	 reads : ${params.reads}
         outdir : ${params.outdir}
         species : ${params.species}
         """
         .stripIndent()

include { var_call_flow } from './var_call_flows.nf'

workflow {
        read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists:true )
        var_call_flow(params.genome, params.genome_index, params.genome_dict, read_pairs_ch, params.repeat_bed, params.species)
}


