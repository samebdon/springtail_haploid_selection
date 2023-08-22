params.genome = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.fasta"
params.genome_index = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.fasta.fai"
params.genome_dict = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.dict"
params.reads = "$launchDir/data/raw_data/reseq/WW1-2/LR09_EDSW200011406-1a_HJ5JVDSXY_L2-trimmed-pair{1,2}.fastq.gz"
params.outdir = "$launchDir/data/results/var_call/allacma_fusca/"
params.species = "allacma_fusca"

log.info """\
         V A R  C A L L   N F   P I P E L I N E    
         ===================================
         genome : ${params.genome}
         genome index : ${params.genome_index}
         genome dict : ${params.genome_dict}
	 reads : ${params.reads}
         outdir : ${params.outdir}
         species : ${params.species}
         """
         .stripIndent()

include { var_call_flow } from './var_call_flows.nf'

workflow {
        read_pairs_ch = Channel.fromFilePairs( params.reads, checkIfExists:true )
        var_call_flow(params.genome, params.genome_index, params.genome_dict, read_pairs_ch, params.species)
}