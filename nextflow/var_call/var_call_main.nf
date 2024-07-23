// could prob just put index and dict in as processes
// need to figure out config files

params.genome = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.earlGrey_masked.fasta"
params.genome_index = "$launchDir/data/results/allacma_fusca/GCA_947179485.1.simple_header.earlGrey_masked.fasta.fai"
params.genome_dict = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.earlGrey_masked.dict"
params.reads = "$launchDir/data/raw_data/reseq/allacma_fusca/*/*.{1,2}.fastq.gz"
params.outdir = "$launchDir/data/results/var_call/allacma_fusca_sdcovmax/"
params.species = "allacma_fusca"

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
        var_call_flow(params.genome, params.genome_index, params.genome_dict, read_pairs_ch, params.species)
}


