params.genome = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.fasta"
params.reads = "$launchDir/data/raw_data/rnaseq/AF_[MF]_*/*.{1,2}.fastq.gz"
params.outdir = "$launchDir/data/results/rnaseq_aln/allacma_fusca/"

log.info """\
         R N A S E Q   N F   P I P E L I N E    
         ===================================
         genome	: ${params.genome}
         reads	: ${params.reads}
         outdir	: ${params.outdir}
         """
         .stripIndent()

include { raw_qc_flow; rnaseq_aln_flow } from './rnaseq_aln_flows.nf'

workflow {
        read_pairs_ch = Channel .fromFilePairs( params.reads, checkIfExists:true )
        raw_qc_flow( read_pairs_ch)
        rnaseq_aln_flow(params.genome, read_pairs_ch)
}