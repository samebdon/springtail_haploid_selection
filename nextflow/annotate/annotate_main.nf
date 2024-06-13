// MAKE SURE FASTA HEADERS ARE SIMPLE FOR EARLGREY AND BRAKER
// PIPELINE WILL REMOVE EXISTING SOFTMASKING

params.meta = "dicyrtomina_minuta"
params.genome = "$launchDir/data/results/genomes/dicyrtomina_minuta/qeDicMinu4_1.earlGrey_masked.fa"
params.prot_seq = "$launchDir/data/results/braker2/dbs/Arthropoda.allacmaFusca.fa"
params.outdir = "$launchDir/data/results/braker2/${params.meta}_nf"

log.info """\
         O R T H O L O G  P O P  G E N   N F   P I P E L I N E    
         ===================================
         species : ${params.meta}
         genome :  ${params.genome}
         prot db : ${params.prot_seq}
         outdir :  ${params.outdir}
         """
         .stripIndent()

include { braker2_flow; braker2_only_flow } from './annotate_flows.nf'

workflow {
        braker2_flow(params.meta, params.genome, params.prot_seq)
        // braker2_only_flow(params.meta, params.genome, params.prot_seq)
}