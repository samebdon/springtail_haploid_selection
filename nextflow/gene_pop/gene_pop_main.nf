params.genome = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.fasta"
params.genome_index = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.fasta.fai"
params.genome_dict = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.dict"
params.vcf = "$launchDir/data/results/lg_het/allacma_fusca/vcfs/AF_M_6.coord_sorted.RG.deduped.GCA_947179485.1.simple_header.sorted.filtered.vcf.gz"
params.annotation = "$launchDir/data/results/braker/allacma_fusca/hint_support/tsebra.augustus.gtf"
params.outdir = "$launchDir/data/results/var_call/allacma_fusca/"
params.species = "allacma_fusca"

log.info """\
         V A R  C A L L   N F   P I P E L I N E    
         ===================================
         vcf : ${params.vcf}
	 annotation : ${params.annotation}
         outdir : ${params.outdir}
         species : ${params.species}
         """
         .stripIndent()

include { gene_pop_flow } from './gene_pop_flows.nf'

workflow {
        gene_pop_flow(params.genome, params.vcf, params.annotation, params.species)
}
