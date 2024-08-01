include { gene_pop_flow_SFS } from './gene_pop_flows.nf'

workflow {
        gene_pop_flow_SFS(params.genome, params.genome_dict, params.vcf,params.vcf_index, params.annotation, params.species)
}
