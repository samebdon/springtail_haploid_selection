include {degenotate; } from './gene_pop_tasks.nf'

workflow gene_pop_flow {
        take:
        	genome
        	vcf 
        	annotation
        	species
        main:
        	degenotate(genome, annotation, species)

}