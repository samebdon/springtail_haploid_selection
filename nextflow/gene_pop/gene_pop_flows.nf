include {getGeneBed; splitBed; degenotate; filterBed; subsetVCF; calculatePiBed; join; concat_all} from './gene_pop_tasks.nf'

workflow gene_pop_flow {
        take:
        	genome
        	vcf 
        	annotation
        	species
        main:
        	// getGeneBed(annotation, species)
        	// splitBed(annotation.out)
        	degenotate(genome, annotation, species)
        	// filterBed(degenotate.out)
        	// subsetVCF(filterBed.out, vcf)
        	// calculatePiBed(subsetVCF.out, Channel.of(splitBed.out))
        	// join(calculatePiBed.out)
        	// concat_all(join.out.collect())
}