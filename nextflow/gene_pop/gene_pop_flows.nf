include {getGeneBed; splitBed; degenotate; filterBed; subsetVCF; calculatePiBed; joinPi; concat_all} from './gene_pop_tasks.nf'

workflow gene_pop_flow {
        take:
        	genome
        	vcf 
        	annotation
        	species
        main:
        	getGeneBed(annotation, species)
        	splitBed(getGeneBed.out)
        	degenotate(genome, annotation, species)
        	filterBed(degenotate.out)
        	subsetVCF(filterBed.out, vcf)
        	calculatePiBed(subsetVCF.out, filterBed.out, Channel.of(splitBed.out))
        	joinPi(calculatePiBed.out)
        	concat_all(joinPi.out.collect())
}