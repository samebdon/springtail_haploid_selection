include {makeGenomeFile; getGeneBedGTF; getGeneBedGFF; splitBed; degenotate; filterBed; subsetVCF; calculatePiBed; joinPi; concat_all} from './gene_pop_tasks.nf'

workflow gene_pop_flow_GTF {
        take:
        	genome
        	vcf 
        	vcf_index
        	annotation
        	species
        main:
        	getGeneBedGTF(annotation, species)
        	splitBed(getGeneBedGTF.out)
        	degenotate(genome, annotation, species)
        	//filterBed(degenotate.out)
        	//subsetVCF(filterBed.out, vcf)
        	//calculatePiBed(subsetVCF.out, filterBed.out, Channel.of(splitBed.out))
        	//joinPi(calculatePiBed.out)
        	//concat_all(joinPi.out.collect())
}

workflow gene_pop_flow_GFF {
        take:
        	genome
        	genome_dict
        	vcf 
        	vcf_index
        	annotation
        	species
        main:
        	makeGenomeFile(genome_dict, species)
        	getGeneBedGFF(annotation, species)
        	splitBed(getGeneBedGFF.out)
        	bam_ch = splitBed.out.flatten()
        	degenotate(genome, annotation, species)
        	filterBed(degenotate.out.degen, degenotate.out.longest_isoforms)
        	// subsetVCF(filterBed.out, vcf, vcf_index)
        	calculatePiBed(vcf, vcf_index, filterBed.out, bam_ch, makeGenomeFile.out)
        	joinPi(calculatePiBed.out)
        	concat_all(joinPi.out.collect(), species)
}