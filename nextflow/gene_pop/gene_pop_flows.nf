include {agatAnnotation; makeGenomeFile; getGeneBedGTF; getGeneBedGFF; getExonBedGFF; splitBed; degenotate; filterBed; subsetVCF; calculatePiBed; mergePi; concat_all; concat_SFS} from './gene_pop_tasks.nf'

workflow gene_pop_flow_SFS {
        take:
                genome
                genome_dict
                vcf 
                vcf_index
                annotation
                species
        main:
                agatAnnotation(annotation)
                makeGenomeFile(genome_dict, species)
                getGeneBedGFF(agatAnnotation.out, species)
                splitBed(getGeneBedGFF.out)
                bed_ch = splitBed.out.flatten()
                degenotate(genome, agatAnnotation.out, species)
                filterBed(degenotate.out.degen, degenotate.out.longest_isoforms)
                calculatePiBed(vcf, vcf_index, filterBed.out, bed_ch, makeGenomeFile.out)
                mergePi(calculatePiBed.out.pi)
                concat_all(mergePi.out.collect(), species)
                concat_SFS(calculatePiBed.out.sfs.collect(), species)
        emit:
                concat_all.out
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
        	bed_ch = splitBed.out.flatten()
        	degenotate(genome, annotation, species)
        	filterBed(degenotate.out.degen, degenotate.out.longest_isoforms)
        	// subsetVCF(filterBed.out, vcf, vcf_index)
        	calculatePiBed(vcf, vcf_index, filterBed.out, bed_ch, makeGenomeFile.out)
        	mergePi(calculatePiBed.out)
        	concat_all(mergePi.out.collect(), species)
}

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