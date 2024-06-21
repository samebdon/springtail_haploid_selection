include {agatAnnotation; makeGenomeFile; getGeneBedGFF; splitBed; degenotate; filterBed; calculatePiBed; mergePi; concat_all} from './gene_pop_tasks.nf'

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
}