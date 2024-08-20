include {getLongestIsoformAGAT; makeGenomeFile; getGeneBedAGAT; splitBed; degenotate; filterBed; calculatePiBed; mergePi; concat_all; concat_SFS} from './gene_pop_tasks.nf'

workflow gene_pop_flow_SFS {
        take:
                genome
                genome_dict
                vcf 
                vcf_index
                annotation
                species
                sex_linked_contigs
        main:
                getLongestIsoformAGAT(annotation)
                makeGenomeFile(genome_dict, species)
                getGeneBedAGAT(getLongestIsoformAGAT.out, species)
                splitBed(getGeneBedAGAT.out)
                bed_ch = splitBed.out.flatten()
                degenotate(genome, getLongestIsoformAGAT.out, species)
                filterBed(degenotate.out.degen, degenotate.out.longest_isoforms)
                calculatePiBed(vcf, vcf_index, filterBed.out, bed_ch, makeGenomeFile.out)
                mergePi(calculatePiBed.out.pi)
                concat_all(mergePi.out.collect(), species)
                concat_SFS(calculatePiBed.out.sfs.collect(), species, sex_linked_contigs)
        emit:
                concat_all.out
}