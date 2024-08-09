include {getLongestIsoformAGAT; makeGenomeFile; getGeneBed; splitBed; degenotate; filterBed; subsetVCF; calculatePiBed; mergePi; concat_all; concat_SFS} from './gene_pop_tasks.nf'

workflow gene_pop_flow_SFS {
        take:
                genome
                genome_dict
                vcf 
                vcf_index
                annotation
                species
        main:
                getLongestIsoformAGAT(annotation)
                makeGenomeFile(genome_dict, species)
                getGeneBed(getLongestIsoformAGAT.out, species)
                splitBed(getGeneBed.out)
                bed_ch = splitBed.out.flatten()
                degenotate(genome, getLongestIsoformAGAT.out, species)
                filterBed(degenotate.out.degen, degenotate.out.longest_isoforms)
                calculatePiBed(vcf, vcf_index, filterBed.out, bed_ch, makeGenomeFile.out)
                mergePi(calculatePiBed.out.pi)
                concat_all(mergePi.out.collect(), species)
                concat_SFS(calculatePiBed.out.sfs.collect(), species)
        emit:
                concat_all.out
}