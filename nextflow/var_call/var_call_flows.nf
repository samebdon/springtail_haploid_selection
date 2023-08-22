include { bwaIndex; bwaMem; sortBam; markDupes; indexBam; mosdepth; bedtoolsIntersect; samtoolsMerge; freebayes; bcftools_filter; generate_fail_bed; generate_pass_vcf; bedtools_subtract; bcftools_sort; bcftools_index} from './var_call_tasks.nf'
include { indexBam as indexMergedBam } from './var_call_tasks.nf'

workflow var_call_flow {
        take:
	  genome
	  genome_index
	  genome_dict
	  read_files
	  species
        main:
          bwaIndex(genome)
          bwaMem(genome, bwaIndex.out, read_files)
	  sortBam(bwaMem.out)
          markDupes(sortBam.out)
	  indexBam(markDupes.out)
          mosdepth(markDupes.out.join(indexBam.out))
          callables = mosdepth.out.collect()
          bedtoolsIntersect(callables.last(), callables.until(callables.last()), species)
	  samtoolsMerge(markDupes.out.join(indexBam.out).collect(), species)
	  indexMergedBam(samtoolsMerge.out)
	  freebayes(genome, samtoolsMerge.out.join(indexMergedBam.out), bedtoolsIntersect.out)
	  bcftools_filter(freebayes.out)
	  generate_fail_bed(bcftools_filter.out)
	  generate_pass_vcf(bcftools_filter.out)
	  bedtools_subtract(bedtoolsIntersect.out, generate_fail_bed.out)
	  bcftools_sort(generate_pass_vcf.out)
	  bcftools_index(generate_pass_vcf.out)
}

