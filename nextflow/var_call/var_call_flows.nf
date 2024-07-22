include { bwaIndex; bwaMem; sortBamSambamba; markDupesSambamba; indexBamSambamba; mosdepth; intersectBeds; sambambaMerge; freebayes; freebayesParallel; bcftools_filter; generate_fail_bed; generate_pass_vcf; bedtools_subtract; bcftools_sort; bcftools_index} from './var_call_tasks.nf'

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
	  sortBamSambamba(bwaMem.out)
          markDupesSambamba(sortBamSambamba.out)
	  indexBamSambamba(markDupesSambamba.out.meta_bam)
          mosdepth(markDupesSambamba.out.meta_bam.join(indexBamSambamba.out), 8, 4)
          intersectBeds(mosdepth.out.collect(), species)
	  sambambaMerge(markDupesSambamba.out.bam_only.collect(), species)
	  freebayes(genome, genome_index, sambambaMerge.out, intersectBeds.out.overlap)
	  bcftools_filter(genome, freebayes.out)
	  generate_fail_bed(bcftools_filter.out)
	  generate_pass_vcf(bcftools_filter.out)
	  bedtools_subtract(intersectBeds.out.overlap, generate_fail_bed.out)
	  bcftools_sort(generate_pass_vcf.out)
	  bcftools_index(bcftools_sort.out)
}

workflow filter_vcf {

	take:
	genome
	vcf
	callable_bed

	main:
	bcftools_filter(genome, vcf)
	generate_fail_bed(bcftools_filter.out)
	generate_pass_vcf(bcftools_filter.out)
	bedtools_subtract(callable_bed, generate_fail_bed.out)
	bcftools_sort(generate_pass_vcf.out)
	bcftools_index(bcftools_sort.out)
}
