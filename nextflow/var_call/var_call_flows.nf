include { bwaIndex; bwaMem; sortBam; sortBamSambamba; markDupes; markDupesSambamba; indexBam; indexBamSambamba; mosdepth; intersectBeds; samtoolsMerge; sambambaMerge; freebayes; bcftools_filter; generate_fail_bed; generate_pass_vcf; bedtools_subtract; bcftools_sort; bcftools_index} from './var_call_tasks.nf'

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
          mosdepth(markDupes.out.join(indexBam.out), 8, 5)
          intersectBeds(mosdepth.out.collect(), species)
	  samtoolsMerge(markDupes.out.join(indexBam.out).collect(), species)
	  freebayes(genome, samtoolsMerge.out, intersectBeds.out)
	  bcftools_filter(genome, freebayes.out)
	  generate_fail_bed(bcftools_filter.out)
	  generate_pass_vcf(bcftools_filter.out)
	  bedtools_subtract(intersectBeds.out, generate_fail_bed.out)
	  bcftools_sort(generate_pass_vcf.out)
	  bcftools_index(generate_pass_vcf.out)
}

workflow var_call_flow_sambamba {
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
	  indexBamSambamba(markDupesSambamba.out)
          mosdepth(markDupesSambamba.out.join(indexBamSambamba.out), 8, 5)
          intersectBeds(mosdepth.out.collect(), species)
	  sambambaMerge(markDupesSambamba.out.collect().filter( .*\.bam$ ), species)
	  freebayes(genome, sambambaMerge.out, intersectBeds.out)
	  bcftools_filter(genome, freebayes.out)
	  generate_fail_bed(bcftools_filter.out)
	  generate_pass_vcf(bcftools_filter.out)
	  bedtools_subtract(intersectBeds.out, generate_fail_bed.out)
	  bcftools_sort(generate_pass_vcf.out)
	  bcftools_index(generate_pass_vcf.out)
}