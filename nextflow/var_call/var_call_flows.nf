include { sortBam; addRG; markDupes;  indexBam; mosdepth; bedtoolsIntersect; freebayesPop; bcftools; } from './var_call_tasks.nf'
include { indexBam; } } from './var_call_tasks.nf' as indexMergedBam
workflow lg_het_flow {
        take:
	  genome
	  genome_index
	  genome_dict
	  bam_files
	  species
        main:
	  sortBam(bam_files)
	  addRG(sortBam.out)
          markDupes(addRG.out)
	  indexBam(markDupes.out)
          mosdepth(markDupes.out.join(indexBam.out))
          callables = mosdepth.out.collect().first()
          bedtoolsIntersect(callables.last(), callables.until(callables.last()), species)
	  samtoolsMerge(markDupes.out.join(indexBam.out).collect(), species)
	  indexMergedBam(samtoolsMerge.out)
	  freebayesPop(genome, samtoolsMerge.out.join(indexMergedBam.out), bedtoolsIntersect.out)
	  bcftools(freebayes.out)
	 
}
