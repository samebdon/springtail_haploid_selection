include { bwaIndex; bwaMem; sortBam; markDupes;  indexBam; mosdepth; bedtoolsIntersect; freebayes; bcftools; } from './var_call_tasks.nf'
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
          callables = mosdepth.out.collect().first()
          bedtoolsIntersect(callables.last(), callables.until(callables.last()), species)
	  samtoolsMerge(markDupes.out.join(indexBam.out).collect(), species)
	  indexMergedBam(samtoolsMerge.out)
	  freebayes(genome, samtoolsMerge.out.join(indexMergedBam.out), bedtoolsIntersect.out)
	  bcftools(freebayes.out)
}
