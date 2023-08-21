include { bwaIndex; bwaMem; sortBam; markDupes;  indexBam; mosdepth; bedtoolsIntersect; freebayes; bcftools; } from './var_call_tasks.nf'
include { indexBam; } } from './var_call_tasks.nf' as indexMergedBam
workflow lg_het_flow {
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
	  addRG(sortBam.out)
          markDupes(addRG.out)
	  indexBam(markDupes.out)
          mosdepth(markDupes.out.join(indexBam.out))
          callables = mosdepth.out.collect().first()
          bedtoolsIntersect(callables.last(), callables.until(callables.last()), species)
	  samtoolsMerge(markDupes.out.join(indexBam.out).collect(), species)
	  indexMergedBam(samtoolsMerge.out)
	  freebayes(genome, samtoolsMerge.out.join(indexMergedBam.out), bedtoolsIntersect.out)
	  bcftools(freebayes.out)
}
