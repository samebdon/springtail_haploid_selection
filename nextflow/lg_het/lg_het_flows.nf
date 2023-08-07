include { markDupes; sortBam; indexBam; mosdepth; freebayes; bcftools; getHet; alleleCounter } from './lg_het_tasks.nf'

workflow lg_het_flow {
        take:
	  genome
	  bam_files
        main:
          markDupes(bam_files)
	  sortBam(markDupes.out)
	  indexBam(sortBam.out)
          mosdepth(sortBam.out.join(indexBam.out))
	  freebayes(genome, sortBam.out.join(indexBam.out).join(mosdepth.out))
	  bcftools(freebayes.out)
	  getHet(bcftools.out.join(mosdepth.out))
	  alleleCounter(genome, bcftools.out.join(sortBam.out).join(indexBam.out))
}
