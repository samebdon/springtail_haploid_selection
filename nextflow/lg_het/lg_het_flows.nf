/* 
 * include requires tasks 
 */
include { markDupes; sortBam; indexBam; mosdepth; freebayes; bcftools  } from './lg_het_tasks.nf'

/* 
 * define the data analysis workflow 
 */

workflow lg_het_flow {
        // required inputs
        take:
	  genome
	  bam_files
        // workflow implementation
        main:
          markDupes(bam_files)
	  sortBam(markDupes.out)
	  indexBam(sortBam.out)
          mosdepth(sortBam.out, indexBam.out)
	  freebayes(genome, sortBam.out, indexBam.out, mosdepth.out)
	  bcftools(freebayes.out)
}
