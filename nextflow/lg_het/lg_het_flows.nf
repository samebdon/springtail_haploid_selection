include { sortBam; addRG; markDupes;  indexBam; mosdepth; freebayes; bcftools; getHet; alleleCounter ; gatk_aseReadCount } from './lg_het_tasks.nf'

workflow lg_het_flow {
        take:
	  genome
	  genome_index
	  genome_dict
	  bam_files
        main:
	  sortBam(bam_files)
	  addRG(sortBam.out)
          markDupes(addRG.out)
	  indexBam(markDupes.out)
          mosdepth(markDupes.out.join(indexBam.out))
	  freebayes(genome, markDupes.out.join(indexBam.out).join(mosdepth.out))
	  bcftools(freebayes.out)
	  getHet(bcftools.out.join(mosdepth.out))
	  alleleCounter(genome, bcftools.out.join(markDupes.out).join(indexBam.out))
	  // gatk_aseReadCount(genome, genome_index, genome_dict, bcftools.out.join(markDupes.out).join(indexBam.out))
}
