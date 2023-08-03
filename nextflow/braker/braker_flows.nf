
include { braker } from './braker_tasks.nf'

workflow braker_flow {
	take:
	  genome
	  rnaseq_alignments
	  species
	main:
	  braker(genome, rnaseq_alignments, species)
}
