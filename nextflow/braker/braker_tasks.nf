process braker {
	publishDir params.outdir, mode:'move'

	input:
	path genome
	val rnaseq_alignments
	val species

	output:
	path *

	script:
	"""
	braker.pl --genome=${genome} --bam=${rnaseq_alignments} --species=${species}
	"""
}
