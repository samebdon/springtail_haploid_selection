process featureCounts {
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(meta), path(bam_f)
	path(annotation_f)

        output:
        tuple val(meta), path("featureCounts/${meta}.featureCounts.txt")

        script:
        """
        mkdir featureCounts
        featureCounts -p -B -T 4 -a ${annotation_f} -o featureCounts/${meta}.featureCounts.txt ${bam_f}
        """
}

