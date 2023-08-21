include { featureCounts; } from './diff_expr_tasks.nf'

workflow diff_expr_flow {
        take:
          bam_files
	  annotation_file
        main:
	  featureCounts(bam_files, annotation_file)
}
