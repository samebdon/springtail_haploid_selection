params.bams = '/lustre/scratch126/tol/teams/jaron/projects/springtails_haploid_selection/3_hisat/bams/*.sorted.bam'
bam_channel = Channel .fromPath(params.bams, checkIfExists:true )
bam_list =  bam_channel.collect().map { it.join(',') }
bam_list.view()
