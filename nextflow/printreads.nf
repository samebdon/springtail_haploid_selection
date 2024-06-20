params.reads = "$launchDir/data/raw_data/rnaseq/AF_[MF]_*/*.{1,2}.fastq.gz"
read_pairs_ch = Channel .fromFilePairs(params.reads, checkIfExists:true)
read_pairs_ch.view()
