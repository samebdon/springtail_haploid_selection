process herhoTally {
	cpus 6
	memory '40G'
	queue 'long'

        input:
        path(vcf)
	path(vcf_index)
	path(bed)
	val(species)

        output:
        tuple val(species), path("${species}.heRho_tally_per_chromosome.tsv")

        script:
        """
        python /software/team360/heRho/heRho/heRho_tally_pairwise_counts_vcf.py -v ${vcf} -b ${bed} -d 25000 -t ${task.cpus} -c OX359245.1,OX359246.1,OX359247.1,OX359248.1,OX359249.1,OX359250.1 -f ${species}
        """
}

process herhoStandAlone {

        input:
        tuple val(species), path(herho_tally)

        output:
        path("")

        script:
        """
        python /software/team360/heRho/heRho/heRho_stand_alone.py -i ${herho_tally} > ${species}.herho_results.txt
        """
}
