process herhoTally {
	cpus 5

        input:
        path(vcf)
	path(bed)
	val(species)

        output:
        tuple val(species), path("${species}.heRho_tally_per_chromosome.tsv")

        script:
        """
        python /software/team360/heRho/heRho/heRho_tally_pairwise_counts_vcf.py -v ${vcf} -b ${bed} -d 1000 -t ${task.cpus} -c OX359345.1,OX359346.1,OX359347.1,OX359348.1,OX359349.1,OX359350.1 -f ${species}
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
