process unmask_genome {

        input:
        val(meta)
        path(genome)

        output:
        tuple val(meta), path("${meta}.unmasked.fa")

        script:
        """
        awk '/^>/ {print(\$0)}; /^[^>]/ {print(toupper(\$0))}' ${genome} > ${meta}.unmasked.fa
        """
}

process earlGrey {
        memory '200G'
        cpus 64
        queue 'basement'

        input:
        tuple val(meta), path(genome)

        output:
        path("results/*/*_summaryFiles/*.filteredRepeats.bed")

        script:
        """
        earlGrey -g ${unmasked_genome} -s ${meta} -o ./results -t ${task.cpus}
        """
}

process mask_genome {

        input:
        tuple val(meta), path(genome)
        path(bed)

        output:
        tuple val(meta), path("${meta}.masked.fa")

        script:
        """
        bedtools maskfasta -soft -fi ${genome} -bed ${bed} -fo ${meta}.masked.fa
        """
}

process braker2 {
        publishDir params.outdir, mode:'copy'
        memory '40G'
        cpus 32
        queue 'long'

        input:
        tuple val(meta), path(genome)
        path(prot_seq)

        output:
        tuple val(meta), path("*")

        script:
        """
        braker.pl --genome=${genome} --softmasking --workingdir=. --threads ${task.cpus} --species=${meta} --gff3 --prot_seq=${prot_seq} --useexisting
        """
}

// if this doesnt work I could probably get it working using the singularity image