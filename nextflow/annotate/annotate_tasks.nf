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
        publishDir params.outdir, mode:'copy'
        memory '200G'
        cpus 64
        queue 'basement'
        conda '/software/treeoflife/conda/users/envs/team360/se13/earlgrey'

        input:
        tuple val(meta), path(genome)

        output:
        path("./results/*_EarlGrey"), emit: all
        path("./results/*/*_summaryFiles/*.filteredRepeats.bed"), emit: repeat_bed
	tuple val(meta), path("./results/*/*_summaryFiles/*.softmasked.fasta"), emit: softmasked_genome
        
	script:
        """
        earlGrey -g ${genome} -s ${meta} -o ./results -t ${task.cpus} -d yes
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
        tuple val(meta), path("wdir/*")

        script:
        """
        mkdir wdir
        braker.pl \
                --genome=${genome} \
                --softmasking \
                --workingdir=wdir \
                --threads ${task.cpus} \
                --species=${meta} \
                --gff3 \
                --prot_seq=${prot_seq} \
                --useexisting
        """
}
