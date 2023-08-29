process degenotate {
        publishDir params.outdir, mode:'copy'

        input:
        path(genome_f)
        path(annotation_f)
        path(species)

        output:
        tuple val(species), path("degenotate/*")

        script:
        """
        mkdir degenotate
        degenotate.py -g ${genome_f} -a ${annotation_f} -o degenotate -c ${species}.cds-nt.fa -ca ${species}.cds-aa.fa -l ${species}.cds-nt-longest.fa -la ${species}.cds-aa-longest.fa
        """
}