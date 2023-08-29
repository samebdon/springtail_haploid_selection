process degenotate {
        cpus 1
        memory = '2 GB'
        publishDir params.outdir, mode:'copy'

        input:
        path(genome_f)
        path(annotation_f)
        val(species)

        output:
        tuple val(species), path("degenotate/*")

        script:
        """
        degenotate.py -g ${genome_f} -a <(cat ${annotation_f} | awk -v OFS='\\t' '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}') -o degenotate -x 04 -c ${species}.cds-nt.fa -ca ${species}.cds-aa.fa -l ${species}.cds-nt-longest.fa -la ${species}.cds-aa-longest.fa
        """
}

process filterBed{
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path("degenotate/*")

        output:
        tuple val(species), path("${species}.degeneracy.longest_isoforms.bed")

        script:
        """
        get transcript names from longest isoform file
        filter degeneracy all sites bed by longest isoforms
        output scaffold, start, end, transcript, degeneracy
        """
}

process subsetVCF{
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path(bed_f)
        path(vcf_f)

        output:
        tuple val(species), path("${species}.04D.longest_isoforms.vcf.gz")

        script:
        """
        bcftools view --threads ${task.cpus} -Oz -R <(cat ${bed_f} | cut -f1-3) -o ${species}.04D.longest_isoforms.vcf.gz ${vcf_f}
        """
}

process generateSFS{

        input:
        tuple val(species), path(vcf_f), path(bed_f)

        output:
        tuple val(species), path("${species}.04D.longest_isoforms.sfs.tsv")

        script:
        """
        custom python script i think
        load degeneracy bed, create data frame. generate gene name from transcript name
        load vcf, for each variant (should be just 0D and 4D SNPs in genes) calculate SFS or pi?
        make dataframe with chrom, pos, geneid, SFS
        merge with transcript id, degeneracy
        write a tsv

        could probably generate sfs and diversity to two diff files in same script

        trick is probably to get to the appropriate genotype array to perform operations on
        get pipeline running up to subset vcf, then download the two files locally to work on this script
        might need to write a process to write a transcript to gene mapping function
        """
}