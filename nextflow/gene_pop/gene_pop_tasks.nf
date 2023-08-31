process getGeneBed {
        cpus 1
        memory = '2 GB'
        publishDir params.outdir, mode:'copy'

        input:
        path(gtf_f)
        val(species)

        output:
        tuple val(species), path("${species}.gene.bed")

        script:
        """
        cat ${gtf_f} |  awk 'OFS="\t" {if (\$3=="gene") {print \$1,\$4-1,\$5,\$10}}' | tr -d '";'> ${species}.gene.bed
        """
}


process splitBed {
        cpus 1
        memory = '2 GB'

        input:
        tuple val(species), path(bed_f)

        output:
        path("beds/*")

        script:
        """
        bash ${launchDir}/scripts/splitBed.sh ${bed_f}
        """
}

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
        cat ${annotation_f} | awk -v OFS='\\t' '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' > tmp.gtf
        degenotate.py -g ${genome_f} -a tmp.gtf -o degenotate -x 04 -c ${species}.cds-nt.fa -ca ${species}.cds-aa.fa -l ${species}.cds-nt-longest.fa -la ${species}.cds-aa-longest.fa
        """
}

process filterBed{
        cpus 1
        memory = '1 GB'
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path("degenotate/*")

        output:
        tuple val(species), path("${species}.0D.longest_isoforms.bed"), path("${species}.4D.longest_isoforms.bed")

        script:
        """
        join degenotate/degeneracy-all-sites.bed <(grep '>' degenotate/${species}.cds-nt-longest.fa | cut -f1- -d'>') -1 4 -2 1 > ${species}.04D.longest_isoforms.bed
        split ${species}.04D.longest_isoforms.bed into 0D and 4D
        """
}

process subsetVCF{

        input:
        tuple val(species), path(zero_bed_f), path(four_bed_f)
        path(vcf_f)

        output:
        tuple val(species), path("${species}.0D.longest_isoforms.vcf.gz"), path("${species}.0D.longest_isoforms.vcf.gz")

        script:
        """
        bcftools view --threads ${task.cpus} -Oz -R <(cat ${zero_bed_f} | cut -f1-3) -o ${species}.0D.longest_isoforms.vcf.gz ${vcf_f}
        bcftools view --threads ${task.cpus} -Oz -R <(cat ${four_bed_f} | cut -f1-3) -o ${species}.4D.longest_isoforms.vcf.gz ${vcf_f}
        """
}

process calculatePiBed{
        cpus 1
        memory = '2 GB'

        input:
        tuple val(species), path(zero_vcf_f), path(four_vcf_f)
        tuple val(species), path(zero_bed_f), path(four_bed_f)
        path(gene_bed_f)

        output:
        tuple val(species), val(gene_bed_f.simpleName), path("${species}.0D.longest_isoforms.theta.tsv"), path("${species}.4D.longest_isoforms.theta.tsv")

        script:
        """
        python ${launchDir}/scripts/calculatePiBed.py -v ${zero_vcf_f} -b ${gene_bed_f} -a <(grep ${gene_bed_f.simpleName} ${zero_bed_f} | cut -f1-3) -n ${gene_bed_f.simpleName} -o ${species}.${gene_bed_f.simpleName}.0D.longest_isoforms.pi.tsv -l 0D_pi
        python ${launchDir}/scripts/calculatePiBed.py -v ${four_vcf_f} -b ${gene_bed_f} -a <(grep ${gene_bed_f.simpleName} ${four_bed_f} | cut -f1-3) -n ${gene_bed_f.simpleName} -o ${species}.${gene_bed_f.simpleName}.4D.longest_isoforms.pi.tsv -l 4D_pi
        """
}

process joinPi{
        cpus 1
        memory = '2 GB'

        input:
        tuple val(species), val(chrom), path(zero_f), path(four_f)

        output:
        tuple val(species), path("${species}.${chrom}.longest_isoforms.pi.tsv")

        script:
        """
        join -1 2 -2 1 ${zero_f} <(cat ${four_F} | cut -f2-3) > ${species}.${chrom}.longest_isoforms.pi.tsv
        """
}

process concat_all{
        cpus 1
        memory = '2 GB'
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path(files, stageAs: "inputs/*")

        output:
        tuple val(species), path(output_f)
        
        script:
        """
        awk '(NR == 1) || (FNR > 1)' inputs/* > ${species}.longest_isoforms.pi.tsv
        """
}