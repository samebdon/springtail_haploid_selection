process getGeneBedGTF {
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
        cat ${gtf_f} |  awk 'OFS="\t" {if (\$3=="gene") {print \$1,\$4-1,\$5,\$9}}' | tr -d '";' > ${species}.gene.bed
        """
}

process getGeneBedGFF {
        cpus 1
        memory = '2 GB'
        publishDir params.outdir, mode:'copy'

        input:
        path(gff_f)
        val(species)

        output:
        tuple val(species), path("${species}.gene.bed")

        script:
        """
        cat ${gff_f}| awk 'OFS="\t" {if (\$3=="gene") {print \$1,\$4-1,\$5,\$9}}' | tr -d '";' | sed 's/ID=//'> ${species}.gene.bed
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
        tuple val(species), path("degenotate/*"), emit: degen
        tuple val(species), path("${species}.cds-nt.fa"), emit: nt_fa
        tuple val(species), path("${species}.cds-nt-longest.fa"), emit: longest_isoforms

        script:
        """
        degenotate.py -g ${genome_f} -a ${annotation_f} -x 04 -o degenotate
        degenotate.py -g ${genome_f} -a ${annotation_f} -x 04 -c ${species}.cds-nt.fa -l ${species}.cds-nt-longest.fa
        """
}

process filterBed{
        cpus 1
        memory = '5 GB'
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path("degenotate/*")
        tuple val(species), path(longest_isoform_fasta)

        output:
        tuple val(species), path("${species}.longest_isoforms.0D.bed"), path("${species}.longest_isoforms.4D.bed")

        script:
        """
        python ${launchDir}/scripts/filterDegenotate.py -b <(cat degenotate/degeneracy-all-sites.bed | tr : \$'\\t' | awk -v FS='\\t' -v OFS='\\t' '{print \$1, \$2, \$3, \$4, \$6}') -t <(grep '>' ${longest_isoform_fasta}| cut -f2- -d'>') -o ${species}.longest_isoforms
        """
}

process subsetVCF{

        input:
        tuple val(species), path(zero_bed_f), path(four_bed_f)
        path(vcf_f)
        path(vcf_index)

        output:
        tuple val(species), path("${species}.0D.longest_isoforms.vcf.gz"), path("${species}.4D.longest_isoforms.vcf.gz")

        script:
        """
        bcftools view --threads ${task.cpus} -Oz -R ${zero_bed_f} -o ${species}.0D.longest_isoforms.vcf.gz ${vcf_f}
        bcftools view --threads ${task.cpus} -Oz -R ${four_bed_f} -o ${species}.4D.longest_isoforms.vcf.gz ${vcf_f}
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