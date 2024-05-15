process makeGenomeFile{
        publishDir params.outdir, mode:'copy'

        input:
        path(genome_dict)
        val(species)

        output:
        tuple val(species), path("${species}.genomefile")

        script:
        """
        cat ${genome_dict} | grep '^@SQ' | cut -f2-3 | perl -lpe 's/SN:|LN://g' > ${species}.genomefile
        """
}

process getGeneBedGTF {
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

process getExonBedGFF {
        publishDir params.outdir, mode:'copy'

        input:
        path(gff_f)
        val(species)

        output:
        tuple val(species), path("${species}.exon.bed")

        script:
        """
        cat ${gff_f}| awk 'OFS="\t" {if (\$3=="exon") {print \$1,\$4-1,\$5,\$9}}' | tr -d '";' | sed 's/ID=//'> ${species}.exon.bed
        """
}

process splitBed {

        input:
        tuple val(species), path(bed_f)

        output:
        path("beds/OX*[567890].1.bed")

        script:
        """
        bash splitBed.sh ${bed_f}
        """
}

process degenotate {
        memory '8G'
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
        memory '8G'
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path("degenotate/*")
        tuple val(species), path(longest_isoform_fasta)

        output:
        tuple val(species), path("${species}.longest_isoforms.0D.bed"), path("${species}.longest_isoforms.4D.bed")

        script:
        """
        filterDegenotate.py -b <(cat degenotate/degeneracy-all-sites.bed | tr : \$'\\t' | awk -v FS='\\t' -v OFS='\\t' '{print \$1, \$2, \$3, \$4, \$6}') -t <(grep '>' ${longest_isoform_fasta}| cut -f2- -d'>') -o ${species}.longest_isoforms
        """
}

process subsetVCF{

        input:
        tuple val(species), path(zero_bed_f), path(four_bed_f)
        path(vcf_f)
        path(vcf_index)

        output:
        tuple val(species), path("${species}.0D.longest_isoforms.vcf.gz"), path("${species}.4D.longest_isoforms.vcf.gz"), path("${species}.0D.longest_isoforms.vcf.gz.tbi"), path("${species}.4D.longest_isoforms.vcf.gz.tbi")

        script:
        """
        bcftools view --threads ${task.cpus} -Ov -R ${zero_bed_f} -o ${species}.0D.longest_isoforms.vcf ${vcf_f}
        bcftools view --threads ${task.cpus} -Ov -R ${four_bed_f} -o ${species}.4D.longest_isoforms.vcf ${vcf_f}
        bgzip ${species}.0D.longest_isoforms.vcf && tabix ${species}.0D.longest_isoforms.vcf.gz
        bgzip ${species}.4D.longest_isoforms.vcf && tabix ${species}.4D.longest_isoforms.vcf.gz
        """
}

process calculatePiBed{

        input:
        path(vcf_f)
        path(vcf_index)
        tuple val(species), path(zero_bed_f), path(four_bed_f)
        path(gene_bed_f)
        tuple val(species), path(genome_file)

        output:
        tuple val(species), val(gene_bed_f.simpleName), path("${species}.${gene_bed_f.simpleName}.0D.longest_isoforms.pi.tsv"), path("${species}.${gene_bed_f.simpleName}.4D.longest_isoforms.pi.tsv")

        script:
        """
        calculatePiBed.py -v ${vcf_f} -b ${gene_bed_f} -a <(grep ${gene_bed_f.simpleName} ${zero_bed_f} ) -n ${gene_bed_f.simpleName}.1 -o ${species}.${gene_bed_f.simpleName}.0D.longest_isoforms.pi.tsv -l 0D_pi -g ${genome_file}
        calculatePiBed.py -v ${vcf_f} -b ${gene_bed_f} -a <(grep ${gene_bed_f.simpleName} ${four_bed_f} ) -n ${gene_bed_f.simpleName}.1 -o ${species}.${gene_bed_f.simpleName}.4D.longest_isoforms.pi.tsv -l 4D_pi -g ${genome_file}
        """
}

process mergePi{

        input:
        tuple val(species), val(chrom), path(zero_f), path(four_f)

        output:
        path("${species}.${chrom}.longest_isoforms.pi.tsv")

        script:
        """
        mergePi.py -z ${zero_f} -f ${four_f} -o ${species}.${chrom}.longest_isoforms.pi.tsv
        """
}

process concat_all{
        publishDir params.outdir, mode:'copy'

        input:
        path(files, stageAs: "inputs/*")
        val(species)

        output:
        tuple val(species), path("${species}.longest_isoforms.pi.tsv")
        
        script:
        """
        awk '(NR == 1) || (FNR > 1)' inputs/* > ${species}.longest_isoforms.pi.tsv
        """
}
