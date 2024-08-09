process getLongestIsoformAGAT{
        memory '4G'

        input:
        path(annotation)

        output:
        path("${annotation.simpleName}.agat.longest_isoform.gff3")

        script:
        """
        agat_sp_keep_longest_isoform.pl -gff ${annotation} -o ${annotation.simpleName}.agat.longest_isoform.gff3
        """
}

process makeGenomeFile{
        publishDir params.outdir, mode:'copy'
        conda '/software/treeoflife/conda/users/envs/team360/se13/gene_pop'

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

// cat ${annotation} |  awk 'OFS="\t" {if (\$3=="gene") {print \$1,\$4-1,\$5,\$9}}' | tr -d '";' > 

process getGeneBedAGAT {
        memory '4G'
        publishDir params.outdir, mode:'copy'

        input:
        path(annotation)
        val(species)

        output:
        tuple val(species), path("${species}.gene.AGAT.bed")

        script:
        """
        agat_convert_sp_gff2bed.pl -gff ${annotation} -o ${species}.gene.bed --sub gene -o ${species}.gene.bed
        """
}

process splitBed {
        conda '/software/treeoflife/conda/users/envs/team360/se13/gene_pop'

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
        conda '/software/treeoflife/conda/users/envs/team360/se13/gene_pop'

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
        conda '/software/treeoflife/conda/users/envs/team360/se13/gene_pop'

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

// if the zero_bed_f and four_bed_f are per chromosome, how is the gene_bed_f and zero_bed_f lining up?
// is getting gene names from zero/four bed files enough?
process calculatePiBed{
        conda '/software/treeoflife/conda/users/envs/team360/se13/gene_pop'

        input:
        path(vcf_f)
        path(vcf_index)
        tuple val(species), path(zero_bed_f), path(four_bed_f)
        path(gene_bed_f)
        tuple val(species), path(genome_file)

        output:
        tuple val(species), val(gene_bed_f.simpleName), path("${species}.${gene_bed_f.simpleName}.0D.longest_isoforms.pi.tsv"), path("${species}.${gene_bed_f.simpleName}.4D.longest_isoforms.pi.tsv"), emit: pi
        path("*.sfs.txt"), emit: sfs

        script:
        """
        calculatePiBed.py -v ${vcf_f} -b ${gene_bed_f} -a <(grep ${gene_bed_f.simpleName} ${zero_bed_f} ) -n ${gene_bed_f.simpleName}.1 -o ${species}.${gene_bed_f.simpleName}.0D.longest_isoforms.pi.tsv -l 0D_pi -g ${genome_file}
        calculatePiBed.py -v ${vcf_f} -b ${gene_bed_f} -a <(grep ${gene_bed_f.simpleName} ${four_bed_f} ) -n ${gene_bed_f.simpleName}.1 -o ${species}.${gene_bed_f.simpleName}.4D.longest_isoforms.pi.tsv -l 4D_pi -g ${genome_file}
        """
}

process mergePi{
        conda '/software/treeoflife/conda/users/envs/team360/se13/gene_pop'

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
        conda '/software/treeoflife/conda/users/envs/team360/se13/gene_pop'

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

process concat_SFS{
        publishDir params.outdir, mode:'copy'
        input:
        path(files, stageAs: "inputs/*")
        val(species)

        output:
        path('*')

        script:
        """
        agg_sfs.py -i inputs -x OX359249,OX359250 -o ${species}
        """
}