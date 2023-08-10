process sortBam {

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${meta}.coord_sorted.bam")

        script:
        """
        samtools sort -m 4G -O bam -o ${meta}.coord_sorted.bam ${bam_f}
        """
}

process addRG {
	
	input:
	tuple val(meta), path(bam_f)

	output:
	tuple val(meta), path("${bam_f.baseName}.RG.bam")

	script:
	"""
	java -Xms4G -Xmx4G -jar /software/team360/picard.jar AddOrReplaceReadGroups I=${bam_f} O=${bam_f.baseName}.RG.bam RGLB=${meta} RGPL=ILLUMINA RGPU=${meta} RGSM=${meta}
	"""

}

process markDupes {

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${bam_f.baseName}.deduped.bam")

        script:
        """
        java -Xms4G -Xmx4G -jar /software/team360/picard.jar MarkDuplicates I=${bam_f} O=${bam_f.baseName}.deduped.bam M=${bam_f.baseName}.metrics.txt ASO=queryname
	"""
}

process indexBam{

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${bam_f}.bai")

        script:
        """
	samtools index ${bam_f}           
        """
}

process mosdepth {
	publishDir params.outdir, mode:'copy'

        input:
	tuple val(meta), path(bam_f), path(bam_index)

        output:
        tuple val(meta), path("mosdepth/${bam_f.baseName}.CALLABLE.bed")

        script:
        """
        mkdir mosdepth
	mosdepth -n --quantize 0:1:10:150: ${meta} ${bam_f}
	zcat ${meta}.quantized.bed.gz | grep 'CALLABLE' > mosdepth/${bam_f.baseName}.CALLABLE.bed
        """
}

process freebayes {

        input:
	path(genome_f)
        tuple val(meta), path(bam_f), path(bam_index), path(bed_f)

        output:
        tuple val(meta), path("${bam_f.baseName}.${genome_f.baseName}.vcf")        

        script:
        """
        freebayes -f ${genome_f} -b ${bam_f} -t ${bed_f} --strict-vcf -v ${bam_f.baseName}.${genome_f.baseName}.vcf -T 0.01 -p 2 -i -X -u -E 0 -n 4 -m 20 --min-coverage 10 -C 2
        """
}

process bcftools {
	publishDir params.outdir, mode:'copy'

        input:
        tuple val(meta), path(vcf_f)

        output:
        tuple val(meta), path("vcfs/${vcf_f.baseName}.sorted.filtered.vcf.gz"), path("vcfs/${vcf_f.baseName}.sorted.filtered.vcf.gz.csi")
	
        script:
        """
	mkdir vcfs
        bcftools sort ${vcf_f} -O z > ${vcf_f.baseName}.sorted.vcf.gz
	bcftools filter -O z --include "RPL >=1 && RPR>=1 & SAF>=1 && SAR>=1 && N_MISSING=0" ${vcf_f.baseName}.sorted.vcf.gz > vcfs/${vcf_f.baseName}.sorted.filtered.vcf.gz
        bcftools index vcfs/${vcf_f.baseName}.sorted.filtered.vcf.gz -o vcfs/${vcf_f.baseName}.sorted.filtered.vcf.gz.csi
	"""
}

process getHet {
	publishDir params.outdir, mode:'copy'

        input:
        tuple val(meta), path(vcf_f), path(vcf_index), path(bed)

        output:
        path("lg_het_tsvs/${meta}.lg_het.tsv")

        script:
        """
	mkdir lg_het_tsvs
        python ${launchDir}/scripts/getHet.py -v ${vcf_f} -b ${bed} -o lg_het_tsvs/${meta}
        """
}

process alleleCounter {
	publishDir params.outdir, mode:'copy'

        input:
	path(genome)
        tuple val(meta), path(vcf_f), path(vcf_index), path(bam_f), path(bam_index)

        output:
        tuple val(meta), path("alleleCounter/${meta}.allelecounter.tsv")

        script:
        """
	mkdir alleleCounter
        python ${launchDir}/scripts/allelecounter.py --vcf ${vcf_f} --sample ${meta} --bam ${bam_f} --ref ${genome} --min_cov 10 --min_baseq 1 --min_mapq 1 --out allelecounter/${meta}.allelecounter.tsv
        """
}

process gatk_aseReadCount {
	publishDir params.outdir, mode:'copy'

	input:
        path(genome)
	path(genome_index)
	path(genome_dict)
        tuple val(meta), path(vcf_f), path(vcf_index), path(bam_f), path(bam_index)

	output:
        tuple val(meta), path("aseReadCount/${meta}.allelecounter.tsv")

	script:
	"""
	mkdir aseReadCount
	gatk ASEReadCounter -R ${genome} -I ${bam_f} -V {vcf_f} -O aseReadCount/${meta}.aseReadCount.tsv
	"""
}
