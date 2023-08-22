process bwaIndex {

        input:
        path(genome_f)

        output:
        path('${genome_f}.*')

        script:
        """
        bwa index ${genome_f}
        """
}

process bwaMem {
        publishDir params.outdir, mode:'copy'

        input:
        path(genome_f)
        path(genome_index)
        tuple val(meta), path(reads)

        output:
        tuple val(meta), path("bwamem/${meta}.${genome.baseName}.bam")

        script:
        """
        mkdir bwamem
        bwa mem -t 4 -R "@RG\tID:${meta}\tSM:${meta}\tPL:ILLUMINA\tPU:${meta}\tLB:${meta}\tDS:${meta}" ${genome_f} ${reads[0]} ${reads[1]} > bwamem/${meta}.${genome.baseName}.bam
        """
}

process sortBam {

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${bam_f.baseName}.coord_sorted.bam")

        script:
        """
        samtools sort -@ 4 -m 4G -O bam -o ${bam_f.baseName}.coord_sorted.bam ${bam_f}
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

        input:
	tuple val(meta), path(bam_f), path(bam_index)

        output:
        tuple val(meta), path("mosdepth/${bam_f.baseName}.CALLABLE.bed")

        script:
        """
        mkdir mosdepth
        mosdepth --fast-mode -t 4 tmp ${bam_f}
        MAX_DEPTH="$(python scripts/mean_depth.py -b tmp.per-base.bed.gz -m 2)" 
	mosdepth -t 4 -n --quantize 0:1:8:\${MAX_DEPTH}: ${meta} ${bam_f}
	zcat ${meta}.quantized.bed.gz | grep 'CALLABLE' > mosdepth/${bam_f.baseName}.callable.bed
        """
}

process bedtoolsIntersect{
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(meta), path(a_bed)
        tuple val(meta), path('*.bed')
        val(species)

        output:
        tuple val(species), path("${species}.callable.bed")

        script:
        """
        bedtools intersect -a ${a_bed} -b *.bed > ${species}.callable.bed 
        """
}

process samtoolsMerge {

        input:
        tuple val(meta), path('*.bam'), path(bam_index)
        val(species)

        output:
        tuple val(species), path('${species}.bam')

        script:
        """
        samtools merge -o ${species}.bam *.bam
        """
}

process freebayes {

        input:
	path(genome_f)
        tuple val(meta), path(bam_f), path(bam_index)
        tuple val(species), path(bed_f)

        output:
        tuple val(species), path("${species}.vcf")        

        script:
        """
        freebayes -f ${genome_f} -b ${bam_f} -t ${bed_f} --strict-vcf -v ${species}.vcf -T 0.01 -p 2 -i -X -u -E 0 -n 4 -m 20 --min-coverage 10 --min-mapping-quality 1 -C 2
        """
}

process bcftools {
	publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path(vcf_f)

        output:
        tuple val(species), path("${vcf_f.baseName}.sorted.filtered.vcf.gz"), path("vcfs/${vcf_f.baseName}.sorted.filtered.vcf.gz.csi")
	
        script:
        """
        bcftools sort --threads 4 -Oz ${vcf_f} | \
        bcftools filter --threads 4 -Oz -s Qual -m+ -e 'QUAL<1' | \
        bcftools filter --threads 4 -Oz -s Balance -m+ -e 'RPL<1 | RPR<1 | SAF<1 | SAR<1' | \
        bcftools filter --threads 4 -Oz -m+ -s+ --SnpGap 2 | \
        bcftools filter --threads 4 -Oz -e 'TYPE!="snp"' -s NonSnp -m+ > ${vcf_f.baseName}.sorted.filtered.vcf.gz
        bcftools index ${vcf_f.baseName}.sorted.filtered.vcf.gz -o ${vcf_f.baseName}.sorted.filtered.vcf.gz.csi
        """
}