process bwaIndex {

        input:
        path(genome_f)

        output:
        path("${genome_f}.*")

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
        path("mosdepth/${bam_f.baseName}.CALLABLE.bed")

        script:
        """
        mkdir mosdepth
        mosdepth --fast-mode -t 4 tmp ${bam_f}
        MAX_DEPTH="\$(python scripts/mean_depth.py -b tmp.per-base.bed.gz -m 2)" 
	mosdepth -t 4 -n --quantize 0:1:8:\${MAX_DEPTH}: ${meta} ${bam_f}
	zcat ${meta}.quantized.bed.gz | grep 'CALLABLE' > mosdepth/${bam_f.baseName}.callable.bed
        """
}

process intersectBeds{

        input:
        path(beds, stageAs: "inputs/*")
        val(species)

        output:
        tuple val(species), path("${species}.intersect.bed")

        script:
        """
        bedtools multiinter -i $beds | cut -f1-5 > ${species}.intersect.bed
        """
}

process samtoolsMerge {

        input:
        tuple val(meta), path('*.bam'), path(bam_index)
        val(species)

        output:
        tuple val(species), path('${species}.bam'), path("${species}.bam.bai")

        script:
        """
        samtools merge -o ${species}.bam *.bam
        samtools index ${bam_f}           
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

process bcftools_filter {

        input:
        tuple val(species), path(vcf_f)

        output:
        tuple val(species), path("${vcf_f.baseName}.soft_filtered.vcf.gz")
	
        script:
        """
        bcftools filter --threads 4 -Oz -s Qual -m+ -e 'QUAL<1' ${vcf_f} | \
        bcftools filter --threads 4 -Oz -s Balance -m+ -e 'RPL<1 | RPR<1 | SAF<1 | SAR<1' | \
        bcftools filter --threads 4 -Oz -m+ -s+ --SnpGap 2 | \
        bcftools filter --threads 4 -Oz -e 'TYPE!="snp"' -s NonSnp -m+ > ${vcf_f.baseName}.soft_filtered.vcf.gz
        """
}

process generate_fail_bed {
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path(vcf_f)

        output:
        tuple val(species), path("${species}.vcf_filter_fails.bed")

        script:
        """
        bcftools view --threads 4 -H -i "%FILTER!='PASS'" ${vcf_f} | \
        perl -lane '$pad=0; print($F[0]."\t".($F[1]-1)."\t".(($F[1]-1)+length($F[3]))."\t".$F[6])' | \
        bedtools sort | \
        bedtools merge > ${species}.vcf_filter_fails.bed
        """
}

process generate_pass_vcf {

        input:
        tuple val(species), path(vcf_f)

        output:
        tuple val(species), path("${vcf_f.baseName}.sorted.hard_filtered.vcf.gz")

        script:
        """
        bcftools view --threads 4 -Oz -f "PASS" ${vcf_f} > ${vcf_f.baseName}.sorted.hard_filtered.vcf.gz
        """
}

process bedtools_subtract {

        input:
        tuple val(species), path(a.bed)
        tuple val(species), path(b.bed)

        output:
        tuple val(species), path("{species}.callable.bed")

        script:
        """
        bedtools subtract -a ${a.bed} -b ${b.bed} > ${species}.callable.bed
        """
}

process bcftools_sort {
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path(vcf_f)

        output:
        tuple val(species), path("${species}.hard_filtered.sorted.vcf.gz")
        
        script:
        """
        bcftools sort --threads 4 -Oz ${vcf_f} > ${species}.hard_filtered.sorted.vcf.gz
        """
}        

process bcftools_index {
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(meta), path(vcf_f)

        output:
        tuple val(meta), path("${vcf_f.baseName}.csi")

        script:
        """
        bcftools index -t ${vcf_f} -o ${vcf_f.baseName}.csi
        """
}
