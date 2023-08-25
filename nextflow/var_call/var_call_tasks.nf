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
        tuple val(meta), path("bwamem/${meta}.${genome_f.baseName}.bam")

        script:
        """
        mkdir bwamem
        bwa mem -t ${task.cpus} -R "@RG\\tID:${meta}\\tSM:${meta}\\tPL:ILLUMINA\\tPU:${meta}\\tLB:${meta}\\tDS:${meta}" ${genome_f} ${reads[0]} ${reads[1]} | \
        sambamba view -t ${task.cpus} -S -f bam /dev/stdin > bwamem/${meta}.${genome_f.baseName}.bam
        """
}

process sortBam {

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${bam_f.baseName}.coord_sorted.bam")

        script:
        """
        samtools sort -@ ${task.cpus} -O bam -o ${bam_f.baseName}.coord_sorted.bam ${bam_f}
        """
}

process sortBamSambamba {

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${bam_f.baseName}.coord_sorted.bam")

        script:
        avail_mem = (task.memory.mega*1).intValue()
        """
        sambamba sort -t ${task.cpus} -m ${avail_mem}MB -o ${bam_f.baseName}.coord_sorted.bam ${bam_f}
        """
}

process markDupes {

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${bam_f.baseName}.deduped.bam")

        script:
        avail_mem = (task.memory.mega*0.8).intValue()
        """
        java -Xmx${avail_mem}M -jar /software/team360/picard.jar MarkDuplicates I=${bam_f} O=${bam_f.baseName}.deduped.bam M=${bam_f.baseName}.metrics.txt ASO=queryname
        """
}

process markDupesSambamba {

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${bam_f.baseName}.deduped.bam")

        script:
        """
        sambamba markdup -t ${task.cpus} ${bam_f} ${bam_f.baseName}.deduped.bam
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

process indexBamSambamba{

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${bam_f}.bai")

        script:
        """
        sambamba index -t ${task.cpus} ${bam_f}             
        """
}

process mosdepth {

        input:
	tuple val(meta), path(bam_f), path(bam_index)
        val(min_depth)
        val(max_depth_factor)

        output:
        path("mosdepth/${bam_f.baseName}.callable.bed")

        script:
        """
        mkdir mosdepth
        mosdepth --fast-mode -t ${task.cpus} tmp ${bam_f}
        MAX_DEPTH="\$(python ${launchDir}/scripts/mean_depth.py -b tmp.per-base.bed.gz -m ${max_depth_factor})" 
	mosdepth -t ${task.cpus} -n --quantize 0:1:${min_depth}:\${MAX_DEPTH}: ${meta} ${bam_f}
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
        samtools index ${species}.bam        
        """
}

process sambambaMerge {

        input:
        path(bams)
        val(species)

        output:
        tuple val(species), path('${species}.bam'), path("${species}.bam.bai")

        script:
        """
        sambamba merge -t ${task.cpus} ${species}.bam ${bams.join(" ")}
        sambamba index -t ${task.cpus} ${species}.bam         
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
        freebayes -f ${genome_f} -b ${bam_f} -t ${bed_f} --strict-vcf -v ${species}.vcf -T 0.01 -k -w -j -E 1
        """
}

process bcftools_filter {

        input:
        path(genome)
        tuple val(species), path(vcf_f)

        output:
        tuple val(species), path("${vcf_f.baseName}.soft_filtered.vcf.gz")
	
        script:
        """
        bcftools norm --threads ${task.cpus} -Ov -f ${genome} ${vcf_f} | \
        vcfallelicprimitives --keep-info --keep-geno -t decomposed | \
        bcftools plugin fill-AN-AC --threads ${task.cpus} -Oz | \
        bcftools filter --threads ${task.cpus} -Oz -s Qual -m+ -e 'QUAL<10' | \
        bcftools filter --threads ${task.cpus} -Oz -s Balance -m+ -e 'RPL<1 | RPR<1 | SAF<1 | SAR<1' | \
        bcftools filter --threads ${task.cpus} -Oz -m+ -s+ --SnpGap 2 | \
        bcftools filter --threads ${task.cpus} -Oz -e 'TYPE!="snp"' -s NonSnp -m+ > ${vcf_f.baseName}.soft_filtered.vcf.gz
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
        bcftools view --threads ${task.cpus} -H -i "%FILTER!='PASS'" ${vcf_f} | \
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
        bcftools view --threads ${task.cpus} -Oz -f "PASS" ${vcf_f} > ${vcf_f.baseName}.sorted.hard_filtered.vcf.gz
        """
}

process bedtools_subtract {

        input:
        tuple val(species), path(a_bed)
        tuple val(species), path(b_bed)

        output:
        tuple val(species), path("{species}.callable.bed")

        script:
        """
        bedtools subtract -a ${a_bed} -b ${b_bed} > ${species}.callable.bed
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
        bcftools sort --threads ${task.cpus} -Oz ${vcf_f} > ${species}.hard_filtered.sorted.vcf.gz
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
