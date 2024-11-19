process trimReads {
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(sample_id), path(reads)

        output:
        tuple val(sample_id), path('fastp/*.fastp.fastq.gz') , optional:true, emit: reads

        script:
        """
        mkdir fastp
        fastp -i ${reads[0]} -I ${reads[1]} -o fastp/${sample_id}.1.fastp.fastq.gz -O fastp/${sample_id}.2.fastp.fastq.gz --length_required 33 --cut_front --cut_tail --cut_mean_quality 20 --thread ${task.cpus}
        """
}

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
        queue 'long'
        memory '8G'

        input:
        path(genome_f)
        path(genome_index)
        tuple val(meta), path(reads)

        output:
        tuple val(meta), path("bwamem/${meta}.${genome_f.baseName}.bam")

        script:
        """
        mkdir -p bwamem
        bwa mem -t ${task.cpus} -R "@RG\\tID:${meta}\\tSM:${meta}\\tPL:ILLUMINA\\tPU:${meta}\\tLB:${meta}\\tDS:${meta}" ${genome_f} ${reads[0]} ${reads[1]} | \
        sambamba view -t ${task.cpus} -S -f bam /dev/stdin > bwamem/${meta}.${genome_f.baseName}.bam
        """
}

process sortBamSambamba {
        publishDir params.outdir, mode:'copy'
        memory '8G'

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

process markDupesSambamba {
        memory '4G'

        input:
        tuple val(meta), path(bam_f)

        output:
        tuple val(meta), path("${bam_f.baseName}.deduped.bam"), emit: meta_bam
        path("${bam_f.baseName}.deduped.bam"), emit: bam_only

        script:
        """
        sambamba markdup -t ${task.cpus} ${bam_f} ${bam_f.baseName}.deduped.bam
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

// mosdepth could do with having a fail retry clause
process mosdepth {
        memory '15G'

        input:
	tuple val(meta), path(bam_f), path(bam_index)
        val(min_depth)

        output:
        path("mosdepth/${bam_f.baseName}.callable.bed")

        script:
        """
        mkdir -p mosdepth
        mosdepth --fast-mode -t ${task.cpus} tmp ${bam_f}
        MAX_DEPTH="\$(max_depth.py -b tmp.per-base.bed.gz)" 
	mosdepth -t ${task.cpus} -n --quantize 0:1:${min_depth}:\${MAX_DEPTH}: ${meta} ${bam_f}
	zcat ${meta}.quantized.bed.gz | grep 'CALLABLE' > mosdepth/${bam_f.baseName}.callable.bed
        """
}

process intersectBeds{
        publishDir params.outdir, mode:'copy'
        memory '4G'

        input:
        path(beds, stageAs: "inputs/*")
        path(repeat_bed)
        path(genome_index)
        val(species)

        output:
        tuple val(species), path("${species}.callable.all.bed"), emit: all
        tuple val(species), path("${species}.callable.freebayes.norepeats.bed"), emit: freebayes

        script:
        """
        N_FILES="\$(ls inputs/*.bed | wc -l)"
        bedtools multiinter -i $beds | cut -f1-5 | bedtools sort -faidx ${genome_index} > ${species}.callable.all.bed
        cat ${species}.callable.all.bed | awk -v var=\$N_FILES '\$4==var'  | cut -f1-3 > ${species}.callable.freebayes.bed
        bedtools subtract -a ${species}.callable.freebayes.bed -b ${repeat_bed} | bedtools sort -faidx ${genome_index} | bedtools merge > ${species}.callable.freebayes.norepeats.bed
        """
}

process intersectBed{
        publishDir params.outdir, mode:'copy'
        memory '4G'

        input:
        path(bed)
        path(repeat_bed)
        path(genome_index)
        val(species)

        output:
        tuple val(species), path("${species}.callable.freebayes.norepeats.bed")

        script:
        """
        bedtools subtract -a ${bed} -b ${repeat_bed} | bedtools sort -faidx ${genome_index} | bedtools merge > ${species}.callable.freebayes.norepeats.bed
        """
}

// could make a separate process to remove repeats which is optional for if you have a repeat annotation

process sambambaMerge {
        memory '8G'

        input:
        path(bams)
        val(species)

        output:
        tuple val(species), path("${species}.bam"), path("${species}.bam.bai")

        script:
        """
        sambamba merge -t ${task.cpus} ${species}.bam ${bams.join(" ")}
        sambamba index -t ${task.cpus} ${species}.bam         
        """
}

process freebayes {
        queue 'basement'
        memory '20G'
        cpus 1

        input:
	path(genome_f)
        path(genome_index)
        tuple val(meta), path(bam_f), path(bam_index)
        tuple val(species), path(bed_f)

        output:
        tuple val(species), path("${species}.vcf")        

        script:
        """
        freebayes -f ${genome_f} -b ${bam_f} -t ${bed_f} --strict-vcf -v ${species}.vcf -T 0.01 -k -w -j -E 1
        """
}

process freebayesParallel {
        queue 'long'
        memory '30G'

        input:
        path(genome_f)
        path(genome_index)
        tuple val(meta), path(bam_f), path(bam_index)
        tuple val(species), path(bed_f)

        output:
        tuple val(species), path("${species}.vcf")        

        script:
        """
        freebayes-parallel <(fasta_generate_regions.py ${genome_f} 10000000) ${task.cpus}  -f ${genome_f} -b ${bam_f} -t ${bed_f} --strict-vcf -T 0.01 -k -w -j -E 1 > ${species}.vcf
        """
}

// vcftools=1.0.3
process bcftools_filter {
        publishDir params.outdir, mode:'copy'

        input:
        path(genome)
        tuple val(species), path(vcf_f)

        output:
        tuple val(species), path("${vcf_f.baseName}.soft_filtered.vcf.gz")
	
        script:
        """
        bcftools norm --threads ${task.cpus} -Ov -f ${genome} ${vcf_f} | \
        vcfallelicprimitives --keep-info --keep-geno -t decomposed | \
        bcftools +fill-tags --threads ${task.cpus} -Oz -- -t AN,AC,F_MISSING | \
        bcftools filter --threads ${task.cpus} -Oz -s Qual -m+ -e 'QUAL<10' | \
        bcftools filter --threads ${task.cpus} -Oz -s Balance -m+ -e 'RPL<1 | RPR<1 | SAF<1 | SAR<1' | \
        bcftools filter --threads ${task.cpus} -Oz -m+ -s+ --SnpGap 2 | \
        bcftools filter --threads ${task.cpus} -Oz -e 'TYPE!="snp"' -s NonSnp -m+ > ${vcf_f.baseName}.soft_filtered.vcf.gz
        """
}

process generate_fail_bed {
        publishDir params.outdir, mode:'copy'
        memory '4G'

        input:
        tuple val(species), path(vcf_f)
        path(genome_index)

        output:
        tuple val(species), path("${species}.vcf_filter_fails.bed")

        script:
        """
        bcftools view --threads ${task.cpus} -H -i "FILTER!='PASS'" ${vcf_f} | \
        perl -lane '\$pad=0; print(\$F[0]."\\t".(\$F[1]-1)."\\t".((\$F[1]-1)+length(\$F[3]))."\\t".\$F[6])' | \
        bedtools sort -faidx ${genome_index} | \
        bedtools merge > ${species}.vcf_filter_fails.bed
        """
}

process generate_pass_vcf {

        input:
        tuple val(species), path(vcf_f)

        output:
        tuple val(species), path("${vcf_f.baseName}.hard_filtered.vcf.gz")

        script:
        """
        bcftools view --threads ${task.cpus} -Oz -f "PASS" ${vcf_f} > ${vcf_f.baseName}.hard_filtered.vcf.gz
        """
}

process bedtools_subtract {
        cpus 1
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path(a_bed)
        tuple val(species), path(b_bed)
        path(genome_index)

        output:
        tuple val(species), path("${species}.callable.bed")

        script:
        """
        bedtools subtract -a ${a_bed} -b ${b_bed} | bedtools sort -faidx ${genome_index} | bedtools merge > ${species}.callable.bed
        """
}

process bcftools_sort {
        cpus 1
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(species), path(vcf_f)

        output:
        tuple val(species), path("${species}.hard_filtered.sorted.vcf.gz")
        
        script:
        """
        bcftools sort -Oz ${vcf_f} > ${species}.hard_filtered.sorted.vcf.gz
        """
}        

process bcftools_index {
        cpus 1
        publishDir params.outdir, mode:'copy'

        input:
        tuple val(meta), path(vcf_f)

        output:
        tuple val(meta), path("${vcf_f}.csi")

        script:
        """
        bcftools index -c ${vcf_f} -o ${vcf_f}.csi
        """
}
