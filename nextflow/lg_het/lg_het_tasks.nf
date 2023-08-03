process markDupes {

        input:
        path bam_f

        output:
        path ("${bam_f.baseName}.marked.bam")

        script:
        """
        java -Xms4G -Xmx4G -jar /software/team360/picard.jar MarkDuplicates I=${bam_f} O=${bam_f.baseName}.marked.bam M=${bam_f.baseName}.metrics.txt ASO=queryname
	"""
}

process sortBam{

        input:
        path bam_f

        output:
        path("${bam_f.baseName}.sorted.bam")

        script:
        """
        samtools sort ${bam_f} > ${bam_f.baseName}.sorted.bam
        """
}

process indexBam{

        input:
        path bam_f

        output:
        path("${bam_f}.bai"), optional:true, emit: bai

        script:
        """
	samtools index ${bam_f}           
        """
}

process mosdepth {
	publishDir params.outdir, mode:'copy'

        input:
        path bam_f
	path bam_index

        output:
        path("mosdepth/${bam_f.baseName}.CALLABLE.bed")

        script:
        """
        mkdir mosdepth
	mosdepth -n --quantize 0:1:10:150: ${bam_f.baseName} ${bam_f}
	zcat ${bam_f.baseName}.quantized.bed.gz | grep 'CALLABLE' > mosdepth/${bam_f.baseName}.CALLABLE.bed
        """
}

process freebayes {

        input:
        path genome_f
	path bam_f
	path bam_index
	path bed_f

        output:
        path("${bam_f.baseName}.${genome_f.baseName}.vcf")        

        script:
        """
        freebayes -f ${genome_f} -b ${bam_f} -t ${bed_f} --strict-vcf -v ${bam_f.baseName}.${genome_f.baseName}.vcf -T 0.01 -p 2 -i -X -u -E 0 -n 4 -m 20 --min-coverage 10 -C 2
        """
}

process bcftools {
	publishDir params.outdir, mode:'copy'

        input:
        path vcf_f

        output:
        path("vcfs/${vcf_f.baseName}.sorted.filtered.vcf.gz")
	path("vcfs/vcfs/${vcf_f.baseName}.sorted.filtered.vcf.gz.csi")
	
        script:
        """
	mkdir vcfs
        bcftools sort ${vcf_f} -O z > ${vcf_f.baseName}.sorted.vcf.gz
	bcftools filter -O z --include "RPL >=1 && RPR>=1 & SAF>=1 && SAR>=1 && N_MISSING=0" ${vcf_f.baseName}.sorted.vcf.gz > vcfs/${vcf_f.baseName}.sorted.filtered.vcf.gz
        bcftools index vcfs/${vcf_f.baseName}.sorted.filtered.vcf.gz
	"""
}

process getHet {
	publishDir params.outdir, mode:'move'

        input:
        path vcf_f

        output:
        path("lg_het_tsvs/${vcf_f.baseName}.lg_het.tsv")

        script:
        """
	mkdir lg_het_tsvs
        python getHet.py ${vcf_f} ${lg_het_tsvs/${vcf_f.baseName}.lg_het.tsv}
        """
}
