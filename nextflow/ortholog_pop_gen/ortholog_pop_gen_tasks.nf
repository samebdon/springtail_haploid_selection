// Output of the pairwise comparisons could be a matrix of sample comparisons at each gene at 0d and 4d sites, where the diagonal would be diversity


process get_best_cds_bed {

        input:
        val(meta)
        val(annotation)

        output:
        tuple val(meta), path("${meta}.longest_isoform.gff3"), emit: gff
        tuple val(meta), path("${meta}.cds.bed"), emit: bed

        // The bed file at the end of this produces one line per transcript with all the
        // exon start positions relative to the transcript start in the last column
        script:
        """
        agat_sp_keep_longest_isoform.pl -gff ${annotation} -o ${meta}.longest_isoform.gff3
        awk '\$3=="exon"' ${meta}.longest_isoform.gff3 > ${meta}.longest_isoform.cds.gff3
        agat_convert_sp_gff2bed.pl --gff ${meta}.longest_isoform.cds.gff3 -o ${meta}.cds.bed
        """
}

process get_best_pep_fasta {

        input:
        tuple val(meta), path(cds_bed)
        path(pep_fasta)

        output:
        tuple val(meta), path("${meta}.best.pep.fasta")

        // Should be fine, check best.pep.fasta looks ok
        script:
        """
        cut -f4 ${cds_bed} > prots.lst
        fastaqual_select.pl -f ${pep_fasta} | cut -f1 -d" " > selected.pep
        grep --no-group-separator -A1 -wFf prots.lst selected.pep | sed "s/>/&${meta}./g" > ${meta}.best.pep.fasta
        """
}

process get_callable_cds_bed {

        input:
        tuple val(meta), path(cds_bed)
        path(callable_bed)

        output:
        tuple val(meta), path("${meta}.callable.cds.bed")

        // -split selects exons from the bed12 only
        script:
        """
        bedtools intersect -split -a ${cds_bed} -b ${callable_bed} > ${meta}.callable.cds.bed
        """
}

process make_genome_file {

        input:
        path(genome_fasta)

        output:
        path("${genome_fasta.baseName}.genomefile")

        script:
        """
        samtools faidx ${genome_fasta} && cut -f1,2 ${genome_fasta}.fai | sort -Vk1 > ${genome_fasta.baseName}.genomefile
        """
}

process get_mask_bed {

        input:
        tuple val(meta), path(callable_bed)
        path(genome_file)

        output:
        tuple val(meta), path("${meta}.mask.bed")

        // Using all callable loci, worth knowing another option would be to try use the callable loci for each genome
        // I think this should just be the fasta used for variant calling 
        // i think complement is not using the exons in the bed12 to mask just the intervals
        // can i convert the bed12 here to normal bed ?
        script:
        """
        sort -Vk1 ${callable_bed} > sorted.bed
        bedtools complement -i sorted.bed -g ${genome_file} > ${meta}.mask.bed
        """
}

process get_samples{

        input:
        path(vcf)

        output:
        stdout

        script:
        """
        bcftools query -l ${vcf} 
        """
}

process remove_missing_vcf {

        input:
        val(meta)
        path(vcf)

        output:
        tuple path("${meta}.no_missing.vcf.gz"), path("${meta}.no_missing.vcf.gz.csi")

        script:
        """
        bcftools index -c ${vcf} -o ${vcf}.csi
        bcftools filter -O z --include "N_MISSING=0" ${vcf} > ${meta}.no_missing.vcf.gz
        bcftools index -c ${meta}.no_missing.vcf.gz -o ${meta}.no_missing.vcf.gz.csi
        """
}

process generate_loci {

        input:
        val(meta)
        tuple val(bed_meta), path(mask_bed)
        path(fasta)
        tuple path(vcf), path(vcf_index)

        output:
        tuple val(meta), val(bed_meta), path("${meta}.${bed_meta}.snp.1.fasta"), path("${meta}.${bed_meta}.snp.2.fasta")

        script:
        """
        bcftools consensus -f ${fasta} -m ${mask_bed} -o ${meta}.${bed_meta}.snp.1.fasta -H 1 -s ${meta} ${vcf}
        bcftools consensus -f ${fasta} -m ${mask_bed} -o ${meta}.${bed_meta}.snp.2.fasta -H 2 -s ${meta} ${vcf}
        """
}

process generate_effective_fasta_AGAT {

        input:
        tuple val(meta), val(fasta_meta), path(consensus_fasta_1), path(consensus_fasta_2)
        tuple val(gff_meta), path(gff)

        output:
        tuple val(meta), val(fasta_meta), path("${meta}.${fasta_meta}.snp.1.cds.fasta"), path("${meta}.${fasta_meta}.snp.2.cds.fasta")

        script:
        """
        agat_sp_extract_sequences.pl --gff ${gff} --fasta ${consensus_fasta_1} -t exon --merge -o ${meta}.${fasta_meta}.snp.1.cds.fasta
        agat_sp_extract_sequences.pl --gff ${gff} --fasta ${consensus_fasta_2} -t exon --merge -o ${meta}.${fasta_meta}.snp.2.cds.fasta
        """
}

// I think it seems ok for now, its at least given me individuals with two haplotypes different from the reference
// but i havent seen a heterozygous individual yet
// I might go with this for now and if in the future it hasnt applied the haplotypes right i can fix it then

process orthofinder {
        cpus 16

        input:
        tuple val(meta), path(protein_fastas, stageAs: "fastas/*")

        output:
        path("fastas/OrthoFinder/Results_*/Single_Copy_orthologue_sequences/*")

        script:
        """
        orthofinder -f fastas -t ${task.cpus} -a ${task.cpus}
        """
}

// in the flows need to split out the orthofinder dir and get it to run mafft on each file
process mafft {

        input:
        path(fasta)

        output:
        path("")

        script:
        """
        mafft ${fasta} > ${fasta.baseName}.mafft.fa
        """
}

// here i might need a script to rename the aligned protein fasta headers and 
// rename the haplotyped fasta headers to match
// then i need to take each protein alignment and grep all the sample haplotyped fastas
// to get all the translatorx input fastas
// then give these to translatorx
// might be easier to make the files up to this stage and then look at them

process translatorx {

        input:
        path(prot_fasta)
        path(nuc_fasta)

        output:
        path("${nuc_fasta.baseName}.tlx.fa")

        script:
        """
        translatorx -i ${nuc_fasta} -a ${prot_fasta} -o ${nuc_fasta.baseName}.tlx.fa
        """
}

process orthodiver {

        input:
        path(SCO_dir)

        output:
        path("results_dir")

        script:
        """
        orthodiver.py -d ${SCO_dir} -A -B -o results_dir
        """
}