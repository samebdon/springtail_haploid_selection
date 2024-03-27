// Output of the pairwise comparisons could be a matrix of sample comparisons at each gene at 0d and 4d sites, where the diagonal would be diversity


process get_best_cds_bed {

        input:
        val(meta)
        path(annotation)

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
        tuple path("${meta}.${fasta_meta}.snp.1.cds.fasta"), path("${meta}.${fasta_meta}.snp.2.cds.fasta")

        script:
        """
        agat_sp_extract_sequences.pl --gff ${gff} --fasta ${consensus_fasta_1} -t exon --merge -o ${meta}.${fasta_meta}.snp.1.cds.fasta.tmp
        awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < ${meta}.${fasta_meta}.snp.1.cds.fasta.tmp > ${meta}.${fasta_meta}.snp.1.cds.fasta
        agat_sp_extract_sequences.pl --gff ${gff} --fasta ${consensus_fasta_2} -t exon --merge -o ${meta}.${fasta_meta}.snp.2.cds.fasta.tmp
        awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < ${meta}.${fasta_meta}.snp.2.cds.fasta.tmp > ${meta}.${fasta_meta}.snp.2.cds.fasta
        """
}
// I've got a mixture of soft and hard masking in here, Should double check what to do and go with one


// I think it seems ok for now, its at least given me individuals with two haplotypes different from the reference
// but i havent seen a heterozygous individual yet
// I might go with this for now and if in the future it hasnt applied the haplotypes right i can fix it then

process orthofinder {
        cpus 16

        input:
        tuple val(meta_1), path(prot_fasta_1, stageAs: "fastas/*")
        tuple val(meta_2), path(prot_fasta_2, stageAs: "fastas/*")

        output:
        path("fastas/OrthoFinder/Results_*/Single_Copy_Orthologue_Sequences/*")

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
        tuple val("${fasta.baseName}"), path("${fasta.baseName}.mafft.fa")

        script:
        """
        mafft ${fasta} > ${fasta.baseName}.mafft.fa
        """
}

process dupe_prot_fasta {

        input:
        tuple val(meta), path(prot_fasta)

        output:
        tuple val(meta), path("${meta}.mafft.happed.fa")

        script:
        """
        duplicate_prot_aln.sh ${prot_fasta}
        """
}

// Irritatingly it looks like the protein file needs identical fasta headers to the
// nucleotide file
// Maybe the way to do this is in the translatorX process take the sample file I have
// and impute the names into the protein file
// Not sure how to do this
// Could get sample 1 and 2 and sed it in before running

// Once renamed, generate all pairwise combinations of haplotypes for alignment with translatorx
// feels liike a python job to do this step

// forward prot fasta and each nuc fasta to channels

// Not sure how to sort out the nuc channel to give to this, need all samples at once

// This should output all pairwise combinations of sample haplotypes per orthogroup


// run this now and then add the bit that gets the pairwise combos
process get_orthogroup_haps {

        input:
        tuple val(meta), path(prot_fasta)
        path(sp1_fastas, stageAs: "sp1_fastas/*")
        path(sp2_fastas, stageAs: "sp2_fastas/*")

        output:
        tuple val(meta), path("*.unaln.fa")

        script:
        """
        SP1_PROT="\$(cat ${prot_fasta} | grep '>' | head -n 1 | cut -d '>' -f2- | cut -d'.' -f2-)"
        SP2_PROT="\$(cat ${prot_fasta} | grep '>' | tail -n 1 | cut -d '>' -f2- | cut -d'.' -f2-)"

        mkdir hap_fastas
        mkdir hap_fastas_rn

        get_hap.sh sp1_fastas/* \$SP1_PROT hap_fastas ${meta}
        get_hap.sh sp2_fastas/* \$SP2_PROT hap_fastas ${meta}

        rename_hap_fastas.sh hap_fastas/* hap_fastas_rn
        """

}

process translatorx {

        input:
        tuple val(meta), path(prot_fasta), path(nuc_fasta)

        output:
        path("${nuc_fasta.baseName}.tlx.fa")

        script:
        """
        translatorx -i ${nuc_fasta} -a ${prot_fasta} -o ${meta}.tlx.fa
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