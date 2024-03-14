// Output of the pairwise comparisons could be a matrix of sample comparisons at each gene at 0d and 4d sites, where the diagonal would be diversity


process get_best_cds_bed {

        input:
        val(meta)
        val(annotation)

        output:
        tuple val(meta), path("${meta}.cds.bed")

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
        val(meta)
        path(callable_bed)
        path(genome_file)

        output:
        tuple val(meta), path("${meta}.mask.bed")

        // Now masking only callable and at the end cutting out CDS
        // I think im using wrong callable file, check whats used for vcf
        // Im going to try use all callable loci, worth knowing another option would be to try use the callable loci for each genome
        // I think this should just be the fasta used for variant calling 

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

        // gets sample names in vcf to create a channel
        script:
        """
        bcftools query -l ${vcf} 
        """
}

process generate_loci {

        input:
        val(meta)
        tuple val(bed_meta), path(mask_bed)
        path(fasta)
        path(vcf)

        output:
        tuple val(meta), val(bed_meta), path("${meta}.${bed_meta}.snp.1.callable.fasta"), path("${meta}.${fasta_meta}.snp.2.callable.fasta")

        // this is tabixing the vcf each time a sample is run should make its own process really
        // or include index

        // GENERATE LOCI FROM WHOLE GENOME MASKED CDS
        // Got an error because of variants outside of CDS
        // I could give this all callable loci because thats what the vcf went with
        // Surely though i can get it to ignore masked sites though
        // Should also double check i guess which callable sites file is definitely correct
        script:
        """
        tabix -p vcf ${vcf}
        bcftools consensus -f ${fasta} -m ${mask_bed} -o ${meta}.${bed_meta}.snp.1.fasta -H 1 -s ${meta} ${vcf}
        bcftools consensus -f ${fasta} -m ${mask_bed} -o ${meta}.${bed_meta}.snp.2.fasta -H 2 -s ${meta} ${vcf}
        """
}

process generate_effective_fastas {

        input:
        tuple val(meta), val(fasta_meta), path(consensus_fasta_1), path(consensus_fasta_2)
        tuple val(cds_meta), path(cds_bed)

        output:
        tuple val(meta), val(fasta_meta), path("${meta}.${fasta_meta}.snp.1.effective.cds.fasta"), path("${meta}.${fasta_meta}.snp.2.effective.cds.fasta")

        // removed interleaved since i dont think i need it. if i do include here at a future date
        // EXTRACT CDS USING CDS BED FROM WHOLE GENOME HAPLOTYPES
        // The other option is instead of doing this on the whole genome is to apply
        // the cds location to fasta headers
        // this way is more generalisable though i could make it for intergenic regions quite easily
        script:
        """
        bedtools getfasta -name -s -fi ${consensus_fasta_1} -bed ${cds_bed} -fo ${meta}.${fasta_meta}.snp.1.cds.fasta
        bedtools getfasta -name -s -fi ${consensus_fasta_2} -bed ${cds_bed} -fo ${meta}.${fasta_meta}.snp.2.cds.fasta
        """
}

process orthofinder {
        cpus 64

        input:
        path(fasta_dir)

        output:
        path("${fasta_dir}/OrthoFinder/Results_*/Single_Copy_orthologue_sequences/*")

        script:
        """
        orthofinder -f ${fasta_dir} -t ${task.cpus} -a ${task.cpus}
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