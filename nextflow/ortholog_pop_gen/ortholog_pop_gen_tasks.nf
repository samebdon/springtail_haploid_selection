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

process seqtk_get_callable_cds {
        // get cds in callable regions

        input:
        tuple val(meta), path(callable_cds_bed)
        path(cds_fasta)

        output:
        tuple val(meta), path("${meta}.callable.cds.fasta")

        // This just gets the names of transcripts that are callable for further operations
        // Later need to get the CDS out of the transcripts?
        // Or is it already cds because of the braker fasta so i skip the last step?
        // I THINK CAN REMOVE THIS
        script:
        """
        cut -f4 ${callable_cds_bed} | sort | uniq > callable_cds.lst
        seqtk subseq ${cds_fasta} callable_cds.lst > ${meta}.callable.cds.fasta
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

process mask_fasta {

        input:
        tuple val(meta), path(callable_cds_bed)
        path(genome_fasta)
        path(genome_file)

        output:
        tuple val(meta), path("${meta}.callable.masked.fasta")

        // This takes those cds's with callable regions from seqtk and
        // masks the non callable regions using the callable bed
        // is nothing getting masked? need to check
        // this should probably mask the whole genome rather than the callable cds
        // fasta though right? If so do i need run seqtk on the whole genome or just the callable regions?
        // can generate the genomefile in the script if im giving the genomefile anyway

        // RUN THIS ON THE GENOME SO MASK EVERYTHING THAT ISNT CALLABLE CDS
        script:
        """
        sort -Vk1 ${callable_cds_bed} > sorted.bed
        bedtools complement -i sorted.bed -g ${genome_file} > complement.bed
        bedtools maskfasta -fi ${genome_fasta} -bed complement.bed -fo ${meta}.callable.masked.fasta -mc - 
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
        tuple val(fasta_meta), path(masked_fasta)
        path(vcf)

        output:
        tuple val(meta), val(fasta_meta), path("${meta}.${fasta_meta}.snp.1.callable.fasta"), path("${meta}.${fasta_meta}.snp.2.callable.fasta")

        // this is tabixing the vcf each time a sample is run should make its own process really
        // or include index
        // This takes the CDS fasta with masked non callable sites and creates each haplotype based on the vcf
        // Im a bit confused how it can know where the variants are with the information given?
        // Need to look at the previous files to see what they looked like
        // It could be that this should actually be applied to the reference genome and then cut the CDS out of it in the next step?

        // GENERATE LOCI FROM WHOLE GENOME MASKED CDS
        script:
        """
        tabix -p vcf ${vcf}
        bcftools consensus -f ${masked_fasta} -o ${meta}.${fasta_meta}.snp.1.callable.fasta -H 1 -s ${meta} ${vcf}
        bcftools consensus -f ${masked_fasta} -o ${meta}.${fasta_meta}.snp.2.callable.fasta -H 2 -s ${meta} ${vcf}
        """
}

// This bit isnt working, producing no output
// I still dont exactly get what this adds from the previous one
// I should probably look at a positive example to figure it out and see exactly if its necessary
// Looks like its looking for chromosome in fasta file, not the transcript
// Double check if i need to put chromosome name in the fasta file to make this work or if i can just skip?
// I've taken the longest transcripts, in the annotation does this include introns or do i need to select the CDSs
// So should this last step be cutting the CDS's out of the transcripts? I think so

// In theory, given we use the CDS file from braker we shouldnt need to cut the cds out of the haplotypes
// Previously i guess this might have been transcripts and then needed to get out the CDS's
// Figure out if this is the case
process generate_effective_fastas {

        input:
        tuple val(meta), val(fasta_meta), path(consensus_fasta_1), path(consensus_fasta_2)
        tuple val(cds_meta), path(cds_bed)

        output:
        tuple val(meta), val(fasta_meta), path("${meta}.${fasta_meta}.snp.1.effective.cds.fasta"), path("${meta}.${fasta_meta}.snp.2.effective.cds.fasta")

        // removed interleaved since i dont think i need it. if i do include here at a future date
        // EXTRACT CDS USING CDS BED FROM WHOLE GENOME HAPLOTYPES
        script:
        """
        bedtools getfasta -name -s -fi ${consensus_fasta_1} -bed ${cds_bed} -fo ${meta}.${fasta_meta}.snp.1.effective.cds.fasta
        bedtools getfasta -name -s -fi ${consensus_fasta_2} -bed ${cds_bed} -fo ${meta}.${fasta_meta}.snp.2.effective.cds.fasta
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