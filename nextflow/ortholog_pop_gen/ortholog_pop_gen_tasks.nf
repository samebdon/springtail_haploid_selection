// Output of the pairwise comparisons could be a matrix of sample comparisons at each gene at 0d and 4d sites, where the diagonal would be diversity

process seqtk_get_callable_cds {
        // get cds in callable regions

        input:
        val(meta)
        path(cds_fasta)
        path(callable_bed)

        output:
        tuple val(meta), path("${meta}.callable.cds.fasta")

        script:
        """
        seqtk subseq ${cds_fasta} <(cut -f1 ${callable_bed} | sort | uniq) > ${meta}.callable.cds.fasta
        """
}

process mask_fasta {

        input:
        tuple val(meta), path(callable_cds_fasta)
        path(callable_bed)
        path(genome_file)

        output:
        tuple val(meta), path("${meta}.callable.cds.masked.fasta")

        script:
        """
        bedtools maskfasta -fi ${callable_cds_fasta} -bed <(bedtools complement -i <(sort -Vk1 ${callable_bed}) -g ${genome_file}) -fo ${meta}.callable.cds.masked.fasta -mc - 
        """
}

process get_samples{

        input:
        path(vcf)

        output:
        path("${vcf.baseName}.sample_IDs.txt")

        script:
        """
        bcftools query -l ${vcf} > ${vcf.baseName}.sample_IDs.txt
        """
}

process generate_loci {

        input:
        val(meta)
        tuple val(fasta_meta), path(masked_fasta)
        path(vcf)

        output:
        tuple val(meta), val(fasta_meta), path("${meta}.${fasta_meta}.snp.1.callable.fasta"), path("${meta}.${fasta_meta}.snp.2.callable.fasta")

        script:
        """
        bcftools consensus -f ${masked_fasta} -o ${meta}.${fasta_meta}.snp.1.callable.fasta -H 1 -s ${meta} ${vcf}
        bcftools consensus -f ${masked_fasta} -o ${meta}.${fasta_meta}.snp.2.callable.fasta -H 2 -s ${meta} ${vcf}
        """
}

process generate_effective_fastas {

        input:
        tuple val(meta), val(fasta_meta), path(consensus_fasta_1), path(consensus_fasta_2)
        path(cds_bed)

        output:
        tuple val(meta), val(fasta_meta), path("${meta}.${fasta_meta}.snp.1.effective.cds.fasta"), path("${meta}.${fasta_meta}.snp.2.effective.cds.fasta")

        // removed interleaved since i dont think i need it. if i do include here at a future date
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