// Output of the pairwise comparisons could be a 
// matrix of sample comparisons at each gene at 0d and 4d sites, 
// where the diagonal would be diversity?


process get_best_cds_bed {
        publishDir params.outdir, mode:'copy'
        memory '4G'
        scratch true

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
        path("${meta}.best.pep.fasta")

        // Should be fine, check best.pep.fasta looks ok
        script:
        """
        cut -f4 ${cds_bed} > prots.lst
        fastaqual_select.pl -f ${pep_fasta} | cut -f1 -d" " > selected.pep
        grep --no-group-separator -A1 -wFf prots.lst selected.pep | sed "s/>/&${meta}./g" > ${meta}.best.pep.fasta
        """
}

process get_callable_cds_bed {
        memory '4G'

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
        // see if you can run tiny jobs like this on the execution node

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
        bcftools index -c ${vcf}
        bcftools plugin fill-tags -Oz ${vcf} -t 'F_MISSING' > tagged.vcf.gz
        bcftools filter -O z --include "F_MISSING=0" tagged.vcf.gz > ${meta}.no_missing.vcf.gz
        bcftools index -c ${meta}.no_missing.vcf.gz
        """
}

process generate_loci {

        input:
        val(meta)
        tuple val(bed_meta), path(mask_bed)
        path(fasta)
        path(vcf)

        output:
        tuple val(meta), val(bed_meta), path("${meta}.${bed_meta}.snp.1.fasta"), path("${meta}.${bed_meta}.snp.2.fasta")

        script:
        """
        bcftools index -c ${vcf}
        bcftools consensus -f ${fasta} -m ${mask_bed} -o ${meta}.${bed_meta}.snp.1.fasta -H 1 -s ${meta} ${vcf}
        bcftools consensus -f ${fasta} -m ${mask_bed} -o ${meta}.${bed_meta}.snp.2.fasta -H 2 -s ${meta} ${vcf}
        """
}

process generate_effective_fasta_AGAT {
        memory '4G'
        scratch true

        input:
        tuple val(meta), val(fasta_meta), path(consensus_fasta_1), path(consensus_fasta_2)
        tuple val(gff_meta), path(gff)

        output:
        tuple path("${meta}.${fasta_meta}.snp.1.cds.fasta"), path("${meta}.${fasta_meta}.snp.2.cds.fasta")

        script:
        """
        agat_sp_extract_sequences.pl --gff ${gff} --fasta ${consensus_fasta_1} -t exon --merge -o ${meta}.${fasta_meta}.snp.1.cds.fasta.tmp
        awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < ${meta}.${fasta_meta}.snp.1.cds.fasta.tmp | tail -n +2 > ${meta}.${fasta_meta}.snp.1.cds.fasta
        agat_sp_extract_sequences.pl --gff ${gff} --fasta ${consensus_fasta_2} -t exon --merge -o ${meta}.${fasta_meta}.snp.2.cds.fasta.tmp
        awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < ${meta}.${fasta_meta}.snp.2.cds.fasta.tmp | tail -n +2 > ${meta}.${fasta_meta}.snp.2.cds.fasta
        """
}

// I think it seems ok for now, its at least given me individuals with two haplotypes different from the reference
// but i havent seen a heterozygous individual yet
// I might go with this for now and if in the future it hasnt applied the haplotypes right i can fix it then

process orthofinder {
        cpus 16
        scratch true

        // include all springtails to improve inference

        input:
        path(prot_fastas, stageAs: "fastas/*")

        output:
        path("results/*"), emit: all
        path("results/*/Single_Copy_Orthologue_Sequences/*"), emit: sco

        script:
        """
        orthofinder -f fastas -t ${task.cpus} -a ${task.cpus} -o results
        """
}

// works with 2 sp for now
process get_SCO_genes{
        publishDir params.outdir, mode:'copy'

        input:
        path(orthofinder_results, stageAs: "orthofinder_results/*")

        output:
        path('SC_orthogroups.txt'), emit: sc_orthogroups
        tuple path('*.sp1.SCO_genes.txt'), path('*.sp2.SCO_genes.txt'), emit: genes

        script:
        """
        parallel -j1 'grep "{}" orthofinder_results/*/Orthogroups/Orthogroups.txt >> SC_orthogroups.txt' :::: orthofinder_results/*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt
        
        SP1="\$(cat SC_orthogroups.txt | cut -f2-2 -d' ' | cut -f-1 -d'.' | head -n1)"
        cat SC_orthogroups.txt | cut -f2-2 -d' ' | cut -f2- -d'.' > \$SP1.sp1.SCO_genes.txt

        SP2="\$(cat SC_orthogroups.txt | cut -f3-3 -d' ' | cut -f-1 -d'.' | head -n1)"
        cat SC_orthogroups.txt | cut -f3-3 -d' ' | cut -f2- -d'.' > \$SP2.sp2.SCO_genes.txt
     
        """
}

process filter_annotation{

        input:
        tuple path(sp1_genes), path(sp2_genes)
        path(annotation_1, stageAs: "sp1_annotation/*")
        path(annotation_2, stageAs: "sp2_annotation/*")

        output:
        path("${sp1_genes.baseName}.SCOs.gtf"), emit: sp1
        path("${sp2_genes.baseName}.SCOs.gtf"), emit: sp2

        script:
        """
        parallel -j1 'grep "{}" sp1_annotation/annotation.gtf >> ${sp1_genes.baseName}.SCOs.gff' :::: ${sp1_genes}
        parallel -j1 'grep "{}" sp2_annotation/annotation.gtf >> ${sp2_genes.baseName}.SCOs.gff' :::: ${sp2_genes}
        """
}

// TO DO
process filter_orthogroups{

        input:
        val(species_1)
        val(species_2)
        path(orthogroup_fastas)

        output:

        script:
        """    
        echo ${species_1}

        """
}

process mafft {

        input:
        path(fasta)

        output:
        tuple val("${fasta.baseName}"), path("${fasta.baseName}.mafft.fa")

        script:
        """
        mafft --thread ${task.cpus} ${fasta} > ${fasta.baseName}.mafft.fa
        """
}

process dupe_prot_fasta {

        input:
        tuple val(meta), path(prot_fasta)

        output:
        tuple val(meta), path("${meta}.mafft.happed.fa")

        script:
        """
        awk '/^>/ {printf("\\n%s\\n",\$0);next; } { printf("%s",\$0);}  END {printf("\\n");}' < ${prot_fasta} | tail -n +2 > ${meta}.mafft.single_line.fa
        duplicate_prot_aln.sh ${meta}.mafft.single_line.fa
        """
}

process mafft_batch {
        cpus 4
        scratch true

        input:
        path(fastas, stageAs: "fastas/*")

        output:
        path("out/*.mafft.happed.fa")

        // could change parallelising to be per process 
        script:
        """
        mkdir out
        for fasta in fastas/*
        do
                ORTHOGROUP="\$(echo \$fasta | cut -f2- -d'/' | cut -f-1 -d'.')"
                mafft --thread ${task.cpus} \$fasta > \$ORTHOGROUP.mafft.fa
                cat \$ORTHOGROUP.mafft.fa | seqtk seq > \$ORTHOGROUP.mafft.single_line.fa
                duplicate_prot_aln.sh \$ORTHOGROUP.mafft.single_line.fa
        done
        """
}

// can I remove non informative genes before running get hap, renaming, and pair fastas?

process get_orthogroup_haps {
        scratch true

        input:
        tuple val(meta), path(prot_fasta)
        path(sp1_fastas, stageAs: "sp1_fastas/*")
        path(sp2_fastas, stageAs: "sp2_fastas/*")

        output:
        tuple val(meta), path("hap_fasta_pairs/*.unaln.fa"), optional: true

        script:
        """
        SP1="\$(cat ${prot_fasta} | grep '>' | head -n 1 | cut -d '>' -f2- | cut -d'.' -f-1)"
        SP2="\$(cat ${prot_fasta} | grep '>' | tail -n 1 | cut -d '>' -f2- | cut -d'.' -f-1)"
        SP1_PROT="\$(cat ${prot_fasta} | grep '>' | head -n 1 | cut -d '>' -f2- | cut -d'.' -f2-)"
        SP2_PROT="\$(cat ${prot_fasta} | grep '>' | tail -n 1 | cut -d '>' -f2- | cut -d'.' -f2-)"

        mkdir hap_fastas
        mkdir hap_fastas_rn

        get_hap.sh sp1_fastas \$SP1_PROT hap_fastas ${meta}
        get_hap.sh sp2_fastas \$SP2_PROT hap_fastas ${meta}

        rename_hap_fastas.sh hap_fastas hap_fastas_rn

        pair_fastas.py -i hap_fastas_rn -o hap_fasta_pairs -a \$SP1 -b \$SP2
        """
}

// could even combine this with translatorx pair to really reduce symlinking
// this bit is slow be nice to parallelise the loop
process get_orthogroup_haps_batch {
        scratch true

        input:
        path(prot_fastas, stageAs: "prot_fastas/*")
        path(sp1_fastas, stageAs: "sp1_fastas/*")
        path(sp2_fastas, stageAs: "sp2_fastas/*")

        output:
        path("out.*"), optional: true

        script:
        """

        mkdir hap_fastas
        mkdir hap_fastas_rn

        for prot_fasta in prot_fastas/*
        do
                ORTHOGROUP="\$(echo \$prot_fasta | cut -f2- -d'/' | cut -f-1 -d'.')"
                SP1="\$(cat \$prot_fasta | grep '>' | head -n 1 | cut -d '>' -f2- | cut -d'.' -f-1)"
                SP2="\$(cat \$prot_fasta | grep '>' | tail -n 1 | cut -d '>' -f2- | cut -d'.' -f-1)"
                SP1_PROT="\$(cat \$prot_fasta | grep '>' | head -n 1 | cut -d '>' -f2- | cut -d'.' -f4-)"
                SP2_PROT="\$(cat \$prot_fasta | grep '>' | tail -n 1 | cut -d '>' -f2- | cut -d'.' -f4-)"

                get_hap.sh sp1_fastas \$SP1_PROT hap_fastas \$ORTHOGROUP
                get_hap.sh sp2_fastas \$SP2_PROT hap_fastas \$ORTHOGROUP

                rename_hap_fastas.sh hap_fastas hap_fastas_rn

                pair_fastas.py -i hap_fastas_rn -o out -a \$SP1 -b \$SP2

                rm hap_fastas/*
                rm hap_fastas_rn/*
        done
        """
}

// currently with pair fastas im just comparing different species
// 0d diversity matrix is not that interesting, it shouldnt vary that much
// what would be better is to get the SFS. how do I do that?
// be good to get the sfs for genes and 0d and 4d sites but also whole genome intergenically
// dont need any alignments for this, should just be able to for genes continue off the initial workflow
// then just use the vcf and beds to subset gene third codon positions like in LG het
// this could be same but with intergenic regions to get a genome wide sfs rather than per gene?

// can i run all of this in get orthogroup haps?

process translatorx {
        scratch true

        input:
        tuple path(hap_fasta), path(prot_fasta)

        output:
        tuple val("\$SAMPLE_1_\$SAMPLE_2"), path("*.tlx.fa")

        script:
        """
        OUT_PREFIX="\$(ls *.unaln.fa| cut -d'.' -f-5)"
        SAMPLE_1="\$(ls *.unaln.fa | cut -d'.' -f2-2)"
        SAMPLE_2="\$(ls *.unaln.fa | cut -d'.' -f4-4)"

        sed -e "s/sample_1/\$SAMPLE_1/g" ${prot_fasta} > ${prot_fasta.baseName}.rn.fa
        sed -i -e "s/sample_2/\$SAMPLE_2/g" ${prot_fasta.baseName}.rn.fa

        translatorx -i ${hap_fasta} -a ${prot_fasta.baseName}.rn.fa -o \$OUT_PREFIX.tlx.fa
        """
}


process translatorx_pair {
        memory '4G'
        scratch true

        input:
        path(prot_fastas, stageAs: "prot_fastas/*")
        path(pair_fastas)

        output:
        tuple val("\$SAMPLE_1.\$SAMPLE_2"), path("tlx_fastas/*.nt_ali.fasta")

        script:
        """
        mkdir tlx_fastas

        for prot_fasta in prot_fastas/*
        do
                ORTHOGROUP="\$(echo \$prot_fasta | cut -f2- -d'/' | cut -f-1 -d'.')"
                PAIR_FASTA="${pair_fastas}/\$ORTHOGROUP*"
                OUT_PREFIX="\$(echo \$PAIR_FASTA | cut -f2- -d'/' | cut -d'.' -f-5)"
                SAMPLE_1="\$(echo \$PAIR_FASTA | cut -f2- -d'/' | cut -d'.' -f2-2)"
                SAMPLE_2="\$(echo \$PAIR_FASTA | cut -f2- -d'/' | cut -d'.' -f4-4)"

                sed -e "s/sample_1/\$SAMPLE_1/g" \$prot_fasta > \$prot_fasta.rn.fa
                sed -i -e "s/sample_2/\$SAMPLE_2/g" \$prot_fasta.rn.fa

                translatorx -i \$PAIR_FASTA -a \$prot_fasta.rn.fa -o tlx_fastas/\$OUT_PREFIX
        done
        """
}

// is it worth including a multiple sequence alignment QC step here? trimAl?

process orthodiver {
        memory '4G'

        input:
        tuple val(meta), path(orthlg_fastas, stageAs: "fastas/*")

        output:
        tuple path("*.0d_pi_by_locus.txt"), path("*.4d_pi_by_locus.txt")

        script:
        """
        SAMPLE_1="\$(ls fastas | head -n 1 | cut -d'.' -f2-2)"
        SPECIES_1="\$(ls fastas/* | head -n 1 | cut -d'.' -f3-3)"
        SAMPLE_2="\$(ls fastas/* | head -n 1 | cut -d'.' -f4-4)"
        SPECIES_2="\$(ls fastas/* | head -n 1 | cut -d'.' -f5-5)"
 
        orthodiver.py -d fastas -A \$SPECIES_1.\$SAMPLE_1 -B \$SPECIES_2.\$SAMPLE_2 -o \$SAMPLE_1.\$SAMPLE_2.results
        """
}

process agg_orthodiver {
        publishDir params.outdir, mode:'copy'
        memory '4G'

        input:
        path(pi_by_locus, stageAs: "orthodiver_results/*")

        output:
        path("*.orthodiver_agg.tsv")

        script:
        """
        agg_orthodiver.py -i orthodiver_results
        """
}
