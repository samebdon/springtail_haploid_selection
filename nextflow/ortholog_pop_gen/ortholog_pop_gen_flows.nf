include { get_best_cds_bed; get_best_pep_fasta; get_callable_cds_bed; make_genome_file; get_mask_bed; get_samples; remove_missing_vcf; generate_loci; generate_effective_fasta_AGAT; orthofinder; mafft; dupe_prot_fasta; get_orthogroup_haps; translatorx; orthodiver} from './ortholog_pop_gen_tasks.nf'

workflow gen_haps_flow {
        take:
          species
          genome_fasta
          vcf
          callable_bed
          annotation
          pep_fasta

        main:
          get_best_cds_bed(species, annotation)
          get_best_pep_fasta(get_best_cds_bed.out.bed, pep_fasta)
          get_callable_cds_bed(get_best_cds_bed.out.bed, callable_bed)
          make_genome_file(genome_fasta)
          get_mask_bed(get_callable_cds_bed.out, make_genome_file.out)
          get_samples(vcf)
          remove_missing_vcf(species, vcf)
          generate_loci(get_samples.out.splitText( by: 1 ).map{it -> it.trim()}, get_mask_bed.out, genome_fasta, remove_missing_vcf.out)
          generate_effective_fasta_AGAT(generate_loci.out, get_best_cds_bed.out.gff)

        emit:
          generate_effective_fasta_AGAT.out.collect()
          get_best_pep_fasta.out
}

// assuming 2 protein files in prot_dir for now. should generalise for any number of samples
workflow orthodiver_flow {
        take:
          hap_fastas_1
          prot_fasta_1
          hap_fastas_2
          prot_fasta_2
        main:
          orthofinder(prot_fasta_1, prot_fasta_2)
          mafft(orthofinder.out.flatten())
          dupe_prot_fasta(mafft.out)
          get_orthogroup_haps(mafft.out, hap_fastas_1, hap_fastas_2)
          tlx_in_ch = dupe_prot_fasta.out.join(get_orthogroup_haps.out).map { it -> [it[2], [it[1]]].combinations() }.flatten().collate(2)
          translatorx(tlx_in_ch)
          // orthodiver(translatorx.out.collect())
}

// should be able to make a workflow that doesnt need generate haplotypes and so just gets dxy not pi
workflow orthodiver_nopi_flow {
        take:
          prot_dir
        main:
          orthofinder(prot_dir)
}