include { seqtk_get_callable_cds; mask_fasta; get_samples; generate_loci; generate_effective_fastas; orthofinder; mafft; translatorx; orthodiver} from './ortholog_pop_gen_tasks.nf'

workflow generate_haplotypes_flow {
        take:
          species
          genome_file
          vcf
          callable_bed
          cds_bed
          cds_fasta

        main:
          seqtk_get_callable_cds(species, cds_fasta, callable_bed)
          mask_fasta(seqtk_get_callable_cds.out, callable_bed, genome_file)
          get_samples(vcf)
          sample_ch = Channel.value(file(get_samples.out).readlines())
          generate_loci(sample_ch, mask_fasta.out, vcf)
          generate_effective_fastas(generate_loci.out, cds_bed)

        emit:
          generate_effective_fastas.out
}

// assuming 2 protein files in prot_dir for now. should generalise for any number of samples
// maybe can just do all pairwise combinations of samples but take proteins by species
workflow ortholog_pop_gen_flow {
        take:
          nuc_fastas_1
          nuc_fastas_2
          prot_fasta_1
          prot_fasta_2
        main:
          orthofinder(prot_fasta_1, prot_fasta_2) // might have to do a stageas for these
          mafft(orthofinder.out) // need to figure out how to split orthofinder out dir to single files for mafft
          translatorx()
          orthodiver()
}

// should be able to make a workflow that doesnt need generate haplotypes and so just gets dxy not pi
workflow ortholog_div_flow {
        take:
          prot_dir
        main:
          orthofinder(prot_dir)
}