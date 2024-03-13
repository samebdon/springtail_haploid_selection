include { get_best_cds_bed; get_best_pep_fasta; get_callable_cds_bed; seqtk_get_callable_cds; mask_fasta; get_samples; generate_loci; generate_effective_fastas; orthofinder; mafft; translatorx; orthodiver} from './ortholog_pop_gen_tasks.nf'

workflow generate_haplotypes_flow {
        take:
          species
          genome_file
          vcf
          callable_bed
          annotation
          cds_fasta
          pep_fasta

        main:
          get_best_cds_bed(species, annotation)
          get_best_pep_fasta(get_best_cds_bed.out, pep_fasta)
          get_callable_cds_bed(get_best_cds_bed.out, callable_bed)
          seqtk_get_callable_cds(get_callable_cds_bed.out, cds_fasta)
          mask_fasta(seqtk_get_callable_cds.out, get_callable_cds_bed.out, genome_file)
          get_samples(vcf)
          generate_loci(get_samples.out.splitText( by: 1 ).map{it -> it.trim()}, mask_fasta.out, vcf)
          generate_effective_fastas(generate_loci.out, get_best_cds_bed.out)

        emit:
          generate_effective_fastas.out.collect()
          get_best_pep_fasta.out
}

// assuming 2 protein files in prot_dir for now. should generalise for any number of samples
// maybe can just do all pairwise combinations of samples but take proteins by species
workflow ortholog_pop_gen_flow {
        take:
          nuc_fastas_1
          prot_fasta_1
          nuc_fastas_2
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