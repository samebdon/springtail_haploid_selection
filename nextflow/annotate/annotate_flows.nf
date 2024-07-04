include {unmask_genome; earlGrey; mask_genome; braker2} from './annotate_tasks.nf'

workflow braker2_flow {
        take:
          meta
          genome
          prot_seq

        main:
          unmask_genome(meta, genome)
          earlGrey(unmask_genome.out)
          mask_genome(genome, earlGrey.out.repeat_bed)
          braker2(mask_genome.out, prot_seq)
}

workflow braker2_only_flow {
        take:
          meta
          genome
          prot_seq

        main:
          braker2(meta.combine(genome), prot_seq)
}