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

//workflow infer_orthology {
//        take:
//          prot_fastas
//
//        main:
//          orthofinder(prot_fastas)
//          select_orthogroups()
//
//        emit:
//          select_orthogroups.out
//}

// assuming 2 protein files in prot_dir for now. should generalise for any number of samples
workflow orthodiver_flow {
        take:
          hap_fastas_1
          prot_fasta_1
          hap_fastas_2
          prot_fasta_2
        main:

          orthofinder(prot_fasta_1, prot_fasta_2) // edit to just take one pile of fastas, join?
          // filter_orthogroups() into mafft (filter_orthogroups.out.flatten())
          mafft(orthofinder.out.flatten())
          dupe_prot_fasta(mafft.out)
          get_orthogroup_haps(mafft.out, hap_fastas_1, hap_fastas_2)
          tlx_in_ch = dupe_prot_fasta.out.join(get_orthogroup_haps.out).map { it -> [it[2], [it[1]]].combinations() }.flatten().collate(2)
          translatorx(tlx_in_ch)
          orthodiver(translatorx.out.groupTuple())
}


// can just join up all protein fasta inputs into one in channel
// should be able to make a workflow that doesnt need generate haplotypes and so just gets dxy not pi
// I could make alternative workflows that run orthofinder with an arbitrary number of protein fastas
// Then could come up with criteria from sampling orthogroups from these and send the overall set to a set of comparisons
// maybe the thing to do is split out the orthology inference from the alignment flow
// then this last bit is only relevant to a pair and most of it can stay the same
// can control orthology selection and comparisons in previous step
// first step of alignment and pop gen can be selecting relevant orthologs from the orthogroups output
// maybe dont want to edit this yet cos its a lot of unnecessary redoing
// can do it when i next need to go through the pipeline as a whole
// this should just make edits and future pipelines more straight forward since the beginning and end should already be stable
// just separate out the unstable bit for easy editing