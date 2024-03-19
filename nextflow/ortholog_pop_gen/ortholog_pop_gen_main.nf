// Species 1 data files
params.species_1 = "allacma_fusca"
params.genome_fasta_1 = "$launchDir/data/results/genomes/allacma_fusca/GCA_947179485.1.simple_header.earlGrey_masked.fasta"
params.vcf_1 = "$launchDir/data/results/var_call/allacma_fusca/allacma_fusca.hard_filtered.sorted.vcf.gz"
params.callable_bed_1 = "$launchDir/data/results/var_call/allacma_fusca/allacma_fusca.all_callable.bed"
params.annot_1 = "$launchDir/data/results/braker3/allacma_fusca/braker.gtf"
params.prot_fasta_1 = "$launchDir/data/results/braker3/allacma_fusca/braker.aa"

// Will have to test this on something I have reads for variant calling for. Is this only Dicyrtomina minuta? Not sure if its annotated either
// Run workflow 1 for testing and then when i have that output and I can work on workflow 2 i can sort out these extra data files
// Species 2 data files
params.species_2 = "sminthurides_aquaticus"
// params.genome_file_2 = ".genomefile"
// params.vcf_file_2 = ".vcf.gz"
// params.callable_bed_2 = ".bed"
// params.cds_bed_2 = ".bed"
// params.cds_fasta_2 = ".fasta"
// params.prot_fasta_2 = ".fasta"

// Output Directory
params.outdir = "$launchDir/ortholog_pop_gen/allacma_fusca.vs.sminthurides_aquaticus"

log.info """\
         O R T H O L O G  P O P  G E N   N F   P I P E L I N E    
         ===================================
         Species 1 : ${params.species_1}
         Species 2 : ${params.species_2}
         outdir : ${params.outdir}
         """
         .stripIndent()

include { gen_haps_flow as gen_haps_flow_1 } from './ortholog_pop_gen_flows.nf'
include { gen_haps_flow as gen_haps_flow_2 } from './ortholog_pop_gen_flows.nf'
include { orthodiver_flow } from './ortholog_pop_gen_flows.nf'

workflow {
        gen_haps_flow_1(params.species_1, params.genome_fasta_1, params.vcf_1, params.callable_bed_1, params.annot_1, params.prot_fasta_1)
        gen_haps_flow_1(params.species_2, params.genome_fasta_2, params.vcf_2, params.callable_bed_2, params.annot_2, params.prot_fasta_2)
        // orthodiver_flow(generate_haplotypes_flow_1.out[0], generate_haplotypes_flow_1.out[1], generate_haplotypes_flow_2.out[0], generate_haplotypes_flow_2.out[1])
}

// mamba activate ortholog_pop_gen
// mamba install -c conda-forge -c bioconda seqtk samtools bedtools bcftools docopt orthofinder mafft translatorx tqdm agat