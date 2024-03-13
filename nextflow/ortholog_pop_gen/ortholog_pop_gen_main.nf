// Genome file generation command
// samtools faidx species.fasta && cut -f1,2 species.fasta.fai | sort -Vk1 > species.genomefile

// Species 1 data files
params.species_1 = "allacma_fusca"
params.genome_file_1 = "$launchDir/data/results/genomes/allacma_fusca/allacma_fusca.genomefile.txt"
params.vcf_file_1 = "$launchDir/data/results/var_call/allacma_fusca/allacma_fusca.hard_filtered.sorted.vcf.gz"
params.callable_bed_1 = "$launchDir/data/results/var_call/allacma_fusca/allacma_fusca.callable.bed"
params.annot_1 = "$launchDir/data/results/braker3/allacma_fusca/braker.gtf"
params.cds_fasta_1 = "$launchDir/data/results/braker3/allacma_fusca/braker.codingseq"
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

include { generate_haplotypes_flow as generate_haplotypes_flow_1 } from './ortholog_pop_gen_flows.nf'
include { generate_haplotypes_flow as generate_haplotypes_flow_2 } from './ortholog_pop_gen_flows.nf'
include { ortholog_pop_gen_flow } from './ortholog_pop_gen_flows.nf'

workflow {
        generate_haplotypes_flow_1(params.species_1, params.genome_file_1, params.vcf_file_1, params.callable_bed_1, params.annot_1 params.cds_bed_1, params.cds_fasta_1)
        // generate_haplotypes_flow_2(params.species_2, params.genome_file_2, params.vcf_file_2, params.callable_bed_2, params.annot_2, params.cds_bed_2, params.cds_fasta_2)
        // ortholog_pop_gen_flow(generate_haplotypes_flow_1.out, generate_haplotypes_flow_2.out, params.prot_fasta_1, params.prot_fasta_2)
}

// mamba activate ortholog_pop_gen
// mamba install -c conda-forge -c bioconda seqtk samtools bedtools bcftools docopt orthofinder mafft translatorx tqdm agat