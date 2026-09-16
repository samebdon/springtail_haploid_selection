calculatePiBed.py -v allacma_fusca.hard_filtered.sorted.vcf.gz -b chromosome_bed_file -a desired part of allacma_fusca.longest_isoforms.0D.bed -n chromosome_name -o resample.n.0D.pi.tsv -l 0D_pi -g allacma_fusca.genomefile
calculatePiBed.py -v allacma_fusca.hard_filtered.sorted.vcf.gz -b chromosome_bed_file -a desired part of allacma_fusca.longest_isoforms.4D.bed -n chromosome_name -o resample.n.4D.pi.tsv -l 4D_pi -g allacma_fusca.genomefile
maybe generate directories with resampled input files for the script and then
run the script on each directory for each chromosome for 0d and 4d 
write a script to sort out the results directories to aggregate the SFSs and reformat to dfe alpha


could be good to compare braker annotation pipeline to repeatmasker
also check how braker2 output looks from nextflow to the manual one just in case


/lustre/scratch126/tol/teams/jaron/data/assemblies_Sanger/arthropods/Sminthurus_viridis/assembly/draft/qlSmiViri2.20241204/qlSmiViri2.20241204.primary.fa.gz
$NF_PATH/annotate/annotate_main.nf

-params-file $NF_PATH/annotate/s_viridis_draft_params.json

module load samtools/1.20--h50ea8bc_0 
samtools sort -n -o data/raw_data/svir/svir.qsort.bam ../../data/assemblies_Sanger/arthropods/Sminthurus_viridis/genomic_data/qlSmiViri2/pacbio/m84098_241129_133649_s2.ccs.bc2012.rmdup.bam
samtools bam2fq -1 svir.1.fastq.gz -2 svir.2.fastq.gz svir.qsort.bam

samtools fastq -@ 14 svir.qsort.bam \
    -1 svir.1.fastq.gz \
    -2 svir.2.fastq.gz \
    -0 /dev/null -s /dev/null -n

rm var call log
bsub var call

nano /lustre/scratch126/tol/teams/jaron/users/sam/nf_pipelines/var_call/svir_params.json

data/results/annotations/sminthurus_viridis/earlgrey/sminthurus_viridis_draft_summaryFiles/sminthurus_viridis_draft.softmasked.fasta 
data/results/annotations/sminthurus_viridis/earlgrey/sminthurus_viridis_draft_summaryFiles/sminthurus_viridis_draft.filteredRepeats.bed

samtools fastq --reference ../../results/annotations/sminthurus_viridis/earlgrey/sminthurus_viridis_draft_summaryFiles/sminthurus_viridis_draft.softmasked.fasta -1 svir.1.fastq.gz -2 svir.2.fastq.gz -0 other.fastq.gz ../../../../../data/assemblies_Sanger/arthropods/Sminthurus_viridis/genomic_data/qlSmiViri2/pacbio/m84098_241129_133649_s2.ccs.bc2012.rmdup.bam