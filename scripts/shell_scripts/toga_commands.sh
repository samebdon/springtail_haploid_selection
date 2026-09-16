toga processing
from toga/processed

grep 'one2one' ../orthology_classification.tsv > SCO.tsv

#Could filter one2ones for longest reference isoforms

cat SCO.tsv | awk '{print $2}' > reference_transcripts.txt
cat SCO.tsv | awk '{print $4}' > query_transcripts.txt

mkdir codon_fastas
mkdir prot_fastas
mkdir nuc_fastas

parallel -j1 "cat ../codon.fasta | grep -A1 '{}' > codon_fastas/{}.codon.fasta; cat ../codon.fasta | grep -A1 '{}' | grep -A1 'REFERENCE' >> one2one.codon.reference.fasta; cat ../codon.fasta | grep -A1 '{}' | grep -A1 'QUERY' >> one2one.codon.query.fasta" :::: query_transcripts.txt

parallel -j1 "cat ../prot.fasta | grep -A1 '{}' > prot_fastas/{}.prot.fasta; cat ../codon.fasta | grep -A1 '{}' | grep -A1 'REFERENCE' >> one2one.prot.reference.fasta; cat ../codon.fasta | grep -A1 '{}' | grep -A1 'QUERY' >> one2one.prot.query.fasta" :::: query_transcripts.txt

parallel -j1 --link "cat ../nucleotide.fasta | grep -A1 '{1}$' > nuc_fastas/{2}.nuc.fasta; cat ../nucleotide.fasta | grep -A1 '{2}' >> nuc_fastas/{2}.nuc.fasta; cat ../nucleotide.fasta | grep -A1 '{1}$' >> one2one.nuc.reference.fasta; cat ../nucleotide.fasta | grep -A1 '{2}' >> one2one.nuc.query.fasta" :::: reference_transcripts.txt query_transcripts.txt

#Get code to calculate identity and coverage for the alignments
#Can filter out crap alignments

#Optional extra alignment block
#Maybe i can skip this by editing the codon alignments which look better
#Calculating script with alignmentscores would be good for filtering out shit ones or knowing whether to improve them

mkdir aln_nuc_fastas
mkdir aln_prot_fastas

#need to change nuc names to same as protein names
#proteins seem already aligned, assume dont need to improve this

translatorx -i nuc.fasta -o nuc.aln.fasta -a prot.aln.fasta
translatorx -i nuc_fastas/g20106.t1.2.nuc.fasta -o aln_nuc_fastas/g20106.t1.2.nuc.aln.fasta -a prot_fastas/g20106.t1.2.prot.fasta


#CALCULATING DNDS
#For now, can do by editing codon alignment files
#For future, better to edit orthodiver script to take 2 haplotypes rather than 4, this can run on anything edited previously

mkdir orthodiver_input
mkdir orthodiver_output

python prep_orthodiver.py

python orthodiver.py -d orthodiver_input -A allacma_fusca.1 -B sminthurus_viridis.1 -o orthodiver_output/v1

#should look at distribution of dn/ds values to see if they make sense across genes before calculating further stuff
#hopefully if i do this can at least get dn/ds estimates where i know alignment quality

write a script to process output dn/ds estimates per orthogroup into one file for plotting