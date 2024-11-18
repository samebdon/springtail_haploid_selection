#BSUB -o logs/fastp.out.%J
#BSUB -e logs/fastp.err.%J
#BSUB -q normal
#BSUB -n 8
#BSUB -M 4096
#BSUB -R "select[mem>4096] rusage[mem=4096]"

fastp -i data/raw_data/reseq/lipothrix_lubbocki/SRR22681197_1.fastq -I data/raw_data/reseq/lipothrix_lubbocki/SRR22681197_2.fastq -o data/raw_data/reseq/lipothrix_lubbocki/SRR22681197_1.fastp.fastq -O data/raw_data/reseq/lipothrix_lubbocki/SRR22681197_1.fastp.fastq --length_required 33 --cut_front --cut_tail --cut_mean_quality 20 --thread 8
