"""describe_gtf.py

Usage:
describe_gtf.py -f <FILE> [-h]

Options:
-f, --file <FILE>	       GTF file

"""

# Example Command
# python describe_gtf.py -f braker.gtf

import pandas as pd
import numpy as np
from docopt import docopt


def get_gene_from_exon(attribute):
    return attribute.split(" ")[-1].split('"')[1]


def get_transcript_from_exon(attribute):
    return attribute.split(" ")[1].split('"')[1]


def get_gene_from_transcript(attribute):
    return attribute.split(".")[0]


def superfeature_array_to_count_dist(arr):
    feature, feature_counts = np.unique(arr, return_counts=True)
    i_arr, n_arr = np.unique(feature_counts, return_counts=True)
    distr_arr = np.zeros(i_arr[-1], dtype=np.uint32)
    for i, n in zip(i_arr, n_arr):
        distr_arr[i - 1] = n
    return distr_arr


def get_attribute_count_distribution(gtf_df, superfeature="gene", subfeature="exon"):
    attributes = df[df["feature"] == subfeature]["attribute"]

    if subfeature == "exon":
        if superfeature == "gene":
            arr = np.array([get_gene_from_exon(attribute) for attribute in attributes])
        if superfeature == "transcript":
            arr = np.array(
                [get_transcript_from_exon(attribute) for attribute in attributes]
            )
    if subfeature == "transcript":
        if superfeature == "gene":
            arr = np.array(
                [get_gene_from_transcript(attribute) for attribute in attributes]
            )

    return superfeature_array_to_count_dist(arr)


def find_median(count_arr):
    n = 0
    for i, count in enumerate(count_arr):
        n += count
        if n >= (np.sum(count_arr) + 1) / 2:
            return i + 1
            break


if __name__ == "__main__":
    args = docopt(__doc__)

    columns = [
        "sequence",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]

    df = pd.read_csv(args["--file"], sep="\t", names=columns, index_col=None)
    df["length"] = df["end"] - df["start"]
    gene_exon_count_distribution = get_attribute_count_distribution(df, "gene", "exon")
    transcript_exon_count_distribution = get_attribute_count_distribution(
        df, "transcript", "exon"
    )
    gene_transcript_count_distribution = get_attribute_count_distribution(
        df, "gene", "transcript"
    )

    n_genes = df["feature"].value_counts()["gene"]
    median_gene_length = int(np.median(df[df["feature"] == "gene"]["length"]))
    max_gene_length = int(np.max(df[df["feature"] == "gene"]["length"]))
    min_gene_length = int(np.min(df[df["feature"] == "gene"]["length"]))
    n_short_genes = df[(df["feature"] == "gene") & (df["length"] < 200)].shape[0]
    n_transcripts = df["feature"].value_counts()["transcript"]
    median_transcript_length = int(
        np.median(df[df["feature"] == "transcript"]["length"])
    )

    n_exons = df["feature"].value_counts()["exon"]
    median_exon_length = int(np.median(df[df["feature"] == "exon"]["length"]))
    n_mono_exonic_genes = gene_exon_count_distribution[0]
    n_multi_exonic_genes = np.sum(gene_exon_count_distribution) - n_mono_exonic_genes
    median_exons_per_transcript = find_median(transcript_exon_count_distribution)
    max_exons_per_transcript = len(transcript_exon_count_distribution)
    median_exons_per_gene = find_median(gene_exon_count_distribution)
    max_exons_per_gene = len(gene_exon_count_distribution)
    median_transcripts_per_gene = find_median(gene_transcript_count_distribution)
    max_transcripts_per_gene = len(gene_transcript_count_distribution)

    print(f"Number of genes: {n_genes}")
    print(f"Median gene length: {median_gene_length}")
    print(f"Max gene length: {max_gene_length}")
    print(f"Min gene length: {min_gene_length}")
    print(f"Number of genes < 200bp: {n_short_genes}")
    print(f"Number of transcripts: {n_transcripts}")
    print(f"Median transcript length: {median_transcript_length}")
    print(f"Median number of transcripts per gene: {median_transcripts_per_gene}")
    print(f"Max number of transcripts per gene: {max_transcripts_per_gene}")
    print(f"Number of exons: {n_exons}")
    print(f"Median exon length: {median_exon_length}")
    print(f"Median number of exons per gene: {median_exons_per_gene}")
    print(f"Max number of exons per gene: {max_exons_per_gene}")
    print(f"Median number of exons per transcript: {median_exons_per_transcript}")
    print(f"Maximum number of exons per transcript: {max_exons_per_transcript}")
    print(
        f"Ratio of monoexonic to multi-exonic genes: {np.round(n_mono_exonic_genes/n_multi_exonic_genes, 5)}"
    )

