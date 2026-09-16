#!/usr/bin/env python3
"""Test whether dN/dS differs between Allacma fusca X1 and X2 genes.

The script joins:
  1. DESeq2 gene-level results (Geneid like g1000)
  2. single-copy orthogroup mappings (allacma_fusca.g1000.t1)
  3. orthogroup-level population genetic summaries (dnds)
  4. syngraph chromosome coordinates for A. fusca orthogroups

It reports X1 vs X2 comparisons for all mapped orthogroups and for the subset
whose A. fusca gene is male-biased by the chosen DESeq2 thresholds.
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Optional


X_CHROMS = {
    "OX359249.1": "X1",
    "OX359250.1": "X2",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare dN/dS between Allacma fusca X1 and X2 overall and in male-biased genes."
    )
    parser.add_argument(
        "--de",
        default="data/results/diff_expr/allacma_fusca/deseq2/v_4/allacma_fusca.DEseq2_results_raw.tsv",
        help="DESeq2 results TSV with Geneid, log2FoldChange, and padj columns.",
    )
    parser.add_argument(
        "--popgen",
        default=(
            "data/results/ortholog_pop_gen/allacma_fusca.vs.sminthurus_viridis/"
            "allacma_fusca.sminthurus_viridis.orthodiver_agg.tsv"
        ),
        help="Orthogroup-level pop-gen TSV containing Orthogroup and dnds columns.",
    )
    parser.add_argument(
        "--orthogroups",
        default="data/results/ortholog_pop_gen/allacma_fusca.vs.sminthurus_viridis/SC_orthogroups.txt",
        help="Single-copy orthogroup mapping file.",
    )
    parser.add_argument(
        "--chrom-table",
        default="Allacma_fusca.tsv",
        help="A. fusca syngraph table: Orthogroup, chromosome, start, end.",
    )
    parser.add_argument(
        "--out-summary",
        default="data/results/ortholog_pop_gen/allacma_fusca.vs.sminthurus_viridis/x1_x2_dnds_tests.tsv",
        help="Output summary TSV.",
    )
    parser.add_argument(
        "--out-joined",
        default="data/results/ortholog_pop_gen/allacma_fusca.vs.sminthurus_viridis/x1_x2_dnds_joined.tsv",
        help="Output joined per-orthogroup TSV.",
    )
    parser.add_argument(
        "--padj-threshold",
        type=float,
        default=0.05,
        help="padj cutoff for male-biased genes. Use 1 to keep all positive-log2FC genes.",
    )
    parser.add_argument(
        "--log2fc-min",
        type=float,
        default=0.0,
        help="Minimum log2FoldChange for male-biased genes.",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=10000,
        help="Permutation count for the two-sided median-difference test.",
    )
    parser.add_argument("--seed", type=int, default=1, help="Random seed for permutations.")
    return parser.parse_args()


def as_float(value: Optional[str]) -> float:
    if value is None:
        return math.nan
    value = value.strip()
    if value == "" or value.upper() == "NA":
        return math.nan
    try:
        return float(value)
    except ValueError:
        return math.nan


def read_tsv_dict(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_de(path: Path, padj_threshold: float, log2fc_min: float) -> dict[str, dict[str, object]]:
    genes: dict[str, dict[str, object]] = {}
    for row in read_tsv_dict(path):
        gene = row["Geneid"]
        lfc = as_float(row.get("log2FoldChange"))
        padj = as_float(row.get("padj"))
        male_biased = math.isfinite(lfc) and lfc > log2fc_min
        if padj_threshold < 1:
            male_biased = male_biased and math.isfinite(padj) and padj <= padj_threshold
        genes[gene] = {
            "log2FoldChange": lfc,
            "padj": padj,
            "male_biased": male_biased,
        }
    return genes


def gene_id_from_member(member: str, species_prefix: str = "allacma_fusca.") -> Optional[str]:
    if not member.startswith(species_prefix):
        return None
    rest = member[len(species_prefix) :]
    match = re.match(r"^(g\d+)(?:\.|$)", rest)
    return match.group(1) if match else None


def read_gene_to_orthogroups(path: Path) -> tuple[dict[str, set[str]], Counter]:
    gene_to_ogs: dict[str, set[str]] = defaultdict(set)
    species_counts: Counter = Counter()
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            og, members_text = line.split(":", 1)
            for member in members_text.split():
                species = member.split(".", 1)[0]
                species_counts[species] += 1
                gene = gene_id_from_member(member)
                if gene:
                    gene_to_ogs[gene].add(og)
    return gene_to_ogs, species_counts


def read_chrom_table(path: Path) -> dict[str, dict[str, str]]:
    chrom_by_og: dict[str, dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 4:
                continue
            og, chrom, start, end = row[:4]
            chrom_by_og[og] = {
                "chrom": chrom,
                "x_chrom": X_CHROMS.get(chrom, ""),
                "start": start,
                "end": end,
            }
    return chrom_by_og


def read_popgen(path: Path) -> dict[str, dict[str, str]]:
    return {row["Orthogroup"]: row for row in read_tsv_dict(path)}


def median(values: list[float]) -> float:
    return statistics.median(values) if values else math.nan


def mean(values: list[float]) -> float:
    return statistics.fmean(values) if values else math.nan


def rankdata(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[indexed[k][0]] = rank
        i = j
    return ranks


def mann_whitney_u(x: list[float], y: list[float]) -> dict[str, float]:
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return {"u": math.nan, "z": math.nan, "p": math.nan}

    combined = x + y
    ranks = rankdata(combined)
    r1 = sum(ranks[:n1])
    u1 = r1 - n1 * (n1 + 1) / 2.0
    u2 = n1 * n2 - u1
    u = min(u1, u2)

    counts = Counter(combined)
    tie_sum = sum(count**3 - count for count in counts.values())
    n = n1 + n2
    variance = n1 * n2 / 12.0 * ((n + 1) - tie_sum / (n * (n - 1))) if n > 1 else math.nan
    if not math.isfinite(variance) or variance <= 0:
        return {"u": u, "z": math.nan, "p": math.nan}

    mean_u = n1 * n2 / 2.0
    correction = 0.5 if u < mean_u else -0.5
    z = (u - mean_u + correction) / math.sqrt(variance)
    p = math.erfc(abs(z) / math.sqrt(2.0))
    return {"u": u, "z": z, "p": p}


def permutation_pvalue(
    x: list[float], y: list[float], permutations: int, rng: random.Random
) -> tuple[float, float]:
    if not x or not y or permutations <= 0:
        return math.nan, math.nan

    observed = median(x) - median(y)
    combined = x + y
    n1 = len(x)
    more_extreme = 0
    for _ in range(permutations):
        rng.shuffle(combined)
        diff = median(combined[:n1]) - median(combined[n1:])
        if abs(diff) >= abs(observed):
            more_extreme += 1
    p = (more_extreme + 1) / (permutations + 1)
    return observed, p


def cliffs_delta(x: list[float], y: list[float]) -> float:
    if not x or not y:
        return math.nan
    greater = 0
    less = 0
    for xv in x:
        for yv in y:
            if xv > yv:
                greater += 1
            elif xv < yv:
                less += 1
    return (greater - less) / (len(x) * len(y))


def fmt(value: object) -> str:
    if isinstance(value, float):
        if math.isnan(value):
            return "NA"
        return f"{value:.8g}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key, "")) for key in fieldnames})


def summarize(label: str, rows: list[dict[str, object]], permutations: int, rng: random.Random) -> dict[str, object]:
    x1 = [row["dnds"] for row in rows if row["x_chrom"] == "X1"]
    x2 = [row["dnds"] for row in rows if row["x_chrom"] == "X2"]
    mw = mann_whitney_u(x1, x2)
    median_diff, perm_p = permutation_pvalue(x1, x2, permutations, rng)
    return {
        "set": label,
        "n_x1": len(x1),
        "n_x2": len(x2),
        "median_x1": median(x1),
        "median_x2": median(x2),
        "median_x1_minus_x2": median_diff,
        "mean_x1": mean(x1),
        "mean_x2": mean(x2),
        "mean_x1_minus_x2": mean(x1) - mean(x2) if x1 and x2 else math.nan,
        "mann_whitney_u": mw["u"],
        "mann_whitney_z": mw["z"],
        "mann_whitney_p_approx": mw["p"],
        "permutation_median_p": perm_p,
        "cliffs_delta_x1_vs_x2": cliffs_delta(x1, x2),
    }


def main() -> None:
    args = parse_args()
    de_path = Path(args.de)
    popgen_path = Path(args.popgen)
    og_path = Path(args.orthogroups)
    chrom_path = Path(args.chrom_table)

    for path in [de_path, popgen_path, og_path, chrom_path]:
        if not path.exists():
            raise SystemExit(f"Missing input file: {path}")

    de = read_de(de_path, args.padj_threshold, args.log2fc_min)
    gene_to_ogs, _ = read_gene_to_orthogroups(og_path)
    popgen = read_popgen(popgen_path)
    chrom_by_og = read_chrom_table(chrom_path)

    og_to_de_genes: dict[str, list[str]] = defaultdict(list)
    for gene in de:
        for og in gene_to_ogs.get(gene, set()):
            og_to_de_genes[og].append(gene)

    joined: list[dict[str, object]] = []
    for og, stats in popgen.items():
        dnds = as_float(stats.get("dnds"))
        chrom_info = chrom_by_og.get(og)
        if not chrom_info or not chrom_info["x_chrom"] or not math.isfinite(dnds):
            continue

        genes = sorted(og_to_de_genes.get(og, []))
        male_biased = any(bool(de[gene]["male_biased"]) for gene in genes)
        best_gene = genes[0] if len(genes) == 1 else ";".join(genes)
        lfc_values = [de[gene]["log2FoldChange"] for gene in genes if math.isfinite(de[gene]["log2FoldChange"])]
        padj_values = [de[gene]["padj"] for gene in genes if math.isfinite(de[gene]["padj"])]

        joined.append(
            {
                "Orthogroup": og,
                "allacma_fusca_gene": best_gene,
                "chrom": chrom_info["chrom"],
                "x_chrom": chrom_info["x_chrom"],
                "start": chrom_info["start"],
                "end": chrom_info["end"],
                "dnds": dnds,
                "male_biased": "TRUE" if male_biased else "FALSE",
                "max_log2FoldChange": max(lfc_values) if lfc_values else math.nan,
                "min_padj": min(padj_values) if padj_values else math.nan,
            }
        )

    rng = random.Random(args.seed)
    summary = [
        summarize("all_x_orthogroups", joined, args.permutations, rng),
        summarize(
            f"male_biased_log2fc_gt_{args.log2fc_min:g}_padj_le_{args.padj_threshold:g}",
            [row for row in joined if row["male_biased"] == "TRUE"],
            args.permutations,
            rng,
        ),
    ]

    write_tsv(
        Path(args.out_joined),
        joined,
        [
            "Orthogroup",
            "allacma_fusca_gene",
            "chrom",
            "x_chrom",
            "start",
            "end",
            "dnds",
            "male_biased",
            "max_log2FoldChange",
            "min_padj",
        ],
    )
    write_tsv(
        Path(args.out_summary),
        summary,
        [
            "set",
            "n_x1",
            "n_x2",
            "median_x1",
            "median_x2",
            "median_x1_minus_x2",
            "mean_x1",
            "mean_x2",
            "mean_x1_minus_x2",
            "mann_whitney_u",
            "mann_whitney_z",
            "mann_whitney_p_approx",
            "permutation_median_p",
            "cliffs_delta_x1_vs_x2",
        ],
    )

    print(f"Wrote joined table: {args.out_joined}")
    print(f"Wrote summary: {args.out_summary}")
    for row in summary:
        print(
            "\t".join(
                fmt(row[key])
                for key in [
                    "set",
                    "n_x1",
                    "n_x2",
                    "median_x1",
                    "median_x2",
                    "median_x1_minus_x2",
                    "mann_whitney_p_approx",
                    "permutation_median_p",
                ]
            )
        )


if __name__ == "__main__":
    main()
