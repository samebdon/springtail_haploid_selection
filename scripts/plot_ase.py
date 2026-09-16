import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

X_CHROMS = {"OX359249.1","OX359250.1"}
MIN_DP = 8

ase = pd.read_csv(
    "ase_raw.tsv",
    sep="\t",
    names=["chrom", "pos", "sample", "gt", "ad"]
)
meta = pd.read_csv("samples.tsv", sep="\t")


ase = ase[ase["gt"] == "0/1"]
ase[["ref", "alt"]] = ase["ad"].str.split(",", expand=True).astype(int)
ase["dp"] = ase["ref"] + ase["alt"]
ase = ase[ase["dp"] >= MIN_DP]

ase["alt_prop"] = ase["alt"] / ase["dp"]

ase["chrom_class"] = ase["chrom"].apply(
    lambda c: "X" if c in X_CHROMS else "Autosome"
)
ase = ase.merge(meta, on="sample", how="left")

sample_means = (
    ase
    .groupby(["sample", "sex", "chrom_class"], as_index=False)
    .agg(mean_alt_prop=("alt_prop", "mean"))
)

sns.set(style="whitegrid")

plt.figure(figsize=(8, 5))

ax = sns.violinplot(
    data=sample_means,
    x="chrom_class",
    y="mean_alt_prop",
    hue="sex",
    cut=0,
    inner=None,
    alpha=0.6
)

sns.stripplot(
    data=sample_means,
    x="chrom_class",
    y="mean_alt_prop",
    hue="sex",
    dodge=True,
    jitter=True,
    size=6,
    linewidth=0.5,
    edgecolor="black"
)

plt.axhline(0.5, linestyle="--", color="black")

plt.ylim(0, 1)
plt.xlabel("")
plt.ylabel("Mean ALT allele proportion")
plt.title("Biallelic expression across chromosomes")

handles, labels = ax.get_legend_handles_labels()
plt.legend(handles[:2], labels[:2], title="Sex")

plt.tight_layout()
plt.savefig("ase_alt_proportion_violin.png", dpi=300)
plt.close()
