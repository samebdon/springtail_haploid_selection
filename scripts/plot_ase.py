import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import joypy

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

#density plot
g = sns.FacetGrid(
    sample_means,
    col="sex",
    hue="chrom_class",
    height=4,
    aspect=1.2,
    sharex=True,
    sharey=True
)

g.map(
    sns.kdeplot,
    "mean_alt_prop",
    fill=True,
    alpha=0.4,
    linewidth=1.5
)

# Reference line
for ax in g.axes.flat:
    ax.axvline(0.5, linestyle="--", color="black")

g.set_axis_labels("Mean ALT allele proportion", "Density")
g.add_legend(title="Chromosome")

plt.tight_layout()
plt.savefig("ase_alt_proportion_density.png", dpi=300)
plt.close()

# ridgeplot
fig, axes = joypy.joyplot(
    sample_means,
    by="chrom_class",
    column="mean_alt_prop",
    hue="sex",
    overlap=0.6,
    figsize=(8, 5),
    alpha=0.6,
    legend=True,
    linewidth=1
)

# Reference line
for ax in axes:
    ax.axvline(0.5, linestyle="--", color="black", linewidth=1)

axes[-1].set_xlabel("Mean ALT allele proportion")

plt.tight_layout()
plt.savefig("ase_alt_proportion_ridge.png", dpi=300)
plt.close()