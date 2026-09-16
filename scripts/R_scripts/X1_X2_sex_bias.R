library(dplyr)
library(readr)
library(ggplot2)

# Read DESeq2 results
de <- read_tsv(
    "data/results/diff_expr/allacma_fusca/deseq2/v_4/allacma_fusca.DEseq2_results_raw.tsv",
    show_col_types = FALSE
)

# Read GFF3
gff <- read_tsv(
    "data/results/braker3/allacma_fusca/braker.gff3",
    comment = "#",
    col_names = c(
        "chr", "source", "type", "start", "end",
        "score", "strand", "phase", "attributes"
    ),
    show_col_types = FALSE
)

# Extract gene locations
genes <- gff %>%
    filter(type == "gene") %>%
    mutate(
        Geneid = sub("ID=([^;]+).*", "\\1", attributes)
    ) %>%
    select(Geneid, chr)

# Join DE results to chromosome
dat <- left_join(de, genes, by = "Geneid")

# Keep the chromosomes of interest
dat <- dat %>%
    filter(chr %in% c("OX359250.1", "OX359249.1"))

# Classify genes
dat <- dat %>%
    mutate(
        bias = case_when(
            padj < 0.05 & log2FoldChange > 0 ~ "male",
            padj < 0.05 & log2FoldChange < 0 ~ "female",
            TRUE ~ "unbiased"
        )
    )

# Counts
tab <- table(dat$chr, dat$bias)
print(tab)

# Proportions within chromosome
prop.table(tab, margin = 1)

# Statistical test
chisq.test(tab)

# If expected counts are small, also run:
fisher.test(tab)

# Plot
plot_dat <- dat %>%
    count(chr, bias) %>%
    group_by(chr) %>%
    mutate(prop = n / sum(n))

p <- ggplot(plot_dat,
            aes(chr, prop, fill = bias)) +
    geom_col() +
    labs(
        x = "Chromosome",
        y = "Proportion of genes",
        fill = "Expression bias"
    ) +
    theme_classic()

ggsave(
    "data/results/diff_expr/allacma_fusca/deseq2/v_4/sex_bias_proportions_by_chromosome.pdf",
    p,
    width = 5,
    height = 4
)

ggsave(
    "data/results/diff_expr/allacma_fusca/deseq2/v_4/sex_bias_proportions_by_chromosome.png",
    p,
    width = 5,
    height = 4,
    dpi = 300
)
