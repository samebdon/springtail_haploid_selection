library(dplyr)
library(readr)
library(ggplot2)

## Files
tpm_file <- "data/results/diff_expr/allacma_fusca/deseq2/v_4/renamed_DESeq2_TPM_matrix.tsv"
gff_file <- "data/results/braker3/allacma_fusca/braker.gff3"
outdir <- "data/results/diff_expr/allacma_fusca/deseq2/v_4"

## Read TPM matrix
tpm <- read_tsv(
    tpm_file,
    show_col_types = FALSE
) |>
    select(-rowname)

## Read GFF
gff <- read_tsv(
    gff_file,
    comment = "#",
    col_names = c(
        "chr","source","type","start","end",
        "score","strand","phase","attributes"
    ),
    show_col_types = FALSE
)

genes <- gff |>
    filter(type == "gene") |>
    mutate(
        Geneid = sub("ID=([^;]+).*", "\\1", attributes)
    ) |>
    select(Geneid, chr)

## Join chromosome information
dat <- left_join(tpm, genes, by = "Geneid")

## Keep chromosomes of interest
dat <- dat |>
    filter(chr %in% c("OX359249.1", "OX359250.1"))

## Sample columns
female_cols <- grep("^AF_F_", names(dat), value = TRUE)
male_cols   <- grep("^AF_M_", names(dat), value = TRUE)

## Mean TPM per sex
dat <- dat |>
    rowwise() |>
    mutate(
        female_mean = mean(c_across(all_of(female_cols))),
        male_mean   = mean(c_across(all_of(male_cols)))
    ) |>
    ungroup()

## Expressed in both sexes
dat <- dat |>
    filter(female_mean > 0, male_mean > 0)

## Remove bottom/top 5%
f_lim <- quantile(dat$female_mean, c(0.05, 0.95))
m_lim <- quantile(dat$male_mean, c(0.05, 0.95))

dat <- dat |>
    filter(
        between(female_mean, f_lim[1], f_lim[2]),
        between(male_mean,   m_lim[1], m_lim[2])
    )

## log2(M/F)
dat <- dat |>
    mutate(
        log2MF = log2((male_mean + 0.01)/(female_mean + 0.01))
    )

## Save filtered table
write_tsv(
    dat,
    file.path(outdir, "X_chromosome_log2MF.tsv")
)

## Statistics
sink(file.path(outdir, "X_chromosome_log2MF_tests.txt"))

cat("Genes retained\n")
print(table(dat$chr))

cat("\nMedian log2(M/F)\n")
print(tapply(dat$log2MF, dat$chr, median))

cat("\nWilcoxon signed-rank tests against 0\n\n")

for(chr in unique(dat$chr)) {
    cat(chr, "\n")
    print(wilcox.test(
        dat$log2MF[dat$chr == chr],
        mu = 0
    ))
    cat("\n")
}

cat("\nWilcoxon rank-sum test between chromosomes\n\n")
print(wilcox.test(log2MF ~ chr, data = dat))

sink()

## Plot
p <- ggplot(dat, aes(chr, log2MF, fill = chr)) +
    geom_violin(trim = FALSE, alpha = 0.5) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    geom_jitter(width = 0.08, alpha = 0.4, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        x = "",
        y = expression(log[2](Male/Female~TPM))
    ) +
    theme_classic() +
    theme(legend.position = "none")

ggsave(
    file.path(outdir, "X_chromosome_log2MF.pdf"),
    p,
    width = 5,
    height = 4
)

ggsave(
    file.path(outdir, "X_chromosome_log2MF.png"),
    p,
    width = 5,
    height = 4,
    dpi = 300
)
