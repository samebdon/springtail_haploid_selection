library(dplyr)
library(readr)
library(ggplot2)

## Files
tpm_file <- "data/results/diff_expr/allacma_fusca/deseq2/v_4/renamed_DESeq2_TPM_matrix.tsv"
gff_file <- "data/results/braker3/allacma_fusca/braker.gff3"
outdir <- "data/results/diff_expr/allacma_fusca/deseq2/v_4"

x_chrs <- c("OX359249.1", "OX359250.1")

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
dat_all <- left_join(tpm, genes, by = "Geneid") |>
    filter(!is.na(chr))

## Sample columns
female_cols <- grep("^AF_F_", names(dat_all), value = TRUE)
male_cols   <- grep("^AF_M_", names(dat_all), value = TRUE)

## Mean TPM per gene by sex
dat_all <- dat_all |>
    rowwise() |>
    mutate(
        female_mean = mean(c_across(all_of(female_cols))),
        male_mean   = mean(c_across(all_of(male_cols)))
    ) |>
    ungroup()

## X/autosome ratios by sex
autosome_means <- dat_all |>
    filter(!chr %in% x_chrs) |>
    summarise(
        female_autosome_mean = mean(female_mean),
        male_autosome_mean   = mean(male_mean),
        autosome_genes       = n()
    )

x_autosome_ratios <- dat_all |>
    filter(chr %in% x_chrs) |>
    group_by(chr) |>
    summarise(
        female_x_mean = mean(female_mean),
        male_x_mean   = mean(male_mean),
        x_genes       = n(),
        .groups = "drop"
    ) |>
    mutate(
        female_autosome_mean = autosome_means$female_autosome_mean,
        male_autosome_mean   = autosome_means$male_autosome_mean,
        autosome_genes       = autosome_means$autosome_genes,
        female_x_autosome_ratio = female_x_mean / female_autosome_mean,
        male_x_autosome_ratio   = male_x_mean / male_autosome_mean
    )

combined_x_autosome_ratios <- dat_all |>
    filter(chr %in% x_chrs) |>
    summarise(
        chr = "combined_X",
        female_x_mean = mean(female_mean),
        male_x_mean   = mean(male_mean),
        x_genes       = n()
    ) |>
    mutate(
        female_autosome_mean = autosome_means$female_autosome_mean,
        male_autosome_mean   = autosome_means$male_autosome_mean,
        autosome_genes       = autosome_means$autosome_genes,
        female_x_autosome_ratio = female_x_mean / female_autosome_mean,
        male_x_autosome_ratio   = male_x_mean / male_autosome_mean
    )

x_autosome_ratios <- bind_rows(
    x_autosome_ratios,
    combined_x_autosome_ratios
)

write_tsv(
    x_autosome_ratios,
    file.path(outdir, "X_autosome_expression_ratios.tsv")
)

## Keep the two X chromosomes and combine them
dat <- dat_all |>
    filter(chr %in% x_chrs) |>
    mutate(chr = "X")

## Keep genes expressed in both sexes
dat <- dat |>
    filter(female_mean > 0, male_mean > 0)

## Remove lower and upper 5% of expression
f_lim <- quantile(dat$female_mean, c(0.05, 0.95))
m_lim <- quantile(dat$male_mean, c(0.05, 0.95))

dat <- dat |>
    filter(
        between(female_mean, f_lim[1], f_lim[2]),
        between(male_mean,   m_lim[1], m_lim[2])
    )

## Calculate log2(M/F)
dat <- dat |>
    mutate(
        log2MF = log2((male_mean + 0.01) / (female_mean + 0.01))
    )

## Save filtered data
write_tsv(
    dat,
    file.path(outdir, "combined_X_log2MF.tsv")
)

## Statistics
sink(file.path(outdir, "combined_X_log2MF_tests.txt"))

cat("Genes retained\n")
print(nrow(dat))

cat("\nMedian log2(M/F)\n")
print(median(dat$log2MF))

cat("\nSummary\n")
print(summary(dat$log2MF))

cat("\nX/autosome expression ratios\n")
print(x_autosome_ratios)

cat("\nWilcoxon signed-rank test against 0\n\n")
print(wilcox.test(dat$log2MF, mu = 0))

sink()

## Plot
p <- ggplot(dat, aes(x = "X", y = log2MF)) +
    geom_violin(fill = "grey80", trim = FALSE) +
    geom_boxplot(width = 0.15, outlier.shape = NA) +
    geom_jitter(width = 0.08, alpha = 0.35, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
        x = "",
        y = expression(log[2](Male/Female~TPM))
    ) +
    theme_classic()

ggsave(
    file.path(outdir, "combined_X_log2MF.pdf"),
    p,
    width = 3.5,
    height = 4
)

ggsave(
    file.path(outdir, "combined_X_log2MF.png"),
    p,
    width = 3.5,
    height = 4,
    dpi = 300
)
