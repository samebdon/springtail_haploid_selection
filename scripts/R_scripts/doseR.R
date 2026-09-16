if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("doseR")

library(doseR)
library(SummarizedExperiment)