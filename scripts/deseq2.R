library("tibble")
library("Rsubread")
library("DESeq2")
library( "gplots" )
library( "RColorBrewer" )
library( "genefilter" )

version = 2
dir_base = 'figures/deseq2/allacma_fusca/v_'
fig_dir = paste(dir_base, version, sep = '')
results_file_base = 'data/results/diff_expr/allacma_fusca/deseq2/v_'
results_dir = paste(results_file_base, version, sep = '')

dir.create(file.path(fig_dir), recursive = TRUE)
dir.create(file.path(results_dir), recursive = TRUE)

featureCounts_file = 'data/results/diff_expr/allacma_fusca/featureCounts/allacma_fusca.braker3.featureCounts.txt' 
sampleInfo_file = 'data/results/diff_expr/allacma_fusca/sampleinfo.tsv'

featureCounts_table <- read.table(featureCounts_file, sep = '\t', header = 1, row.names='Geneid')
sampleInfo <- read.table(sampleInfo_file, sep='\t', header=1)

colnames(featureCounts_table) <- c('Chr','Start','End','Strand','Length',
    'AF_F_1','AF_F_2','AF_F_3','AF_F_4','AF_F_5','AF_F_6','AF_F_7','AF_F_8','AF_F_9','AF_F_10',
    'AF_M_1','AF_M_2','AF_M_3','AF_M_4','AF_M_5','AF_M_6','AF_M_7','AF_M_8','AF_M_9','AF_M_10')

counts <- featureCounts_table[, 6:25]
#counts <- featureCounts(bams, blablabla_restofcommand)$
#do we want a more complicated design, modelling the number of samples per dataset or sampling location?
dds <- DESeqDataSetFromMatrix(countData=counts, colData=sampleInfo, design=~sex)

dds <- DESeq(dds) #, betaPrior = betaPrior)
resultsNames(dds) # lists the coefficients

write.table(as.data.frame(counts(dds)) %>% rownames_to_column('Geneid'), 
    file=paste(results_dir, 'allacma_fusca.DEseq2_counts.tsv', sep = '/'), 
    sep = '\t', 
    row.names = FALSE, 
    quote = FALSE)

res <- results(dds, name="sex_male_vs_female")
resFilt <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
# test: sex male vs female
# 8326
# 5816 positive
# 2510 negative
# whats the polarity of male and female bias? positive male biased negative female biased?

write.table(as.data.frame(res) %>% rownames_to_column('Geneid'), 
    file=paste(results_dir, 'allacma_fusca.DEseq2_results_raw.tsv', sep = '/'), 
    sep = '\t', 
    row.names = FALSE, 
    quote = FALSE)

write.table(as.data.frame(resFilt) %>% rownames_to_column('Geneid'), 
    file=paste(results_dir, 'allacma_fusca.DEseq2_results_filtered.tsv', sep = '/'), 
    sep = '\t', 
    row.names = FALSE, 
    quote = FALSE)

#make a geneid bed file

png(filename=paste(fig_dir, 'MA.png', sep = '/'))
plotMA( res, ylim = c(min(na.omit(res$log2FoldChange)), max(na.omit(res$log2FoldChange))) )
dev.off()

png(filename=paste(fig_dir, 'dispersion_estimates.png', sep = '/'))
plotDispEsts( dds, ylim = c(1e-3, 1e3) )
dev.off()

png(filename=paste(fig_dir, 'p_hist.png', sep = '/'))
hist( res$pvalue, breaks=100, col="grey" )
dev.off()

qs <- c( 0, quantile( res$baseMean[res$baseMean > 0], 0:7/7 ) )
bins <- cut( res$baseMean, qs )
levels(bins) <- paste0("~",round(.5*qs[-1] + .5*qs[-length(qs)]))
ratios <- tapply( res$pvalue, bins, function(p) mean( p < .01, na.rm=TRUE ) ) # plot these ratios
png(filename=paste(fig_dir, 'ratio_of_small_ps.png', sep = '/'))
barplot(ratios, xlab="mean normalized count", ylab="ratio of small $p$ values")
dev.off()

metadata(res)$filterThreshold
##  1.375557% 
## 0.02148373 

png(filename=paste(fig_dir, 'filter_rej_dist.png', sep = '/'))
plot(metadata(res)$filterNumRej,type="b", xlab="quantiles of 'baseMean'", ylab="number of rejections")
dev.off()

# rlog QC
rld <- rlogTransformation(dds)

png(filename=paste(fig_dir, 'rlog_before_after.png', sep = '/'))
par( mfrow = c( 1, 2 ) )
plot( log2( 1+counts(dds, normalized=TRUE)[, 1:2] ), col="#00000020", pch=20, cex=0.3 )
plot( assay(rld)[, 1:2], col="#00000020", pch=20, cex=0.3 )
dev.off()

sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld$id
colnames(sampleDistMatrix) <- rld$id
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
png(filename=paste(fig_dir, 'sample_distances.png', sep = '/'))
heatmap.2( sampleDistMatrix, trace="none", col=colours)
dev.off()

png(filename=paste(fig_dir, 'PCA_sex.png', sep = '/'))
plotPCA(rld, intgroup = c('sex'))
dev.off()

png(filename=paste(fig_dir, 'PCA_ID.png', sep = '/'))
pca <- plotPCA(rld, intgroup = c('id'))
dev.off()
#AF_M_3/4/5 are split byy PC2, maybe could compare PC1 with PC3 to see any female variance?

png(filename=paste(fig_dir, 'top_var_genes.png', sep = '/'))
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row",
     trace="none", dendrogram="column",
     col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))
dev.off()

# or to shrink log fold changes association with condition:
#res <- lfcShrink(dds, coef="sex_male_vs_female", type="apeglm")
