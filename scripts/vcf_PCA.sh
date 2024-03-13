mamba install -c bioconda bioconductor-snprelate

nano pops.txt #in order of sample.id below
WW
WW
WW
WW
WW
WW
WW
WW
WW
WW
BH
WW


R 
library("SNPRelate")

version = 1
dir_base = 'figures/PCA/allacma_fusca/v_'
fig_dir = paste(dir_base, version, sep = '')
dir.create(file.path(fig_dir), recursive = TRUE)


vcf.fn <- 'allacma_fusca.hard_filtered.sorted.vcf'
snpgdsVCF2GDS(vcf.fn, "ccm.gds",  method="biallelic.only")
genofile <- snpgdsOpen("ccm.gds")
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pop_code <- scan("pops.txt", what=character())
cbind( sample.id, pop_code)
ccm_pca<-snpgdsPCA(genofile, autosome.only=FALSE, num.thread=4)
tab <- data.frame(sample.id = ccm_pca$sample.id,
              pop = factor(pop_code)[match(ccm_pca$sample.id, sample.id)],
              EV1 = ccm_pca$eigenvect[,1],    # the first eigenvector
              EV2 = ccm_pca$eigenvect[,2],    # the second eigenvector
              stringsAsFactors = FALSE)

png(filename=paste(fig_dir, 'PCA.png', sep = '/'))
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), 
	xlab=paste("eigenvector 2 (",round(ccm_pca$varprop[2], 2),")"), 
	ylab=paste("eigenvector 1 (",round(ccm_pca$varprop[1],2),")"))
legend("topright", legend=levels(tab$pop), pch="o", col=1:nlevels(tab$pop))
tab$EV2_jitter <- jitter(tab$EV2, 30)
text(tab$EV2_jitter, tab$EV1+0.025, labels=tab$sample.id, cex = 0.5)
dev.off()
