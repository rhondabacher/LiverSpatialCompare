setwd("LiverSpatialCompare")

# Supp. Figure 1 

load("RDATA/dataReady_bothData_genesMapped.RData")

detect.gene.halpern <- apply(halpernUMI.norm[useHalpern.genes,], 1, function(x) mean(x!=0))
detect.gene.morten <- apply(data.norm[useMorten.genes,], 1, function(x) mean(x!=0))

mean.noz.halpern <- apply(halpernUMI.norm[useHalpern.genes,], 1, function(x) mean(x[x!=0]))
mean.noz.morten <- apply(data.norm[useMorten.genes,], 1, function(x) mean(x[x!=0]))

diff_logMean <- (log2((mean.noz.morten[geneSet[,1]])) - log(mean.noz.halpern[geneSet[,2]]))
diff_detect <- (detect.gene.morten[geneSet[,1]] - detect.gene.halpern[geneSet[,2]]) 


## Check GC content and Gene Length differences in the subset of genes in left side of plot.

library(biomaRt)

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
get_gene_props <- getBM(c("mgi_symbol","percentage_gene_gc_content","transcript_length"), mart= mouse)
head(get_gene_props)

diffDetect.neg <- names(which(diff_detect < 0)); length(diffDetect.neg)
diffDetect.pos <- names(which(diff_detect >= 0)); length(diffDetect.pos)

length(diffDetect.neg) / length(diff_detect)

diffMean.neg <- names(which(diff_logMean < 0)); length(diffMean.neg)
diffMean.pos <- names(which(diff_logMean >= 0)); length(diffMean.pos)

length(diffMean.neg) / length(diff_logMean)

diffDetect.neg <- unique(subset(get_gene_props, mgi_symbol %in% diffDetect.neg))
diffDetect.pos <- unique(subset(get_gene_props, mgi_symbol %in% diffDetect.pos))

diffMean.neg <- unique(subset(get_gene_props, mgi_symbol %in% diffMean.neg))
diffMean.pos <- unique(subset(get_gene_props, mgi_symbol %in% diffMean.pos))


pdf("PLOTS/detectionRat_vs_Expr_compareDatasets_on_GC_Fig1Supp.pdf", height=3, width=7)
par(mar=c(5,5,1,1), mgp=c(3,1,0), mfrow=c(1,2))
plot(density(diffDetect.neg$percentage_gene_gc_content), ylim=c(0,.1),
      lwd=1.5, xlab= "GC Content", ylab= "Density", main="", xlim=c(20,80),
      cex.lab=1.5, cex.axis=2, pch=20, col="cadetblue3", bty='n', xaxt='n', yaxt='n')
axis(1, at=seq(20,80, by=10), cex.axis=1.3)
axis(2, at=seq(0,.1, by=.05), cex.axis=1.3)
lines(density(diffDetect.pos$percentage_gene_gc_content), col="gray40", lwd=2)
lines(density(diffMean.neg$percentage_gene_gc_content), col="cadetblue3", lwd=2)
lines(density(diffMean.pos$percentage_gene_gc_content), col="gray40", lwd=2)
legend(20,.099, title="Genes with higher \n detection or mean in: ", 
            c("Full-length Data", "UMI Data"), bty='n', xjust=0,
            lwd=2, col=c("gray40", "cadetblue3"), cex=.5)

plot(density(log2(diffDetect.neg$transcript_length)), ylim=c(0,.6),
      lwd=1.5, xlab= "log2 Gene Length", ylab= "Density", main="", xlim=c(3, 12),
      cex.lab=1.5, cex.axis=2, pch=20, col="cadetblue3", bty='n', xaxt='n', yaxt='n')
axis(1, at=seq(3,12, by=2), cex.axis=1.3)
axis(2, at=seq(0,.6, by=.2), cex.axis=1.3)
lines(density(log2(diffDetect.pos$transcript_length)), col="gray40", lwd=2)
lines(density(log2(diffMean.neg$transcript_length)), col="cadetblue3", lwd=2)
lines(density(log2(diffMean.pos$transcript_length)), col="gray40", lwd=2)
legend(7,.58, title="Genes with higher \n detection or mean in: ", 
            c("Full-length Data", "UMI Data"), bty='n', xjust=0,
            lwd=2, col=c("gray40", "cadetblue3"), cex=.5)
dev.off()