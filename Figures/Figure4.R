setwd("LiverSpatialCompare")

## Genes with zero expression in either dataset

load("RDATA/dataReady_bothData_genesMapped.RData")

geneSums.h <- rowSums(layerMeans)
zeroGenes <- names(which(geneSums.h == 0))
h.gene.means <- rowMeans(halpernUMI)
h.gene.nozeroMean <- rowSums(halpernUMI!=0)


load("RDATA/analysis_WaveCrest_Morten.RData")

data.norm.order <- data.norm[,wc.order]


# Genes that have significant (<.1) zonation in Morten data and expressed in more than 10 cells.

sig.g.morten <- names(which(adjusted.pvals < .1))
sigMorten_zeroHalpern <- intersect(sig.g.morten, zeroGenes)
detect.gene.morten <- apply(data.norm[useMorten.genes,], 1, function(x) sum(x!=0))
sigMorten_zeroHalpern <- intersect(sigMorten_zeroHalpern, names(which(detect.gene.morten >= 10)))
length(sigMorten_zeroHalpern)

sigMorten_zeroHalpern <- subset(geneSet, MortenGene %in% sigMorten_zeroHalpern)

write.table(sigMorten_zeroHalpern, file="OUT/SIGinMorten_zeroInHalpern_Genes_v2.txt", quote=F, row.names=F, col.names=F)

plotGenes <- names(rev(sort(rowMeans(data.norm[sigMorten_zeroHalpern[,1],])))[1:6])
pdf("PLOTS/plotSIGinMorten_zeroInHalpern_scatter_6genes_Fig4.pdf", height=6, width=15)
par(mfrow=c(2,3), mar=c(5,5,2,1))
for(i in 1:6) {

  gene.m <- plotGenes[i]
  currentY.m <- data.norm.order[gene.m,]
  
  getMeans.m <- tapply(currentY.m, layerSplit, mean)
  rescaleY.morten.means <- ((1 - 0)/(max(getMeans.m) - min(getMeans.m)))*(getMeans.m - min(getMeans.m)) + 0

  plot(1:length(currentY.m), log2(currentY.m+1), col="#fc8d59", pch=20, 
       cex.axis=2.2, cex.lab=2.5,bty = 'n', cex.main=2.5,
       cex=1.4, main=gene.m, ylab="Log2 Expression", 
       xlab="Cells")
  
  
  FIT = smooth.spline(1:length(currentY.m), log2(currentY.m+1), df=4)
    lines(FIT$x, FIT$y, lwd=3, col="#fc8d59")
    
}
dev.off()




## How many in the other direction?


sig.g.h <- names(which(layerStatsPvalue[,2] < .01))
zeroGenes <- names(which(rowSums(data.norm.order) == 0))
sigHalpern_zeroM <- intersect(sig.g.h, zeroGenes)
detect.gene.h <- apply(halpernUMI, 1, function(x) sum(x!=0))
sigMorten_zeroHalpern <- intersect(sigHalpern_zeroM, names(which(detect.gene.h >= 10)))


length(sigMorten_zeroHalpern)

## Only 3 such genes.

