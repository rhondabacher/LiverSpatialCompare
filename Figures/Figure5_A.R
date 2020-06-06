setwd("LiverSpatialCompare")

## Genes with zero expression in either dataset

load("RDATA/dataReady_bothData_genesMapped.RData")



load("RDATA/analysis_WaveCrest_Morten.RData")

data.norm.order <- data.norm[,wc.order]

load("RDATA/analysis_Monocle_Droplet.RData")


# Genes that have significant (<.1) zonation in Morten data and expressed in more than 10 cells.
geneSums.h <- rowSums(layerMeans)
zeroGenes <- names(which(geneSums.h == 0))
h.gene.means <- rowMeans(halpernUMI)
h.gene.nozeroMean <- rowSums(halpernUMI!=0)


sig.g.morten <- names(which(adjusted.pvals < .1))
sigMorten_zeroHalpern <- intersect(sig.g.morten, subset(geneSet, HalpernGene %in% zeroGenes)[,1])
length(sigMorten_zeroHalpern)
detect.gene.morten <- apply(data.norm[useMorten.genes,], 1, function(x) sum(x!=0))
sigMorten_zeroHalpern <- intersect(sigMorten_zeroHalpern, names(which(detect.gene.morten >= 10)))
length(sigMorten_zeroHalpern)


sigMorten_zeroHalpern <- subset(geneSet, Gene %in% sigMorten_zeroHalpern)


sig.g.morten <- names(which(adjusted.pvals < .1))
zeroGenes <- names(which(rowSums(droplet.norm) == 0))
sigMorten_zeroDroplet <- intersect(sig.g.morten, zeroGenes)
length(sigMorten_zeroDroplet)
detect.gene.morten <- apply(data.norm[useMorten.genes,], 1, function(x) sum(x!=0))
sigMorten_zeroDroplet <- intersect(sigMorten_zeroDroplet, names(which(detect.gene.morten >= 10)))
length(sigMorten_zeroDroplet)

sigMorten_zeroDroplet <- subset(geneSet, Gene %in% sigMorten_zeroDroplet)


tosave <- list(NotDetected_MARSseq = sigMorten_zeroHalpern[,1], NotDetected_10X = sigMorten_zeroDroplet[,1])

for (i in c(1:2)){
  write.xlsx(tosave[i], file="OUT/SIGinMorten_zeroInHalpern_Genes_v2.xlsx", sheetName=paste(i), append=T)
} 

useg <- intersect(sigMorten_zeroHalpern[,1], sigMorten_zeroDroplet[,1])

plotGenes <- names(rev(sort(rowMeans(data.norm[useg,])))[1:6])

# Calculate means for Morten's data for plotting:
layerSplit <- split(1:ncol(data.norm.order), cut(seq_along(1:ncol(data.norm.order)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))

layerMeans.morten <- t(apply(data.norm.order, 1, function(x) {
  return(tapply(x, layerSplit, median))
}))

pdf("PLOTS/plotSIGinMorten_zeroInHalpern_scatter_6genes_Fig4.pdf", height=4, width=12, useDingbats=FALSE)
par(mfrow=c(2,3), mar=c(5,5,2,1))
for(i in 1:6) {

  gene.m <- plotGenes[i]
  currentY.m <- data.norm.order[gene.m,]
  
  getMeans.m <- tapply(currentY.m, layerSplit, mean)
  rescaleY.morten.means <- ((1 - 0)/(max(getMeans.m) - min(getMeans.m)))*(getMeans.m - min(getMeans.m)) + 0

  plot(1:length(currentY.m), log2(currentY.m+1), col="#fc8d59", pch=20, 
       cex.axis=1.2, cex.lab=1.3,bty = 'n', cex.main=2,
       cex=1.4, main=gene.m, ylab="Log2 Expression", 
       xlab="Cells")
  
  
  FIT = smooth.spline(1:length(currentY.m), log2(currentY.m+1), df=4)
    lines(FIT$x, FIT$y, lwd=3, col="#fc8d59")
    
}
dev.off()




## How many in the other direction? (Halpern)
sig.g.h <- names(which(layerStatsPvalue[,2] < .1))
zeroGenes <- names(which(rowSums(data.norm.order) == 0))
sigHalpern_zeroM <- intersect(subset(geneSet, HalpernGene %in% sig.g.h)[,1], zeroGenes)
length(sigHalpern_zeroM)
detect.gene.h <- apply(halpernUMI[useHalpern.genes,], 1, function(x) sum(x!=0))
sigMorten_zeroHalpern <- intersect(sigHalpern_zeroM, names(which(detect.gene.h >= 10)))

length(sigMorten_zeroHalpern)

## Only 10 such genes.



## How many in the other direction? (Droplet)




sig.g.h <- rownames(subset(de.droplet, qval < .1))
zeroGenes <- names(which(rowSums(data.norm.order) == 0))
sigHalpern_zeroM <- intersect(sig.g.h, zeroGenes)
detect.gene.h <- apply(droplet.norm[useDroplet.genes,], 1, function(x) sum(x!=0))
sigMorten_zeroHalpern <- intersect(sigHalpern_zeroM, names(which(detect.gene.h >= 10)))

length(sigMorten_zeroHalpern)

## Only 0 such genes.





sig.g.morten <- names(which(layerStatsPvalue[,2] < .1))
zeroGenes <- names(which(rowSums(droplet.norm) == 0))
sigMorten_zeroHalpern <- intersect(sig.g.morten, zeroGenes)
length(sigMorten_zeroHalpern)
detect.gene.morten <- apply(halpernUMI[useHalpern.genes,], 1, function(x) sum(x!=0))
sigMorten_zeroHalpern <- intersect(sigMorten_zeroHalpern, names(which(detect.gene.morten >= 10)))
length(sigMorten_zeroHalpern)


sig.g.morten <- rownames(subset(de.droplet, qval < .1))
zeroGenes <- names(which(rowSums(halpernUMI) == 0))
sigMorten_zeroHalpern <- intersect(sig.g.morten, zeroGenes)
length(sigMorten_zeroHalpern)
detect.gene.morten <- apply(droplet.norm[useDroplet.genes,], 1, function(x) sum(x!=0))
sigMorten_zeroHalpern <- intersect(sigMorten_zeroHalpern, names(which(detect.gene.morten >= 10)))
length(sigMorten_zeroHalpern)

