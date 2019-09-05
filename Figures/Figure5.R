setwd("LiverSpatialCompare")

# Genes that were concatonated in UMI dataset:

load("RDATA/dataReady_bothData_genesMapped.RData")

load("RDATA/analysis_WaveCrest_Morten.RData")
data.norm.order <- data.norm[,wc.order]
layerSplit <- split(1:ncol(data.norm.order), cut(seq_along(1:ncol(data.norm.order)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))

layerMeans.morten <- t(apply(data.norm.order, 1, function(x) {
  return(tapply(x, layerSplit, median))
}))


# First, which genes in UMI are represented by separate genes in full-length data?:
catSet <- unique(geneSet[which(duplicated(geneSet$HalpernGene)),2])
length(catSet)

head(catSet)
# How many unique genes is this?


totalgenes <- subset(geneSet, HalpernGene %in% catSet)
dim(totalgenes)


# Plot Ugt --one of the biggest gene sets:



pdf("PLOTS/concatonatedGenesHalpern_UgtGenes.pdf", height=10, width=8)

par(mfrow=c(4,2),  mar=c(5,5,3,1))  

genes <- subset(geneSet, HalpernGene == catSet[grep("Ugt",catSet)])$MortenGene

for(j in 1:length(genes)) {

gene.h <- catSet[grep("Ugt",catSet)]
gene.m <- genes[j]

getMeans.h <- layerMeans[gene.h, ]
if (sum(getMeans.h > 0)) {
  rescaleY.halpern.means <- ((1 - 0)/(max(getMeans.h) - min(getMeans.h)))*(getMeans.h - min(getMeans.h)) + 0
} else {
  rescaleY.halpern.means <- rep(0, 9)
}
getMeans.m <- layerMeans.morten[gene.m,]
if (sum(getMeans.m > 0)) {
  rescaleY.morten.means <- ((1 - 0)/(max(getMeans.m) - min(getMeans.m)))*(getMeans.m - min(getMeans.m)) + 0
} else {
  rescaleY.morten.means <- rep(0, 9)
}

geneX <- paste0(gene.m)
plot(0,0, col="white", pch=20, cex.axis=2, cex.lab=2.3,bty = 'n', xlim=c(1,9),
     cex.main=2.4, main=geneX, ylim=c(0,1), ylab="Scaled Expressed", xlab="Zonation Group", xaxt='n', yaxt='n')
axis(2, at=seq(0,1, by=.2), label=seq(0,1, by=.2), lwd=2,cex.axis=2)
axis(1, at=1:9, label=1:9, cex.axis=2, cex.lab=2.3, lwd=2)


FIT = smooth.spline(1:9, rescaleY.morten.means, df=4)
lines(FIT$x, FIT$y, lwd=3, col="#fc8d59")

FIT = smooth.spline(1:9, rescaleY.halpern.means, df=4)
lines(FIT$x, FIT$y, lwd=3, col="#91bfdb")

points(1:9, rescaleY.morten.means, col="#fc8d59", pch=95, cex=3)
points(1:9, rescaleY.halpern.means, col="#91bfdb", pch=95, cex=3)

}
dev.off()



