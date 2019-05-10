setwd("LiverSpatialCompare")

# Genes that were concatonated in UMI dataset:

load("RDATA/dataReady_bothData_genesMapped.RData")

load("RDATA/analysis_WaveCrest_Morten.RData")
data.norm.order <- data.norm[,wc.order]
layerSplit <- split(1:ncol(data.norm.order), cut(seq_along(1:ncol(data.norm.order)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))

layerMeans.morten <- t(apply(data.norm.order, 1, function(x) {
  return(tapply(x, layerSplit, mean))


# First, which genes in UMI are represented by separate genes in full-length data?:
catSet <- unique(geneSet[which(duplicated(geneSet$HalpernGene)),2])
length(catSet)

head(catSet)
# How many unique genes is this?


totalgenes <- subset(geneSet, HalpernGene %in% catSet)
dim(totalgenes)


# Plot Mup and Ugt specifically (they are the biggest gene sets):

pdf("PLOTS/concatonatedGenesHalpern_MupGenes.pdf", height=6, width=15)

par(mfrow=c(3,4),  mar=c(5,5,3,1))  

genes <- subset(geneSet, HalpernGene == catSet[grep("Mup",catSet)])$MortenGene

for(j in 1:length(genes)) {

gene.h <- catSet[grep("Mup",catSet)]
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
plot(0,0, col="white", pch=20, cex.axis=1.8, cex.lab=1.8,bty = 'n', xlim=c(1,9),
     cex.main=2, main=geneX, ylim=c(0,1), ylab="Scaled Expressed", xlab="Zonation Group", xaxt='n', yaxt='n')
axis(2, at=seq(0,1, by=.2), label=seq(0,1, by=.2), lwd=2,cex.axis=1.8)
axis(1, at=1:9, label=1:9, cex.axis=1.8, lwd=2)


FIT = smooth.spline(1:9, rescaleY.morten.means,
                    control.spar=list(low=.2, high=.5))
lines(FIT$x, FIT$y, lwd=3, col="#fc8d59")

FIT = smooth.spline(1:9, rescaleY.halpern.means, 
                    control.spar=list(low=.2, high=.5))
lines(FIT$x, FIT$y, lwd=3, col="#91bfdb")

points(1:9, rescaleY.morten.means, col="#fc8d59", pch=95, cex=3)
points(1:9, rescaleY.halpern.means, col="#91bfdb", pch=95, cex=3)


}
dev.off()









pdf("PLOTS/concatonatedGenesHalpern_UgtGenes.pdf", height=4, width=15)

par(mfrow=c(2,4),  mar=c(5,5,3,1))  

genes <- subset(geneSet, HalpernGene == catSet[grep("Ugt",catSet)])$MortenGene

for(j in 1:length(genes)) {

gene.h <- catSet[grep("Mup",catSet)]
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
plot(0,0, col="white", pch=20, cex.axis=1.8, cex.lab=1.8,bty = 'n', xlim=c(1,9),
     cex.main=2, main=geneX, ylim=c(0,1), ylab="Scaled Expressed", xlab="Zonation Group", xaxt='n', yaxt='n')
axis(2, at=seq(0,1, by=.2), label=seq(0,1, by=.2), lwd=2,cex.axis=1.8)
axis(1, at=1:9, label=1:9, cex.axis=1.8, lwd=2)


FIT = smooth.spline(1:9, rescaleY.morten.means,
                    control.spar=list(low=.2, high=.5))
lines(FIT$x, FIT$y, lwd=3, col="#fc8d59")

FIT = smooth.spline(1:9, rescaleY.halpern.means, 
                    control.spar=list(low=.2, high=.5))
lines(FIT$x, FIT$y, lwd=3, col="#91bfdb")

points(1:9, rescaleY.morten.means, col="#fc8d59", pch=95, cex=3)
points(1:9, rescaleY.halpern.means, col="#91bfdb", pch=95, cex=3)

}
dev.off()


pdf("PLOTS/concatonatedGenesHalpern_UgtGenes_legendOnly.pdf", height=6, width=15)


plot(0,0, col="white", pch=20, cex.axis=1.8, cex.lab=1.8,bty = 'n', xlim=c(1,9),
     cex.main=2, main=geneX, ylim=c(0,1), ylab="", xlab="", xaxt='n', yaxt='n')

legend('topright', c("Full-length", "UMI"), col=c("#fc8d59", "#91bfdb"), lty=1, lwd=3, bty='n')


dev.off()


