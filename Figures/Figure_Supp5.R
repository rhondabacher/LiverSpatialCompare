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

load("RDATA/analysis_Monocle_Droplet.RData")
# Split 10X into 9 layers just like Halpern data:
layerSplit <- split(1:length(pt_data.droplet), cut(seq_along(1:length(pt_data.droplet)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))
# Use mean because so many zeros.
data.norm.order.droplet <- droplet.norm[,order(pt_data.droplet)]
layerMeans.droplet <- t(apply(data.norm.order.droplet, 1, function(x) {
  return(tapply(x, layerSplit, mean))
}))



# First, which genes in UMI are represented by separate genes in full-length data?:
catSet <- unique(geneSet[which(duplicated(geneSet$HalpernGene)),2])
length(catSet)

head(catSet)
# How many unique genes is this?


totalgenes <- subset(geneSet, HalpernGene %in% catSet)
dim(totalgenes)


# Plot Ugt --one of the biggest gene sets:



pdf("PLOTS/concatonatedGenesHalpern_Ugt1Genes.pdf", height=10, width=8)

par(mfrow=c(4,2),  mar=c(5,5,3,1))  

genes <- subset(geneSet, HalpernGene == catSet[grep("Ugt1",catSet)])$Gene

for(j in 1:length(genes)) {

gene.h <- catSet[grep("Ugt1",catSet)]
gene.m <- genes[j]

# Scale data
getMeans.h <- layerMeans[gene.h, ]
rescaleY.halpern.means <- ((1 - 0)/(max(getMeans.h) - min(getMeans.h)))*(getMeans.h - min(getMeans.h)) + 0

getMeans.m <- layerMeans.morten[gene.m,]
rescaleY.morten.means <- ((1 - 0)/(max(getMeans.m) - min(getMeans.m)))*(getMeans.m - min(getMeans.m)) + 0

getMeans.d <- layerMeans.droplet[gene.m, ]
rescaleY.droplet.means <- ((1 - 0)/(max(getMeans.d) - min(getMeans.d)))*(getMeans.d - min(getMeans.d)) + 0


geneX <- paste0(gene.m)
plot(0,0, col="white", pch=20, cex.axis=1.5, cex.lab=1.5,bty = 'n', xlim=c(1,9),
     cex=1, main=geneX, ylim=c(0,1), ylab="Scaled Expressed", xlab="Zonation Group", xaxt='n', yaxt='n')
axis(2, at=seq(0,1, by=.5), label=seq(0,1, by=.5), lwd=2,cex.axis=1.5)
axis(1, at=1:9, label=1:9, cex.axis=1.5, lwd=2)

if(sum(getMeans.m) != 0 ){
  FIT = smooth.spline(1:9, rescaleY.morten.means, df=4)
  lines(FIT$x, FIT$y, lwd=3, col="#fc8d59")
  points(1:9, rescaleY.morten.means, col="#fc8d59", pch=95, cex=3)
}
if(sum(getMeans.h) != 0 ){
  FIT = smooth.spline(1:9, rescaleY.halpern.means, df=4)
  lines(FIT$x, FIT$y, lwd=3, col="#91bfdb")
  points(1:9, rescaleY.halpern.means, col="#91bfdb", pch=95, cex=3)
}
if(sum(getMeans.d) != 0 ){
  FIT = smooth.spline(1:9, rescaleY.droplet.means, df=4)
  lines(FIT$x, FIT$y, lwd=3, col="chartreuse3")
  points(1:9, rescaleY.droplet.means, col="chartreuse3", pch=95, cex=3)
}

}
dev.off()



