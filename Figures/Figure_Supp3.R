setwd("LiverSpatialCompare")

# Rereate plots from Halperns figures of non-monotonic genes:

load("RDATA/dataReady_bothData_genesMapped.RData")
load("RDATA/analysis_WaveCrest_Morten.RData")

data.norm.order <- data.norm[,wc.order]


## Plot them as an overlay with Halpern but scale both between 0 and 1 since the data are on different scales.

# Split Morten into 9 layers just like Halpern data:
layerSplit <- split(1:ncol(data.norm.order), cut(seq_along(1:ncol(data.norm.order)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))

layerMeans.morten <- t(apply(data.norm.order, 1, function(x) {
  return(tapply(x, layerSplit, median))
}))

allGenes <- c("Hamp", "Igfbp2", "Cyp8b1", "Mup3")


pdf(paste0("PLOTS/morten_geneExp_Ordered_scatterFit_overlayHalpernfromFig3_SuppFig2.pdf"), height=6, width=10, useDingbats=F)
par(mfrow=c(2,2),  mar=c(5,5,2,1))
for(geneX in allGenes) {

  # Scale data
  getMeans.h <- layerMeans[geneX, ]
  rescaleY.halpern.means <- ((1 - 0)/(max(getMeans.h) - min(getMeans.h)))*(getMeans.h - min(getMeans.h)) + 0
  
  getMeans.m <- layerMeans.morten[geneX,]
  rescaleY.morten.means <- ((1 - 0)/(max(getMeans.m) - min(getMeans.m)))*(getMeans.m - min(getMeans.m)) + 0
  
  plot(0,0, col="white", pch=20, cex.axis=1.5, cex.lab=1.5,bty = 'n', xlim=c(1,9),
       cex=1, main=geneX, ylim=c(0,1), ylab="Scaled Expressed", xlab="Zonation Group", xaxt='n', yaxt='n')
  axis(2, at=seq(0,1, by=.5), label=seq(0,1, by=.5), lwd=2,cex.axis=1.5)
  axis(1, at=1:9, label=1:9, cex.axis=1.5, lwd=2)
  
  FIT = smooth.spline(1:9, rescaleY.morten.means, df=4)
  lines(FIT$x, FIT$y, lwd=3, col="#fc8d59")
  
  FIT = smooth.spline(1:9, rescaleY.halpern.means, df=4)
  lines(FIT$x, FIT$y, lwd=3, col="#91bfdb")
  
  points(1:9, rescaleY.morten.means, col="#fc8d59", pch=95, cex=3)
  points(1:9, rescaleY.halpern.means, col="#91bfdb", pch=95, cex=3)
  
  # Decide where to put legenes
  if (rescaleY.morten.means[9] < .5){
    legend('topright', c("Smart-seq", "MARS-seq"), col=c("#fc8d59", "#91bfdb"), lty=1, lwd=3)
  } else {
    legend('bottomright', c("Smart-seq", "MARS-seq"), col=c("#fc8d59", "#91bfdb"), lty=1, lwd=3)
  }

}
dev.off()

 

allGenes <- c("Cyp7a1","Hsd3b7","Cyp8b1", "Cyp27a1", "Baat")



pdf(paste0("PLOTS/morten_geneExp_Ordered_scatterFit_overlayHalpernfromFig4_SuppFig3.pdf"), height=3, width=12, useDingbats=F)
par(mfrow=c(1,5),  mar=c(5,5,2,1))
for(geneX in allGenes) {

  gene.m <- geneX
  gene.h <- geneSet[which(geneSet$MortenGene == geneX), 2]
  
  # Scale data
  getMeans.h <- layerMeans[gene.h, ]
  rescaleY.halpern.means <- ((1 - 0)/(max(getMeans.h) - min(getMeans.h)))*(getMeans.h - min(getMeans.h)) + 0
  
  getMeans.m <- layerMeans.morten[gene.m,]
  rescaleY.morten.means <- ((1 - 0)/(max(getMeans.m) - min(getMeans.m)))*(getMeans.m - min(getMeans.m)) + 0
  
  plot(0,0, col="white", pch=20, cex.axis=1.5, cex.lab=1.5,bty = 'n', xlim=c(1,9),
       cex=1, main=geneX, ylim=c(0,1), ylab="Scaled Expressed", xlab="Zonation Group", xaxt='n', yaxt='n')
  axis(2, at=seq(0,1, by=.5), label=seq(0,1, by=.5), lwd=2,cex.axis=1.5)
  axis(1, at=1:9, label=1:9, cex.axis=1.5, lwd=2)
  
  FIT = smooth.spline(1:9, rescaleY.morten.means, df=4)
  lines(FIT$x, FIT$y, lwd=3, col="#fc8d59")
  
  FIT = smooth.spline(1:9, rescaleY.halpern.means, df=4)
  lines(FIT$x, FIT$y, lwd=3, col="#91bfdb")
  
  points(1:9, rescaleY.morten.means, col="#fc8d59", pch=95, cex=3)
  points(1:9, rescaleY.halpern.means, col="#91bfdb", pch=95, cex=3)
  
  # Decide where to put legenes
  if (rescaleY.morten.means[9] < .5){
    legend('topright', c("Smart-seq", "MARS-seq"), col=c("#fc8d59", "#91bfdb"), lty=1, lwd=3)
  } else {
    legend('bottomright', c("Smart-seq", "MARS-seq"), col=c("#fc8d59", "#91bfdb"), lty=1, lwd=3)
  }

}
dev.off()