setwd("LiverSpatialCompare")

# Individual genes plots that were stained by Morten

load("RDATA/dataReady_bothData_genesMapped.RData")
load("RDATA/analysis_WaveCrest_Morten.RData")

data.norm.order <- data.norm[,wc.order]

allGenes <- c("Cyp1a2", "Cyp2e1", "Oat", "Rgn", "Aldh3a2", "Tbx3", "Cyp2f2", "Hal")


## Plot them as an overlay with Halpern but scale both between 0 and 1 since the data are on different scales.

# Split Morten into 9 layers just like Halpern data:
layerSplit <- split(1:ncol(data.norm.order), cut(seq_along(1:ncol(data.norm.order)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))

layerMeans.morten <- t(apply(data.norm.order, 1, function(x) {
  return(tapply(x, layerSplit, mean))
}))


pdf(paste0("PLOTS/morten_geneExp_Ordered_scatterFit_",i,"_overlayHalpern_ALL_Fig2.pdf"), height=4, width=10)  
par(mfrow=c(2,4),  mar=c(5,5,2,1))
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
  
  FIT = smooth.spline(1:9, rescaleY.morten.means,
                      control.spar=list(low=.2, high=.5))
  lines(FIT$x, FIT$y, lwd=3, col="#fc8d59")
  
  FIT = smooth.spline(1:9, rescaleY.halpern.means, 
                      control.spar=list(low=.2, high=.5))
  lines(FIT$x, FIT$y, lwd=3, col="#91bfdb")
  
  points(1:9, rescaleY.morten.means, col="#fc8d59", pch=95, cex=3)
  points(1:9, rescaleY.halpern.means, col="#91bfdb", pch=95, cex=3)
  
  # Decide where to put legenes
  if (rescaleY.morten.means[1] > rescaleY.morten.means[9]){
    legend('topright', c("Full-length", "UMI"), col=c("#fc8d59", "#91bfdb"), lty=1, lwd=3)
  } else {
    legend('bottomright', c("Full-length", "UMI"), col=c("#fc8d59", "#91bfdb"), lty=1, lwd=3)
  }

}
dev.off()


## Calculate the correlation of gene expression across the 9 layers:
getCorr <- c()
for(i in 1:length(allGenes)) {
  
  gene.m <- allGenes[i]
  gene.h <- allGenes[i]
  
  getMeans.m <- layerMeans.morten[gene.m,]
  getMeans.h <- layerMeans[gene.h,]
  
  rescale.halpern.means <- ((1 - 0)/(max(getMeans.h) - min(getMeans.h)))*(getMeans.h - min(getMeans.h)) + 0
  rescale.morten.means <- ((1 - 0)/(max(getMeans.m) - min(getMeans.m)))*(getMeans.m - min(getMeans.m)) + 0
  
  if(any(is.na(rescale.morten.means)) | any(is.na(rescale.halpern.means))){
    corr.h.m <- NA
  } else{
    corr.h.m <- cor(as.vector(rescale.morten.means), as.vector(rescale.halpern.means), method='spearman')
  }
  getCorr <- c(getCorr, corr.h.m)
}
names(getCorr) <- allGenes


median(getCorr)
