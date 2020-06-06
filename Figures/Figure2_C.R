# Individual genes plots that were stained by Morten

setwd("scSpatialReconstructCompare-Paper")

load("RDATA/dataReady_bothData_genesMapped.RData")
load("RDATA/analysis_WaveCrest_Morten.RData")
data.norm.order <- data.norm[,wc.order]

allGenes <- c("Cyp1a2", "Cyp2e1", "Oat", "Rgn", "Aldh3a2", "Tbx3", "Cyp2f2", "Hal")


for(i in allGenes) {
pdf(paste0("PLOTS/FORFIGS/morten_geneExp_Ordered_scatterFit_",i,"_Fig2.pdf"), height=2, width=5, useDingbats=F)  
par(mfrow=c(1,1),  mar=c(2,5,2,1))
  geneX <- i
  
  currentY.m <- log2(data.norm[geneX,wc.order]+1)
  
  plot(0,0, col="white", pch=20, cex.axis=1.3, cex.lab=1.3,bty = 'n', 
            xlim=c(1,length(wc.order)),
       cex=1, cex.main=2, main=geneX, ylim=c(0,16), 
       ylab="log2 expression", xlab="", xaxt='n', yaxt='n')
  axis(2, at=seq(0,16, by=4),  lwd=2,cex.axis=1.5, cex.lab=2)
  axis(1, at=seq(1,length(wc.order), by=13), cex.axis=1.5, lwd=2)
  
  points(1:length(wc.order), currentY.m, col="#fc8d59", pch=16, cex=1)
  
  FIT = smooth.spline(1:length(wc.order), currentY.m, df=4)
  lines(FIT$x, FIT$y, lwd=3, col="#fc8d59")
  
dev.off()
}
 