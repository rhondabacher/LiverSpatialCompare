setwd("LiverSpatialCompare")

# Compare isoforms in full-length data:

load("RDATA/dataReady_MortenIsoforms.RData")

load("RDATA/analysis_WaveCrest_Morten.RData")


# Apply the same ordering to the isoform processed data.

data.iso.norm.order <- data.iso.norm[,wc.order]

# Calculate means for Morten's data for plotting:
layerSplit <- split(1:ncol(data.iso.norm.order), cut(seq_along(1:ncol(data.iso.norm.order)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))

layerMeans.morten <- t(apply(data.iso.norm.order, 1, function(x) {
  return(tapply(x, layerSplit, mean))
}))


txs <- data.iso.all[,1]

genes <- unique(subset(data.iso.all[,1:2], data.iso.all[,1] %in% txs)[,2])

txsList <-list()
for(i in 1:length(genes)) {
  txsList[[i]] <- intersect(subset(data.iso.all[,1:2], data.iso.all[,2] == genes[i])[,1], txs)
}
names(txsList) <- genes

# Plot Romo1 isoforms in full-length data:
library(ggplot2)
txX <- txsList[["Romo1"]]
txX1 <- txX[1]
txX2 <- txX[2]

getMeans.h <- log2(data.iso.norm.order[txX1,]+1)
getMeans.m <- log2(data.iso.norm.order[txX2,]+1)

MAX <- max(getMeans.h, getMeans.m)+1

pdf("PLOTS/scatterIsoforms_forGene_Romo1.pdf", height=5, width=6)
par(mfrow=c(1,1), mar=c(5,5,2,1))
plot(1:length(wc.order), getMeans.h, col="white", pch=20, cex.axis=1.8, 
     cex.lab=2,bty = 'n', xlim=c(0,length(wc.order)),
     cex=1, main="Romo1", cex.main=2,
     ylim=c(0,MAX), xlab="Cells", 
     ylab="Log2 Expression", xaxt='n', yaxt='n')
axis(2, at=round(seq(0,MAX, length.out = 10)),lwd=2,cex.axis=1.8)
axis(1, at=round(seq(1,length(wc.order), length.out = 10)), cex.axis=1.8, lwd=2)
points(1:length(wc.order), getMeans.h, pch=16, col=alpha("lightgoldenrod3",.6))
points(1:length(wc.order), getMeans.m, pch=16, col=alpha("indianred1",.6))

FIT = smooth.spline(1:length(wc.order), getMeans.m,
                    control.spar=list(low=.5, high=1))
lines(FIT$x, FIT$y, lwd=3, col="indianred1")

FIT = smooth.spline(1:length(wc.order), getMeans.h, 
                    control.spar=list(low=.5, high=1))
lines(FIT$x, FIT$y, lwd=3, col="lightgoldenrod3")
legend('bottomleft', c(txX2, txX1), lwd=2, col=c("indianred1", "lightgoldenrod3"), bty='n')
dev.off()




########################################


# Plot Mup1, Nadk2:

library(ggplot2)
txX <- txsList[["Mup1"]]
txX1 <- txX[1]
txX2 <- txX[2]

getMeans.h <- log2(data.iso.norm.order[txX1,]+1)
getMeans.m <- log2(data.iso.norm.order[txX2,]+1)

MAX <- max(getMeans.h, getMeans.m)+1

pdf("PLOTS/scatterIsoforms_forGene_Mup1.pdf", height=5, width=6)
par(mfrow=c(1,1), mar=c(5,5,2,1))
plot(1:length(wc.order), getMeans.h, col="white", pch=20, cex.axis=1.8, 
     cex.lab=2,bty = 'n', xlim=c(0,length(wc.order)),
     cex=1, main="Mup1", cex.main=2,
     ylim=c(0,MAX), xlab="Cells", 
     ylab="Log2 Expression", xaxt='n', yaxt='n')
axis(2, at=round(seq(0,MAX, length.out = 10)),lwd=2,cex.axis=1.8)
axis(1, at=round(seq(1,length(wc.order), length.out = 10)), cex.axis=1.8, lwd=2)
points(1:length(wc.order), getMeans.h, pch=16, col=alpha("lightgoldenrod3",.6))
points(1:length(wc.order), getMeans.m, pch=16, col=alpha("indianred1",.6))

FIT = smooth.spline(1:length(wc.order), getMeans.m,
                    control.spar=list(low=.5, high=1))
lines(FIT$x, FIT$y, lwd=3, col="indianred1")

FIT = smooth.spline(1:length(wc.order), getMeans.h, 
                    control.spar=list(low=.5, high=1))
lines(FIT$x, FIT$y, lwd=3, col="lightgoldenrod3")
legend('bottomright', c(txX2, txX1), lwd=2, col=c("indianred1", "lightgoldenrod3"), bty='n')
dev.off()





library(ggplot2)
txX <- txsList[["Nadk2"]]
txX1 <- txX[1]
txX2 <- txX[2]

getMeans.h <- log2(data.iso.norm.order[txX1,]+1)
getMeans.m <- log2(data.iso.norm.order[txX2,]+1)

MAX <- max(getMeans.h, getMeans.m)+2

pdf("PLOTS/scatterIsoforms_forGene_Nadk2.pdf", height=5, width=6)
par(mfrow=c(1,1), mar=c(5,5,2,1))
plot(1:length(wc.order), getMeans.h, col="white", pch=20, cex.axis=1.8, 
     cex.lab=2,bty = 'n', xlim=c(0,length(wc.order)),
     cex=1, main="Nadk2", cex.main=2,
     ylim=c(0,MAX), xlab="Cells", 
     ylab="Log2 Expression", xaxt='n', yaxt='n')
axis(2, at=round(seq(0,MAX, length.out = 10)),lwd=2,cex.axis=1.8)
axis(1, at=round(seq(1,length(wc.order), length.out = 10)), cex.axis=1.8, lwd=2)
points(1:length(wc.order), getMeans.h, pch=16, col=alpha("lightgoldenrod3",.6))
points(1:length(wc.order), getMeans.m, pch=16, col=alpha("indianred1",.6))

FIT = smooth.spline(1:length(wc.order), getMeans.m,
                    control.spar=list(low=1, high=1))
lines(FIT$x, FIT$y, lwd=3, col="indianred1")

FIT = smooth.spline(1:length(wc.order), getMeans.h, 
                    control.spar=list(low=1, high=1))
lines(FIT$x, FIT$y, lwd=3, col="lightgoldenrod3")
legend('topleft', c(txX2, txX1), lwd=2, col=c("indianred1", "lightgoldenrod3"), bty='n')
dev.off()

