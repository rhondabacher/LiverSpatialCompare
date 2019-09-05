setwd("LiverSpatialCompare")

load("RDATA/dataReady_Morten.RData")


###  Now do reordering:
library(WaveCrest)

marker.in <- get(load("RDATA/markerGenes.RData"))
markers <- marker.in[[1]]

mk.groups <- factor(marker.in[[2]])
names(mk.groups) <- markers

# Scale data as input to Wave-Crest:
data.norm.scale <- data.norm
data.norm.scale[which(data.norm.scale < 1)] <- 1
data.norm.scale <- log2(data.norm.scale)
data.norm.scale <- Rescale(data.norm.scale)

# Put back in genes that have all zeros which were removed by the Rescale function:
data.zeroGenes <- data.norm[setdiff(rownames(data.norm),rownames(data.norm.scale)),]


# Wave-Crest will reorder based on the marker genes only:
data.markers <- data.norm.scale[markers,]

set.seed(100)
rand.order <- sample(1:ncol(data.norm.scale),ncol(data.norm.scale))

library(gplots)
library(RColorBrewer)
heatcols <- rev(colorRampPalette(brewer.pal(11,"RdBu"))(100))
heatcol2 <- heatcols[c(1:33,seq(34,66,4), 67:100)]
groupColors <- c("red","blue")
cond.num <- rep(1,ncol(data.norm.scale))

pdf("PLOTS/heatmap_markers_Halpern.pdf")
heatmap.2(data.markers, trace="none",col=heatcol2)
dev.off()

##################################################################
##################################################################

# Use degree 1 for linear

initialOrder <- ImpTC(data.markers, rand.order, cond.num, Ndg = 1)
optOrder <- Opt2TC(data.markers, 100000, initialOrder, cond.num, Ndg = 1)

wc.order <- optOrder[[1]]

pdf("PLOTS/heatmap_markers_Halpern_waveCrest-Order.pdf")
heatmap.2(data.markers[,wc.order], Colv=F, dendrogram='row', 
          trace="none",RowSideColors=groupColors[mk.groups],
          col=heatcol2)
dev.off()



# Remove ERCC genes
data.norm.scale <- data.norm.scale[-grep("ERCC-", rownames(data.norm.scale)),]

fitall.log.mse <- sapply(1:nrow(data.norm.scale), function(j){
  tt = lm(data.norm.scale[j,wc.order] ~ poly(1:ncol(data.norm.scale),1))
  t2 = mean(tt$residuals^2)
  return(t2)  
})
names(fitall.log.mse) <- rownames(data.norm.scale)

fitall.log.mse <- sort(fitall.log.mse, decreasing=F)


# Let's also get the p-value and a slope:

fitall.pvals <- sapply(1:nrow(data.norm.scale),function(j){
  tt = lm(data.norm.scale[j,wc.order]~poly(1:ncol(data.norm.scale),1))
  t2 = summary(tt)$coefficients[2,4]  
  return(t2 )
})
names(fitall.pvals) <- rownames(data.norm.scale)

adjusted.pvals <- p.adjust(fitall.pvals, method='fdr')


# Here use a strict cutoff of FDR 1%
SIG <- names(which(adjusted.pvals < .01))
SIG <- names(sort(fitall.log.mse[SIG]))

pdf("PLOTS/SupplementaryFigure_allSIG_WaveCrest_geneData.pdf",height=12, width=10)
par(mfrow=c(5,2),mar=c(5,5,5,5))
for(i in 1:length(SIG)) {
  plot(data.norm.scale[SIG[i], wc.order],main=SIG[i],
       xlab="Recovered order",ylab="log2 expression", pch=20)
  FIT = smooth.spline(1:ncol(data.norm.scale), data.norm.scale[SIG[i],wc.order], 
                      control.spar=list(low=.8))
  lines(FIT$x, FIT$y, lwd=2)
}
dev.off()



# Since poly = 1, we are just fitting straight lines.
fitall.slopes <- sapply(1:nrow(data.norm.scale),function(j){
  tt=lm(data.norm.scale[j,wc.order]~poly(1:ncol(data.norm.scale),1))
  t2=summary(tt)$coefficients[2,1]  
  return(t2)
})
names(fitall.slopes) <- rownames(data.norm.scale)

fitall.std.error <- sapply(1:nrow(data.norm.scale),function(j){
  tt=lm(data.norm.scale[j,wc.order]~poly(1:ncol(data.norm.scale),1))
  t2=summary(tt)$coefficients[2,2]  
  return(t2)
})
names(fitall.std.error) <- rownames(data.norm.scale)


sum(fitall.slopes > 0)

Genes <- rownames(data.norm.scale)
done.data = data.frame(
                       Gene = Genes,
                       MSE = round(fitall.log.mse[Genes],2),
                       Slope = fitall.slopes[Genes],
                       'Std Error' = fitall.std.error[Genes],
                       pvalue = fitall.pvals[Genes],
                       'FDR adjusted pvalue' = adjusted.pvals[Genes])
head(done.data)
done.data$ZonationAxis <- "NotSignificant"

done.data$ZonationAxis[which(done.data$FDR.adjusted.pvalue < .05 & done.data$Slope < 0)] <- "PV"
done.data$ZonationAxis[which(done.data$FDR.adjusted.pvalue < .05 & done.data$Slope >= 0)] <- "PP"


write.csv(done.data, 
          file="OUT/SupplementaryFile1.csv", row.names=F)

save.image("RDATA/analysis_WaveCrest_Morten.RData")


write.csv(data.norm[,wc.order], "OUT/normalized_WaveCrestOrder.csv")


