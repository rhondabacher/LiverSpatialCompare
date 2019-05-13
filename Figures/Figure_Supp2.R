# For any Kegg category, check correlation and heatmap:

setwd("LiverSpatialCompare")

## Load in datasets:
load("RDATA/dataReady_bothData_genesMapped.RData")
rownames(layerMeans) <- gsub(" ", "", rownames(layerMeans), fixed=TRUE)

load("RDATA/analysis_WaveCrest_Morten.RData")
data.norm.order <- data.norm[,wc.order]

# Calculate means for Morten's data for plotting:
layerSplit <- split(1:ncol(data.norm.order), cut(seq_along(1:ncol(data.norm.order)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))

layerMeans.morten <- t(apply(data.norm.order, 1, function(x) {
  return(tapply(x, layerSplit, mean))
}))

scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
  }
## Calculate the correlation of gene expression across the 9 layers:
getCorr <- c()
for(i in 1:nrow(geneSet)) {
  
  gene.m <- geneSet$MortenGene[i]
  gene.h <- geneSet$HalpernGene[i]
  
  getMeans.m <- layerMeans.morten[gene.m,]
  getMeans.h <- layerMeans[gene.h,]
  
  rescale.halpern.means <- scale01(getMeans.h)
  rescale.morten.means <- scale01(getMeans.m)
  
  
  if(any(is.na(rescale.morten.means)) | any(is.na(rescale.halpern.means))){
    corr.h.m <- NA
  } else{
    corr.h.m <- cor(as.vector(rescale.morten.means), as.vector(rescale.halpern.means), method='spearman')
  }
  getCorr <- c(getCorr, corr.h.m)
}
names(getCorr) <- geneSet[,2]

# Functions to draw histograms:
polyCurve <- function(x, y, from, to, n = 100, miny,
                      col = "red", border = col) {
  drawPoly <- function(fun, from, to, n = 100, miny, col, border) {
    Sq <- seq(from = from, to = to, length = n)
    polygon(x = c(Sq[1], Sq, Sq[n]),
            y = c(miny, fun(Sq), miny),
            col = col, border = border)
  }
  lf <- length(from)
  stopifnot(identical(lf, length(to)))
  if(length(col) != lf)
    col <- rep(col, length.out = lf)
  if(length(border) != lf)
    border <- rep(border, length.out = lf)
  if(missing(miny))
    miny <- min(y)
  interp <- approxfun(x = x, y = y)
  mapply(drawPoly, from = from, to = to, col = col, border = border,
         MoreArgs = list(fun = interp, n = n, miny = miny))
  invisible()
}

PushOL<- function (Data, qt1 = 0.05, qt2 = 0.95, PushHigh = T, PushLow = T) 
{
  Q5 = apply(Data, 1, function(i) quantile(i, qt1))
  Q95 = apply(Data, 1, function(i) quantile(i, qt2))
  DataSc2 = Data
  for (i in 1:nrow(Data)) {
    if (PushLow) 
      DataSc2[i, which(DataSc2[i, ] < Q5[i])] = Q5[i]
    if (PushHigh) 
      DataSc2[i, which(DataSc2[i, ] > Q95[i])] = Q95[i]
  }
  DataSc2
}


firstPlot <- function(MAIN, TOP) {
  dens1 <- density(na.omit(getCorr), from=-1, to=1, bw=.3)
  par(mar=c(5,5,2,1))
  plot(dens1$x, dens1$y, cex.axis=1.5, 
       cex.lab=1.5, bty = 'n', xlim=c(-1.1,1.1), lwd=2,
       cex=1, main=MAIN, ylim=c(0,(TOP*1.1)), col="black", type="l", cex.main=1.5,
       xlab="Correlation", ylab="Density", xaxt='n', yaxt='n')
  axis(2, at=round(seq(0,TOP, length.out=4),1), lwd=2,cex.axis=2)
  axis(1, at=seq(-1,1, by=.5), cex.axis=1.5, lwd=2,cex.axis=2)
  polyCurve(dens1$x, dens1$y, from = -1, to = 1,miny=0,
            col = alpha("gray", .2), border = alpha("gray", .2))
  
}
quickScatter <- function(allGenes) {
  for(i in allGenes) {

    geneX <- i
  
    # Scale data
    getMeans.h <- layerMeans[geneX, ]
    rescaleY.halpern.means <- scale01(getMeans.h)
  
    getMeans.m <- layerMeans.morten[geneX,]
    rescaleY.morten.means <- scale01(getMeans.m)
  
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
}

library(scales)
library(viridis)
makeCorPlots <- function(KEGG, NAME, catGenes, pushLow=0, pushHigh=1){

  keggCat1 <- bitr_kegg(KEGG, "Path", "ncbi-geneid", "mmu")[,2]
  
  sub2 <- bitr(keggCat1, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")[,2]
  sub2 <- intersect(geneSet[,1], sub2)
  
  sub3 <- bitr(catGenes, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")[,2]
  
  dens2 <- density(na.omit(getCorr[sub2]), from=-1, to=1, bw=.3)
  dens3 <- density(na.omit(getCorr[sub3]), from=-1, to=1, bw=.3)
  
  TOP <- round(max(dens2$y, dens3$y),1)

  data.norm.order.rmout <- log(PushOL(data.norm.order, qt1 = .02, qt2 = 0.98)+1)

  top.res.fitted <- t(apply(log(data.norm.order.rmout[sub3,]+1), 1, function(x) {
    FIT = smooth.spline(1:ncol(data.norm.order.rmout),x, 
                        control.spar=list(low=.6, high=.7))
    return(FIT$y)
  
  }))
  
  useY <- t(apply(top.res.fitted, 1, function(x) tapply(x, layerSplit, mean)))
  useOrder <- names(sort(apply(useY, 1, which.max)))
  
  dataForPlot <- t(base::scale(t((top.res.fitted))))
  dataForPlot <- PushOL(dataForPlot, qt1 = pushLow, qt2 = pushHigh)
  
  dataForPlot <- dataForPlot[useOrder,]
  MIN = min(dataForPlot)
  MAX = max(dataForPlot)
  par(mar=c(5,2,2,5))
 
  heatcols2 <- cividis(100)
  # heatcols2 <- heatcols2[c(1:35, seq(40,80, by=5), 95:100)]
  breaks <- seq(MIN, MAX, length.out=length(heatcols2)+1)

  # par(mar=c(4,2,1,5))
  #   image(t(dataForPlot[useOrder,]), col = heatcols2,
  #           breaks = breaks,
  #           xaxt='n', yaxt='n')
            
  MAT = matrix(c(1,2,3, 10, 1,2,3, 11, 4,5,6, 11, 7:9, 11), nrow=4, ncol=4, byrow=T)
  layout(MAT, width=c(.5,.5,.5,1), height=c(.5,.5,1,1))
  # layout.show(11)

  
  par(mar=c(5,5,2,7))
  firstPlot(NAME, TOP = TOP)
  segments(0,0, 0,TOP, lty=2, lwd=1.5)
  lines(dens2$x, dens2$y,lwd=2)
  polyCurve(dens2$x, dens2$y, from = -1, to = 1,miny=0,
            col = alpha("maroon1", .4), border = "black")
  lines(dens3$x, dens3$y,lwd=2)
  polyCurve(dens3$x, dens3$y, from = -1, to = 1,miny=0,
            col = alpha("maroon3", .7), border = "black")
  legend('topleft', c("All Genes", "KEGG Genes", "Zonated Genes"), 
                fill=c(alpha("gray", .2),alpha("maroon1", .4),alpha("maroon3", .7)), 
                bty = "n")
                
                
  allGenesToPlot <- names(c(sort(getCorr[sub3], decreasing=T)[1:8]))
  quickScatter(allGenesToPlot)



par(mar=c(1.5,2,1,5))
  image(z = matrix(breaks, ncol = 1), col = heatcols2, 
        breaks = breaks, ylim=c(0,1))
  h <- hist(dataForPlot[useOrder,], plot = FALSE, breaks = breaks)
  hx <- scale01(breaks, 0, 1)
  hy <- c(h$counts, h$counts[length(h$counts)])
  lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
      col = 'black') 
par(mar=c(4,2,1,5))     
  image(t(dataForPlot[useOrder,]), col = heatcols2, 
          breaks = breaks,
          xaxt='n', yaxt='n')   
        
  axis(4, at = round(seq(0, 1, length.out = length(useOrder)),2), labels = useOrder, 
         las=1, cex.axis=1)
  axis(1, at=c(0,1), label=c("PC", "PP"), cex.axis=2, tick = F)


}

# Which genes are significant:
sig.h <- names(which(layerStatsPvalue[,2] < .1))
sig.m <- names(which(adjusted.pvals < .1))

sig.genes <- intersect(subset(geneSet, MortenGene %in% sig.m)[,2], sig.h)
sig.genes <- (subset(geneSet, HalpernGene %in% sig.genes)[,1])

library(clusterProfiler)
sig.genes <- bitr(sig.genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")[,2]

kk <- enrichKEGG(gene  = sig.genes,
                 organism     = 'mmu',
                 pvalueCutoff = 0.1)

barplot(kk, showCategory = 12)

firstPlot <- function(MAIN, TOP) {
  dens1 <- density(na.omit(getCorr), from=-1, to=1, bw=.3)
  par(mar=c(5,5,2,1))
  plot(dens1$x, dens1$y, cex.axis=1.5, 
       cex.lab=1.5, bty = 'n', xlim=c(-1.1,1.1), lwd=2,
       cex=1, main=MAIN, ylim=c(0,(TOP*1.2)), col="black", type="l", 
       cex.main=1.3,
       xlab="Correlation", ylab="Density", xaxt='n', yaxt='n')
  axis(2, at=round(seq(0,TOP, length.out=4),1), lwd=2,cex.axis=2)
  axis(1, at=seq(-1,1, by=.5), cex.axis=1.5, lwd=2,cex.axis=2)
  polyCurve(dens1$x, dens1$y, from = -1, to = 1,miny=0,
            col = alpha("gray", .2), border = alpha("gray", .2))
  
}


grep("amino", kk$Description)
pdf("PLOTS/FORFIGS/subplot1_forSupplementalFig1.pdf", height=6, width=15)
subcat <- kk@result[14,]
catGenes <- strsplit(subcat$geneID, "/", fixed = T)[[1]]
makeCorPlots(subcat$ID, subcat$Description, catGenes, .1, .90)
dev.off()

grep("lipo", kk$Description)
grep("P450", kk$Description)
pdf("PLOTS/FORFIGS/subplot2_forSupplementalFig1.pdf", height=6, width=15)
subcat <- kk@result[5,]
catGenes <- strsplit(subcat$geneID, "/", fixed = T)[[1]]
makeCorPlots(subcat$ID, subcat$Description, catGenes, .05, .90)
dev.off()

pdf("PLOTS/FORFIGS/subplot3_forSupplementalFig1.pdf", height=6, width=15)
subcat <- kk@result[2,]
catGenes <- strsplit(subcat$geneID, "/", fixed = T)[[1]]
makeCorPlots(subcat$ID, subcat$Description, catGenes, .1, .90)
dev.off()

pdf("PLOTS/FORFIGS/subplot4_forSupplementalFig1.pdf", height=6, width=15)
subcat <- kk@result[43,]
catGenes <- strsplit(subcat$geneID, "/", fixed = T)[[1]]
makeCorPlots(subcat$ID, subcat$Description, catGenes, .1, .95)
dev.off()

