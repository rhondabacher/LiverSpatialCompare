setwd("LiverSpatialCompare")

# Comparing expression over the spatial axis

## Load in datasets:
load("RDATA/dataReady_bothData_genesMapped.RData")
load("RDATA/analysis_WaveCrest_Morten.RData")

data.norm.order <- data.norm[,wc.order]

# Calculate means for Morten's data for plotting:
layerSplit <- split(1:ncol(data.norm.order), cut(seq_along(1:ncol(data.norm.order)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))

layerMeans.morten <- t(apply(data.norm.order, 1, function(x) {
  return(tapply(x, layerSplit, mean))
}))


## Calculate the correlation of gene expression across the 9 layers:
getCorr <- c()
for(i in 1:nrow(geneSet)) {
  
  gene.m <- geneSet$MortenGene[i]
  gene.h <- geneSet$HalpernGene[i]
  
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
  plot(dens1$x, dens1$y, cex.axis=1.8, 
       cex.lab=1.8, bty = 'n', xlim=c(-1.1,1.1), lwd=2,
       cex=1, main=MAIN, ylim=c(0,TOP), col="black", type="l", cex.main=2,
       xlab="Correlation", ylab="Density", xaxt='n', yaxt='n')
  axis(2, at=round(seq(0,TOP, length.out=4),1), lwd=2,cex.axis=2)
  axis(1, at=seq(-1,1, by=.5), cex.axis=1.5, lwd=2,cex.axis=2)
  polyCurve(dens1$x, dens1$y, from = -1, to = 1,miny=0,
            col = alpha("gray", .2), border = alpha("gray", .2))
  
}

#######################################################
#######################################################

# Which genes are significant, broadly defined:
sig.h <- names(which(layerStatsPvalue[,2] < .1))
sig.m <- names(which(adjusted.pvals < .1))

# Limit to genes with at least ten non-zero values in -both- datasets.
detect.gene.halpern <- apply(halpernUMI.norm[useHalpern.genes,], 1, function(x) sum(x!=0))
detect.gene.morten <- apply(data.norm[useMorten.genes,], 1, function(x) sum(x!=0))

useg.h <- names(which(detect.gene.halpern >= 10))
useg.m <- names(which(detect.gene.morten >= 10))

useg.h <- subset(geneSet, HalpernGene %in% useg.h)$MortenGene
useg <- intersect(useg.h, useg.m)

########### RUN CODE BELOW ###########################

# Since liver, look at all metabloic pathways

library(clusterProfiler)
library(scales)
KEGG = "mmu01100"
NAME = "Metabolic Pathways"

keggCat1 <- bitr_kegg(KEGG, "Path", "ncbi-geneid", "mmu")[,2]

# Genes in PATH
kegg.genes <- bitr(keggCat1, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")[,2]
kegg.genes <- intersect(useg, kegg.genes)

sig.genes <- intersect(subset(geneSet, MortenGene %in% intersect(useg,sig.m))[,2], sig.h)
sig.genes <- subset(geneSet, HalpernGene %in% sig.genes)[,1]
any.sig <- sig.genes
sig.genes <- intersect(sig.genes, kegg.genes)

all.genes <- geneSet[,2]

median(getCorr[all.genes], na.rm=T)
median(getCorr[kegg.genes], na.rm=T)
median(getCorr[sig.genes], na.rm=T)
median(getCorr[any.sig], na.rm=T)

## Make plot now:
dens2 <- density(na.omit(getCorr[kegg.genes]), from=-1, to=1, bw=.3)
dens3 <- density(na.omit(getCorr[sig.genes]), from=-1, to=1, bw=.3)
  
TOP <- round(max(dens2$y, dens3$y),1)

pdf("PLOTS/correlation_Metabolic_Fig3.pdf", height=5, width=5)
par(mar=c(5,5,2,7))
firstPlot(NAME, TOP = TOP)
segments(0,0, 0,TOP, lty=2, lwd=1.5)
lines(dens2$x, dens2$y,lwd=2)
polyCurve(dens2$x, dens2$y, from = -1, to = 1,miny=0,
          col = alpha("maroon1", .4), border = "black")
lines(dens3$x, dens3$y,lwd=2)
polyCurve(dens3$x, dens3$y, from = -1, to = 1,miny=0,
          col = alpha("maroon3", .7), border = "black")
legend(-1.15,1.1, c("All Genes", "Metabolic Genes", "Zonated Metabolic Genes"), 
              fill=c(alpha("gray", .2),alpha("maroon1", .4),alpha("maroon3", .7)), 
              bty = "n")
dev.off()
  
  
# Make smoothed heatmap
data.norm.order.rmout <- log(data.norm.order+1)

top.res.fitted <- t(apply(log(data.norm.order.rmout[sig.genes,]+1), 1, function(x) {
  FIT = smooth.spline(1:ncol(data.norm.order.rmout),x, 
                      control.spar=list(low=.6, high=.7))
  return(FIT$y)
  
}))

useY <- t(apply(top.res.fitted, 1, function(x) tapply(x, layerSplit, mean)))
useOrder <- names(sort(apply(useY, 1, which.max)))

dataForPlot <- t(base::scale(t((top.res.fitted))))
dataForPlot <- dataForPlot[useOrder,]

dataForPlot <- PushOL(dataForPlot, qt1 = .02, qt2 = 0.98)
useOrder <- names(sort(apply(dataForPlot, 1, function(x)  rev(which(x == max(x)))[1])))

MIN = min(dataForPlot)
MAX = max(dataForPlot)


par(mar=c(5,2,2,5))
library(viridis)
heatcols <- cividis(100)
heatcols2 <- heatcols[c(1:35, seq(36,84, by=5), 85:100)]

image(t(dataForPlot[useOrder,]), col = heatcols2, 
        breaks = seq(MIN, MAX, length.out=length(heatcols2)+1),
        xaxt='n', yaxt='n')      
axis(4, at = round(seq(0, 1, length.out = length(useOrder)),2), labels = useOrder, 
       las=1, cex.axis=1)
axis(1, at=c(0,1), label=c("PC", "PP"), cex.axis=2, tick = F)
  

### Next cluster genes in the pathway group:
library(cluster)
dmat <- dist(dataForPlot)
hmat <- kmeans(dmat, 2)
gLIST <- list(names(which(hmat$cluster == 1)),
              names(which(hmat$cluster == 2)))
gLIST <- lapply(gLIST, function(x) bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")[,2])
names(gLIST) <- as.character(c(1:2))

clustk <- compareCluster(geneCluster = gLIST, fun = "enrichKEGG", organism="mmu")
clustk@compareClusterResult <- clustk@compareClusterResult[-grep("disease", clustk@compareClusterResult$Description),]
clustk@compareClusterResult <- subset(clustk@compareClusterResult, clustk@compareClusterResult$p.adjust < .001)



pdf("PLOTS/keggEnrich_Fig3.pdf", height=7, width=7)
QQ <- dotplot(clustk, showCategory=15, font.size = 12)
print(QQ)
dev.off()


useOrder <- names(sort(apply(dataForPlot[names(which(hmat$cluster == 1)),], 1, function(x)  rev(which(x == max(x)))[1])))
useOrder <- c(useOrder, names(sort(apply(dataForPlot[names(which(hmat$cluster == 2)),], 1, function(x)  rev(which(x == max(x)))[1]))))
pdf("PLOTS/metabKegg_heat_Fig3.pdf", height=11, width=8)
par(mar=c(1,1,1,4))
heatcols <- cividis(100, direction = 1)
heatcols2 <- heatcols[c(1:15, seq(16,94, by=7), 95:100)]  
MIN <- min(dataForPlot[useOrder,])
MAX <- max(dataForPlot[useOrder,])
image(t(dataForPlot[useOrder,]), col = heatcols2, 
        breaks = seq(MIN, MAX, length.out=length(heatcols2)+1),
        xaxt='n', yaxt='n')      
axis(4, at = round(seq(0, 1, length.out = length(useOrder)),4), labels = useOrder, 
       las=1, cex.axis=.3)

dev.off()
