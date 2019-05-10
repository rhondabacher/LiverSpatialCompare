# For any Kegg category, check correlation and heatmap:

setwd("LiverSpatialCompare")

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
scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
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




library(scales)
library(viridis)
library(clusterProfiler)

which(rowSums(data.norm.order) > 0)

all.genes.map <- bitr(geneSet$MortenGene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

allKegg <- enrichKEGG(gene  = all.genes.map[,2],
                 organism     = 'mmu',
                 pvalueCutoff = 1)

# list of all Kegg categories we can look at:
allKK <- subset(allKegg@result, Count >= 3)
allKK <- allKK[,c(1:2,8)]
head(allKK)

allKEGG <- paste0("KEGG: ", sort(allKK$Description))

usekegg <- allKK[1,]

catGenes <- strsplit(usekegg$geneID, "/", fixed = T)

data.norm.order.rmout <- log(PushOL(data.norm.order, qt1 = 0, qt2 = 0.97)+1)

top.res.fitted <- t(apply(data.norm.order.rmout, 1, function(x) {
  FIT = smooth.spline(1:ncol(data.norm.order.rmout),x, 
                      control.spar=list(low=.6, high=.7))
  return(FIT$y)

}))

top.res.fitted <- top.res.fitted[which(rowSums(top.res.fitted) > 0),]

save(getCorr, all.genes.map, layerMeans, layerSplit, layerMeans.morten, layerStatsPvalue, adjusted.pvals,
   all.genes.map, geneSet, top.res.fitted, allKK, allKEGG, file="tempTry.RData")




save.image("RDATA/dataAndFunctions_generateKeggPlots.RData")

