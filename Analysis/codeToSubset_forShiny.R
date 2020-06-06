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
## Calculate the correlation of gene expression across the 9 layers:
geneSet <- geneSet[which(geneSet$Gene %in% intersect(useDroplet.genes, useMorten.genes)),]
getCorr_MH <- c()
getCorr_MD <- c()
getCorr_HD <- c()
for(i in 1:nrow(geneSet)) {
  
  gene.m <- geneSet$Gene[i]
  gene.h <- geneSet$HalpernGene[i]
  
  # Scale data
  getMeans.h <- layerMeans[gene.h, ]
  rescale.halpern.means <- ((1 - 0)/(max(getMeans.h) - min(getMeans.h)))*(getMeans.h - min(getMeans.h)) + 0
  
  getMeans.m <- layerMeans.morten[gene.m,]
  rescale.morten.means <- ((1 - 0)/(max(getMeans.m) - min(getMeans.m)))*(getMeans.m - min(getMeans.m)) + 0
  
  getMeans.d <- layerMeans.droplet[gene.m, ]
  rescale.droplet.means <- ((1 - 0)/(max(getMeans.d) - min(getMeans.d)))*(getMeans.d - min(getMeans.d)) + 0
  
  if(any(is.na(rescale.morten.means)) | any(is.na(rescale.halpern.means))){
    corr.h.m <- NA
  } else{
    corr.h.m <- cor(as.vector(rescale.morten.means), as.vector(rescale.halpern.means), method='pearson')
  }
  if(any(is.na(rescale.halpern.means)) | any(is.na(rescale.droplet.means))){
    corr.h.d <- NA
  } else{
    corr.h.d <- cor(as.vector(rescale.halpern.means), as.vector(rescale.droplet.means), method='pearson')
  }
  if(any(is.na(rescale.morten.means)) | any(is.na(rescale.droplet.means))){
    corr.m.d <- NA
  } else{
    corr.m.d <- cor(as.vector(rescale.morten.means), as.vector(rescale.droplet.means), method='pearson')
  }
  getCorr_MH <- c(getCorr_MH, corr.h.m)
  getCorr_HD <- c(getCorr_HD, corr.h.d)
  getCorr_MD <- c(getCorr_MD, corr.m.d)
}
names(getCorr_MH) <- geneSet$Gene
names(getCorr_HD) <- geneSet$Gene
names(getCorr_MD) <- geneSet$Gene



library(scales)
library(viridis)
library(clusterProfiler)

which(rowSums(data.norm.order) > 0)

all.genes.map <- bitr(geneSet$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

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

data.norm.order.rmout <- log(PushOL(data.norm.order, qt1 = .02, qt2 = 0.98)+1)

top.res.fitted <- t(apply(data.norm.order.rmout, 1, function(x) {
  FIT = smooth.spline(1:ncol(data.norm.order.rmout),x, df=4)
  return(FIT$y)

}))

top.res.fitted <- top.res.fitted[which(rowSums(top.res.fitted) > 0),]

layerSplit <- split(1:ncol(data.norm.order), cut(seq_along(1:ncol(data.norm.order)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))


save(adjusted.pvals, all.genes.map,de.droplet, allKEGG, allKK, geneSet, layerMeans.droplet,
      getCorr_MH, getCorr_MD, getCorr_HD, layerMeans, layerMeans.morten, layerSplit, 
      layerStatsPvalue, top.res.fitted, file="Code/dataAndFunctions_generateKeggPlots.RData")

