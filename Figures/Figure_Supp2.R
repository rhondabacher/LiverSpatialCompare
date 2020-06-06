setwd("scSpatialReconstructCompare-Paper")

##################

load("RDATA/analysis_WaveCrest_Morten.RData")

## Run Monocle on Morten data to show cells get similarly ordered:

library(monocle)
pd <- new("AnnotatedDataFrame", data = data.frame(colnames(data.norm), row.names=colnames(data.norm)))
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name=rownames(data.norm), row.names=rownames(data.norm)))
monoset <- newCellDataSet(data.matrix(data.norm), phenoData = pd, featureData = fd)

monoset <- setOrderingFilter(monoset, markers)
monoset <- estimateSizeFactors(monoset)
monoset <- estimateDispersions(monoset)

monoset <- reduceDimension(monoset, max_components = 2, method = 'DDRTree', norm_method="none")
monoset <- orderCells(monoset)

monocle_time <- monoset@phenoData@data$Pseudotime
names(monocle_time) <- monoset@phenoData@data$colnames.data.norm


comp.wc.order <- 1:66
names(comp.wc.order) <- colnames(data.norm[,wc.order])

m.order <- 1:66
names(m.order) <- names(sort(monocle_time))

pdf("PLOTS/CompareCellOrder_Monocle_Morten.pdf", height=4, width=5)
par(mar=c(5,5,2,1))
plot(comp.wc.order, m.order[names(comp.wc.order)], xlab="WaveCrest order", 
        ylab="Monocle order", main="Cell ordering", pch=20, cex=1.5)
dev.off()


cor(comp.wc.order, m.order[names(comp.wc.order)])
## Correlation is .95