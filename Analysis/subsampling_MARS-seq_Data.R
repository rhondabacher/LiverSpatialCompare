setwd("LiverSpatialCompare")


load("RDATA/dataReady_HalpernPaper.RData")


library(readxl)
cellassign <- read_excel("DATA/nature21065-s3.xlsx", sheet=1, col_names=T, skip=1)

# Formatting:
loc.cell <- data.matrix(cellassign[,2:10])
rownames(loc.cell) <- gsub(" ", "", c(cellassign[,1])[[1]], fixed=T)


# Scale as mentioned in the supplement:
halpernUMI.prop <- t(t(halpernUMI)/colSums(halpernUMI))
loc.cell.prop <- t(t(loc.cell)/colSums(loc.cell))

# Calculate Layer specific means per Gene:
layer.means <- halpernUMI.prop %*% loc.cell.prop

layer.means[1:4,1:3]
layerMeans[1:4,1:3]
## Matches!


#### Now we can subsample cells and recalculate means:


rownames(layerMeans) <- gsub(" ", "", rownames(layerMeans), fixed=TRUE)
layerMeans <- layerMeans[names(which(rowSums(layerMeans) > 0)),]

getMeans.h <- layerMeans
getMeans.h <- t(apply(getMeans.h, 1, function(x) {
        ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))


set.seed(5221)
tots <- c(1400,1000,800,600,400,200,100)
MSE <- list()
for(j in 1:length(tots)) {
  TOTAL = tots[j]
  allMSE <- c()
  for(i in 1:25) {
    useg <- sample(rownames(layerMeans), 500)
    samplecells <- sample(1:ncol(halpernUMI.prop), TOTAL)
    sub.halpern <- halpernUMI.prop[,samplecells]
    sub.cell.prop <- loc.cell.prop[samplecells,]

    sub.means <- sub.halpern %*% sub.cell.prop

    getMeans.s <- sub.means[useg, ]
    getMeans.s <- t(apply(getMeans.s, 1, function(x) {
            ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
    getMeans.s[is.na(getMeans.s)] <- 0

    allMSE <- c(allMSE, mean(rowSums((getMeans.h[useg, ] - getMeans.s)^2)))

  }
  MSE[[j]] <- allMSE
}

save(MSE, file="RDATA/subsampleHalpern_Cells-MSE.RData")

library(ggsci)
library("scales")

useCols.Main <- pal_npg("nrc", alpha = .5)(1)
useCols.Dots <- rep(useCols.Main, length(do.call(c,MSE)))
allMSE <- do.call(c,MSE)
allTot <- rep(1:length(MSE), each=length(MSE[[1]]))


pdf("PLOTS/MSE_subSample_sampledGenes_Halpern.pdf", height=4.5, width=8)
par(mar=c(3,4,2,1))
plot(allTot, allMSE, pch=16, col=useCols.Dots, xaxt='n', yaxt='n',
        xlab="", ylab = "", cex=.8, cex.axis=2, ylim=c(0, max(allMSE)+.1))
axis(2, cex.axis=2)
axis(1, at=c(1:length(MSE)), label=as.character(c(tots)), cex.axis=2)
useCols.Main <- pal_npg("nrc", alpha = .8)(1)
points(unique(allTot), sapply(MSE, mean), pch="-", cex=6, col=useCols.Main)
# abline(a = 0, b = (sapply(MSE, mean)[3] - sapply(MSE, mean)[2]), lwd = 2, col = alpha("black", .5))
abline(h=seq(0,round(max(allMSE)+.1), by=.2), lwd=1, lty=3, col="gray80")
dev.off()








## Subsample reads.....


set.seed(88777)
DS <- c(.1, .25, .5, .75, 1)
MSE <- list()
for(j in 1:length(DS)) {
  TOTAL = DS[j]
  allMSE <- c()
  for(R in 1:25) {
    useg <-  sample(rownames(layerMeans), 500)
    
    simCounts <- list()
    for(i in 1:ncol(halpernUMI)) {
      counts = halpernUMI[,i]
      geneProbs_all <- log(counts / sum(counts))

      downSample <- DS[j]*sum(counts)

      simCounts[[i]] <- rmultinom(n=1, size=downSample, prob=exp(geneProbs_all))
    }
    simCounts <- do.call(cbind, simCounts)
    simCounts <- matrix(simCounts, nrow=nrow(halpernUMI), ncol=ncol(halpernUMI), byrow=FALSE)
    rownames(simCounts) <- rownames(halpernUMI)
    colnames(simCounts) <- colnames(halpernUMI)
 
    # Scale as mentioned in the supplement:
    simCounts.prop <- t(t(simCounts)/colSums(simCounts))

    # Calculate Layer specific means per Gene:
    layer.means.sim <- simCounts.prop %*% loc.cell.prop

    getMeans.s <- layer.means.sim[useg, ]
    getMeans.s <- t(apply(getMeans.s, 1, function(x) {
            ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
    getMeans.s[is.na(getMeans.s)] <- 0

    allMSE <- c(allMSE, mean(rowSums(getMeans.h[useg, ] - getMeans.s)^2))

  }
  MSE[[j]] <- allMSE
}

save(MSE, file="RDATA/subsampleHalpern_Depth-MSE.RData")


library(ggsci)
library("scales")

useCols.Main <- pal_npg("nrc", alpha = .5)(1)
useCols.Dots <- rep(useCols.Main, length(do.call(c,MSE)))
allMSE <- do.call(c,MSE)
allTot <- rep(1:length(MSE), each=length(MSE[[1]]))


pdf("PLOTS/MSE_subSample_Depth_Halpern.pdf", height=4.5, width=7, useDingbats=F)
par(mar=c(3,4,2,1))
plot(allTot, allMSE, pch=16, col=useCols.Dots, xaxt='n', yaxt='n',
        xlab="", ylab = "", cex=.8, cex.axis=2, ylim=c(0, max(allMSE)))
axis(2, cex.axis=2)
axis(1, at=c(1:length(MSE)), label=as.character(c(DS)), cex.axis=2)
useCols.Main <- pal_npg("nrc", alpha = .8)(1)
points(unique(allTot), sapply(MSE, mean), pch="-", cex=6, col=useCols.Main)
# abline(a = 2.2, b = -1*(sapply(MSE, mean)[8] - sapply(MSE, mean)[]), lwd = 2, col = alpha("black", .5))
abline(h=seq(0,round(max(allMSE)+.1), by=.5), lwd=1, lty=3, col="gray80");

dev.off()


mean(colSums(halpernUMI))