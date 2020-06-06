setwd("scSpatialReconstructCompare-Paper")

load("RDATA/analysis_WaveCrest_Morten.RData")


library(WaveCrest)
library(EBSeq)


#### First check subsampling depth for correlation of ordering:

original.order <- wc.order
orig.tocomp <- 1:length(original.order)
names(orig.tocomp) <- original.order

set.seed(5845)

DS <- c(.1,.2,.5,.8, 1, 1.5, 2)
all.orders.sim <- list()
for(XX in 1:length(DS)) {
  orders <- list()
for(j in 1:25){
  simCounts <- list()
  for(i in 1:ncol(data.norm)) {
    counts = data.norm[,i]
    geneProbs_all <- log(counts / sum(counts))
  
    downSample <- DS[XX]*sum(counts)
  
    simCounts[[i]] <- rmultinom(n=1, size=downSample, prob=exp(geneProbs_all))
  }
  simCounts <- do.call(cbind, simCounts)
  simCounts <- matrix(simCounts, nrow=nrow(data.norm), ncol=ncol(data.norm), byrow=FALSE)
  rownames(simCounts) <- rownames(data.norm)

  data.norm.scale <- simCounts
  data.norm.scale[which(data.norm.scale < 1)] <- 1
  data.norm.scale <- log2(data.norm.scale)
  data.norm.scale <- Rescale(data.norm.scale)

  markersUse <- intersect(rownames(data.norm.scale), markers)
  data.markers <- data.norm.scale[markersUse,]

  rand.order <- sample(1:ncol(data.norm.scale),ncol(data.norm.scale))

  a3 <- ImpTC(data.markers, rand.order, cond.num, Ndg = 1)
  b3 <- Opt2TC(data.markers, N = 100000, a3, cond.num, Ndg = 1)

  orders[[j]] <- b3[[1]]
}
all.orders.sim[[XX]] <- orders
}

save.image("RDATA/subsampleMorten.RData")



## To make plot:

load("RDATA/subsampleMorten.RData")



corrs <- list()
for(j in 1:length(DS)) {
  corrs[[j]] <- sapply(all.orders.sim[[j]], function(x) {
    
    orders.comp <- 1:length(x)
    names(orders.comp) <- x
    orders.comp <- orders.comp[names(orig.tocomp)]
    COR <- pmax(cor(orig.tocomp, orders.comp), cor(orig.tocomp, rev(orders.comp)))
    # COR <- pmin(sum((orig.tocomp-orders.comp)^2), sum((orig.tocomp-rev(orders.comp))^2))
    
    return(COR)})
}
boxplot(corrs)


library(ggsci)
library("scales")

useCols.Main <- pal_npg("nrc", alpha = .5)(1)
useCols.Dots <- rep(useCols.Main, length(do.call(c,corrs[1:5])))
allMSE <- do.call(c,corrs[1:5])
allTot <- rep(1:length(corrs[1:5]), each=length(corrs[[1]]))


pdf("PLOTS/MSE_subSampleDepth_correlationOrdering_Morten.pdf", height=4.5, width=7)
par(mar=c(3,4,2,1))
plot(allTot, allMSE, pch=16, col=useCols.Dots, xaxt='n', yaxt='n',
        xlab="", ylab = "", cex=.8, cex.axis=2, ylim=c(.965, 1.00))
axis(2, cex.axis=2)
axis(1, at=c(1:length(corrs[1:5])), label=as.character(c(DS[1:5])), cex.axis=2)
useCols.Main <- pal_npg("nrc", alpha = .8)(1)
points(unique(allTot), sapply(corrs[1:5], mean), pch="-", cex=6, col=useCols.Main)

abline(h=seq(0,round(max(allMSE)), by=.01), lwd=1, lty=3, col="gray80")
dev.off()


data.norm.order <- data.norm[,wc.order]
# Split Morten into 9 layers just like Halpern data:
layerSplit <- split(1:ncol(data.norm.order), cut(seq_along(1:ncol(data.norm.order)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))
layerMeans.morten <- t(apply(data.norm.order, 1, function(x) {
  return(tapply(x, layerSplit, median))
}))


getMeans.m <- layerMeans.morten
getMeans.m <- t(apply(getMeans.m, 1, function(x) {
  ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
getMeans.m[is.na(getMeans.m)] <- 0

set.seed(3211)

DS <- c(.1,.25,.5, .75, 1)

MSE <- list()
for(XX in 1:length(DS)) {
allMSE <- c()
for(R in 1:25) {

useg <-  sample(rownames(layerMeans.morten), 500)
  simCounts <- list()
  for(i in 1:ncol(data.norm)) {
    counts = data.norm[,i]
    geneProbs_all <- log(counts / sum(counts))

    downSample <- DS[XX]*sum(counts)

    simCounts[[i]] <- rmultinom(n=1, size=downSample, prob=exp(geneProbs_all))
  }
  simCounts <- do.call(cbind, simCounts)
  simCounts <- matrix(simCounts, nrow=nrow(data.norm), ncol=ncol(data.norm), byrow=FALSE)
  rownames(simCounts) <- rownames(data.norm)
 
    data.norm.order.sim <- simCounts[,wc.order]
  # }

  layerMeans.morten.sim <- t(apply(data.norm.order.sim, 1, function(x) {
  return(tapply(x, layerSplit, median))}))

  
  getMeans.s <- layerMeans.morten.sim[useg, ]
  getMeans.s <- t(apply(getMeans.s, 1, function(x) {
          ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
  getMeans.s[is.na(getMeans.s)] <- 0

  allMSE <- c(allMSE, mean(rowSums((getMeans.m[useg, ] - getMeans.s)^2)))

}
  MSE[[XX]] <- allMSE
}
  

save.image("RDATA/subsampleMorten_Depth-MSE.RData")
  


library(ggsci)
library("scales")

useCols.Main <- pal_npg("nrc", alpha = .5)(1)
useCols.Dots <- rep(useCols.Main, length(do.call(c,MSE)))
allMSE <- do.call(c,MSE)
allTot <- rep(1:length(MSE), each=length(MSE[[1]]))


pdf("PLOTS/MSE_subSample_sampledGenes_Depth_MSE_Morten.pdf", height=4.5, width=7)
par(mar=c(3,4,2,1))
plot(allTot, allMSE, pch=16, col=useCols.Dots, xaxt='n', yaxt='n',
        xlab="", ylab = "", cex=.8, cex.axis=2, ylim=c(0, max(allMSE)))
axis(2, cex.axis=2)
axis(1, at=c(1:length(MSE)), label=as.character(c(DS)), cex.axis=2)
useCols.Main <- pal_npg("nrc", alpha = .8)(1)
points(unique(allTot), sapply(MSE, mean), pch="-", cex=6, col=useCols.Main)

abline(h=seq(0,round(max(allMSE),1), by=.1), lwd=1, lty=3, col="gray80")
dev.off()


  
  
  
  
################################################################################################
################################################################################################
#### Subsampling number of cells for effect on zonation profiles--calculate MSE:


set.seed(5533)

DS <- c(65, 60, 55, 50, 45)
MSE <- list()
for(XX in 1:length(DS)) {
  allMSE <- c()
  for(j in 1:25){
    data.norm.sim <- data.norm[,sample(1:ncol(data.norm), DS[XX])]
    data.norm.scale <- data.norm.sim
    data.norm.scale[which(data.norm.scale < 1)] <- 1
    data.norm.scale <- log2(data.norm.scale)
    data.norm.scale <- Rescale(data.norm.scale)
    
    markersUse <- intersect(rownames(data.norm.scale), markers)
    data.markers <- data.norm.scale[markersUse,]
    
    rand.order <- sample(1:ncol(data.norm.scale),ncol(data.norm.scale))
    
    a3 <- ImpTC(data.markers, rand.order, cond.num, Ndg = 1)
    b3 <- Opt2TC(data.markers, N = 100000, a3, cond.num, Ndg = 1)
    
    sim.order <- b3[[1]]
    
    layerSplit.sim <- split(1:ncol(data.norm.sim), cut(seq_along(1:ncol(data.norm.sim)), 9, labels = FALSE))
    layerSplit.sim <- do.call(c,lapply(1:9, function(x) rep(x, length(layerSplit.sim[[x]]))))

    useg <-  sample(rownames(layerMeans.morten), 500)

    data.norm.order.sim <- data.norm.sim[,sim.order] 
    layerMeans.morten.sim <- t(apply(data.norm.order.sim, 1, function(x) {
      return(tapply(x, layerSplit.sim, median))
    }))
    getMeans.s <- layerMeans.morten.sim[useg, ]
    getMeans.s <- t(apply(getMeans.s, 1, function(x) {
      ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
    getMeans.s[is.na(getMeans.s)] <- 0
    
    m1 <- mean(rowSums((getMeans.m[useg, ] - getMeans.s)^2))
    
    
    data.norm.order.sim <- data.norm.sim[,rev(sim.order)] 
    layerMeans.morten.sim <- t(apply(data.norm.order.sim, 1, function(x) {
      return(tapply(x, layerSplit.sim, median))
    }))
    getMeans.s <- layerMeans.morten.sim[useg, ]
    getMeans.s <- t(apply(getMeans.s, 1, function(x) {
      ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
    getMeans.s[is.na(getMeans.s)] <- 0
    
    m2 <- mean(rowSums((getMeans.m[useg, ] - getMeans.s)^2))

    allMSE <- c(allMSE, pmin(m1, m2))
    
  }
  MSE[[XX]] <- allMSE
}
## Compare to totally random order MSE
allMSE <- c()
for(j in 1:25){
    
    data.norm.order.sim <- data.norm[,sample(1:66)]
  
    layerSplit.sim <- split(1:ncol(data.norm.order.sim), cut(seq_along(1:ncol(data.norm.order.sim)), 9, labels = FALSE))
    layerSplit.sim <- do.call(c,lapply(1:9, function(x) rep(x, length(layerSplit.sim[[x]]))))
    layerMeans.morten.sim <- t(apply(data.norm.order.sim, 1, function(x) {
      return(tapply(x, layerSplit.sim, median))
    }))
    useg <-  sample(rownames(layerMeans.morten), 500)
    
    getMeans.s <- layerMeans.morten.sim[useg, ]
    getMeans.s <- t(apply(getMeans.s, 1, function(x) {
      ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
    getMeans.s[is.na(getMeans.s)] <- 0
    
    allMSE <- c(allMSE, mean(rowSums((getMeans.m[useg, ] - getMeans.s)^2)))
    
}
MSE[[6]] <- allMSE



save.image("RDATA/subsampleMorten_fewerCells.RData")






library(ggsci)
library("scales")

MSE <- MSE[1:5] # Don't plot the random cells.
useCols.Main <- pal_npg("nrc", alpha = .5)(1)
useCols.Dots <- rep(useCols.Main, length(do.call(c,MSE)))
allMSE <- do.call(c,MSE)
allTot <- rep(1:length(MSE), each=length(MSE[[1]]))


pdf("PLOTS/MSE_subSampleCells_MSE_Morten.pdf", height=4.5, width=7)
par(mar=c(3,4,2,1))
plot(allTot, allMSE, pch=16, col=useCols.Dots, xaxt='n', yaxt='n',
     xlab="", ylab = "", cex=.8, cex.axis=2)
axis(2, cex.axis=2)
axis(1, at=c(1:length(MSE)), label=as.character(c(DS)), cex.axis=2)
useCols.Main <- pal_npg("nrc", alpha = .8)(1)
points(unique(allTot), sapply(MSE, mean), pch="-", cex=6, col=useCols.Main)

abline(h=seq(0,round(max(allMSE)), by=.2), lwd=1, lty=3, col="gray80")
dev.off()
