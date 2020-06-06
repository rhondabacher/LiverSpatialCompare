setwd("scSpatialReconstructCompare-Paper")


load("RDATA/analysis_Monocle_Droplet.RData")


###################

dim(droplet.norm)

layerSplit <- split(1:length(pt_data.droplet), cut(seq_along(1:length(pt_data.droplet)), 9, labels = FALSE))
layerSplit <- do.call(c,sapply(1:9, function(x) rep(x, length(layerSplit[[x]]))))
# Use mean because so many zeros.
data.norm.order.droplet <- droplet.norm[,order(pt_data.droplet)]
layerMeans.droplet <- t(apply(data.norm.order.droplet, 1, function(x) {
  return(tapply(x, layerSplit, mean))
}))

getMeans.d <- layerMeans.droplet
getMeans.d <- t(apply(getMeans.d, 1, function(x) {
  ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
getMeans.d[is.na(getMeans.d)] <- 0

library(monocle)
set.seed(1551)

# Monocle sets a seed in a function which causes the same answer each time. not ideal but setting a pre-seed to avoid.
preSeed <- list(sample(1:100000, 25), sample(1:100000, 25), sample(1:100000, 25), sample(1:100000, 25), sample(1:100000, 25))

tots <- c(600,500,400,300,200) #100 didn't work. failed.
MSE <- list()
for(j in 1:length(tots)) {
  TOTAL = tots[j]
  allMSE <- c()
  for(i in 1:25) {
   samplecells <- sample(1:ncol(droplet.norm), TOTAL)
   data.norm.sim <- droplet.norm[,samplecells]

   useg <-  sample(rownames(layerMeans.droplet), 500)

   use.data <- data.norm.sim[names(which(rowSums(data.norm.sim) > 0)),]  

    pd <- new("AnnotatedDataFrame", 
              data = data.frame(colnames(use.data),
                                row.names=colnames(use.data)))
    fd <- new("AnnotatedDataFrame", 
              data = data.frame(gene_short_name=rownames(use.data), 
                                row.names=rownames(use.data)))
    use.data.monocle.cds <- newCellDataSet(as.matrix(use.data), 
                                           phenoData = pd, 
                                           featureData = fd,
                                           expressionFamily=negbinomial.size())
    # Set-up for trajectory analysis:
    use.data.monocle.cds <- estimateSizeFactors(use.data.monocle.cds)
    use.data.monocle.cds <- estimateDispersions(use.data.monocle.cds)
    disp_table <- dispersionTable(use.data.monocle.cds)
    # Choose cutoff for gene inclusion:
    use.hvg <- rownames(hvg[rev(order(hvg[,3])),])[1:200]
    use.data.monocle.cds <- setOrderingFilter(use.data.monocle.cds, use.hvg)
    
    ############### Reduce Dim #################################
    use.data.monocle.cds <- reduceDimension(use.data.monocle.cds,
                                            max_components = 2, norm_method="log")
    ############### Order Cells ################################
    use.data.monocle.cds <- orderCells(use.data.monocle.cds, reverse = F)
    
    ## Format dataframe for plotting
    pt_data <- use.data.monocle.cds@phenoData@data[,3]
    names(pt_data) <- use.data.monocle.cds@phenoData@data[,1]

    layerSplit.sim <- split(1:ncol(data.norm.sim), cut(seq_along(1:ncol(data.norm.sim)), 9, labels = FALSE))
    layerSplit.sim <- do.call(c,lapply(1:9, function(x) rep(x, length(layerSplit.sim[[x]]))))

    data.norm.order.sim <- data.norm.sim[useg,names(sort(pt_data))] 
    getMeans.s <- t(apply(data.norm.order.sim, 1, function(x) {
      return(tapply(x, layerSplit.sim, mean))
    }))
    getMeans.s <- t(apply(getMeans.s, 1, function(x) {
      ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
    getMeans.s[is.na(getMeans.s)] <- 0
   
    m1 <- mean(rowSums((getMeans.d[useg, ] - getMeans.s[useg, ])^2))
    
    data.norm.order.sim <- data.norm.sim[useg,names(rev(sort(pt_data)))] 
    getMeans.s <- t(apply(data.norm.order.sim, 1, function(x) {
      return(tapply(x, layerSplit.sim, mean))
    }))
    getMeans.s <- t(apply(getMeans.s, 1, function(x) {
      ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
    getMeans.s[is.na(getMeans.s)] <- 0
   
    
    m2 <- mean(rowSums((getMeans.d[useg, ] - getMeans.s[useg, ])^2))
    
    allMSE <- c(allMSE, pmin(m1, m2))
  print(i)    
 
  set.seed(preSeed[[j]][i])
  }
  print(j)
  MSE[[j]] <- allMSE
}

save(MSE, file="RDATA/subsampleDroplet_Cells-MSE.RData")

library(ggsci)
library("scales")

useCols.Main <- pal_npg("nrc", alpha = .5)(1)
useCols.Dots <- rep(useCols.Main, length(do.call(c,MSE)))
allMSE <- do.call(c,MSE)
allTot <- rep(1:length(MSE), each=length(MSE[[1]]))


pdf("PLOTS/MSE_subSample_sampledGenes_Droplet.pdf", height=4.5, width=7)
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


set.seed(41112)
preSeed <- list(sample(1:100000, 25), sample(1:100000, 25), 
sample(1:100000, 25), sample(1:100000, 25), sample(1:100000, 25))

DS <- c(.1, .25, .5, .75, 1)
MSE <- list()
for(j in 1:length(DS)) {
  TOTAL = DS[j]
  allMSE <- c()
  for(R in 1:25) {
    useg <-  sample(rownames(layerMeans.droplet), 500)
    
    simCounts <- list()
    for(i in 1:ncol(droplet.norm)) {
      counts = droplet.norm[,i]
      geneProbs_all <- log(counts / sum(counts))

      downSample <- TOTAL*sum(counts)

      simCounts[[i]] <- rmultinom(n=1, size=downSample, prob=exp(geneProbs_all))
    }
    simCounts <- do.call(cbind, simCounts)
    simCounts <- matrix(simCounts, nrow=nrow(droplet.norm), ncol=ncol(droplet.norm), byrow=FALSE)
    rownames(simCounts) <- rownames(droplet.norm)
    colnames(simCounts) <- colnames(droplet.norm)
    
    use.data <- simCounts[names(which(rowSums(simCounts) > 0)),]  

    pd <- new("AnnotatedDataFrame", 
              data = data.frame(colnames(use.data),
                                row.names=colnames(use.data)))
    fd <- new("AnnotatedDataFrame", 
              data = data.frame(gene_short_name=rownames(use.data), 
                                row.names=rownames(use.data)))
    use.data.monocle.cds <- newCellDataSet(as.matrix(use.data), 
                                           phenoData = pd, 
                                           featureData = fd,
                                           expressionFamily=negbinomial.size())
    # Set-up for trajectory analysis:
    use.data.monocle.cds <- estimateSizeFactors(use.data.monocle.cds)
    use.data.monocle.cds <- estimateDispersions(use.data.monocle.cds)
    disp_table <- dispersionTable(use.data.monocle.cds)
    # Choose cutoff for gene inclusion:
    use.hvg <- rownames(hvg[rev(order(hvg[,3])),])[1:200]
    use.data.monocle.cds <- setOrderingFilter(use.data.monocle.cds, use.hvg)
    ############### Reduce Dim #################################
    use.data.monocle.cds <- reduceDimension(use.data.monocle.cds,
                                            max_components = 2, norm_method="log")
    ############### Order Cells ################################
    use.data.monocle.cds <- orderCells(use.data.monocle.cds, reverse = F)
    
    ## Format dataframe for plotting
    pt_data <- use.data.monocle.cds@phenoData@data[,3]
    names(pt_data) <- use.data.monocle.cds@phenoData@data[,1]
    
    layerSplit.sim <- split(1:ncol(simCounts), cut(seq_along(1:ncol(simCounts)), 9, labels = FALSE))
    layerSplit.sim <- do.call(c,lapply(1:9, function(x) rep(x, length(layerSplit.sim[[x]]))))

    data.norm.order.sim <- simCounts[useg,names(sort(pt_data))] 
    getMeans.s <- t(apply(data.norm.order.sim, 1, function(x) {
      return(tapply(x, layerSplit.sim, mean))
    }))
    getMeans.s <- t(apply(getMeans.s, 1, function(x) {
      ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
    getMeans.s[is.na(getMeans.s)] <- 0
    
    m1 <- mean(rowSums((getMeans.d[useg, ] - getMeans.s[useg,])^2))
    
    
    data.norm.order.sim <- simCounts[useg,names(rev(sort(pt_data)))] 
    getMeans.s <- t(apply(data.norm.order.sim, 1, function(x) {
      return(tapply(x, layerSplit.sim, mean))
    }))
    getMeans.s <- t(apply(getMeans.s, 1, function(x) {
      ((1 - 0)/(max(x) - min(x)))*(x - min(x)) + 0}))
    getMeans.s[is.na(getMeans.s)] <- 0
    
    m2 <- mean(rowSums((getMeans.d[useg, ] - getMeans.s[useg,])^2))
    
    allMSE <- c(allMSE, pmin(m1, m2))
  print(R)    
  set.seed(preSeed[[j]][[R]])
  }
  print(j)
  MSE[[j]] <- allMSE
}

save(MSE, file="RDATA/subsampleDroplet_Depth-MSE.RData")


library(ggsci)
library("scales")

useCols.Main <- pal_npg("nrc", alpha = .5)(1)
useCols.Dots <- rep(useCols.Main, length(do.call(c,MSE)))
allMSE <- do.call(c,MSE)
allTot <- rep(1:length(MSE), each=length(MSE[[1]]))


pdf("PLOTS/MSE_subSample_Depth_Droplet.pdf", height=4.5, width=7, useDingbats=F)
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


mean(colSums(droplet.norm))