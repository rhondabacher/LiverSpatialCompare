setwd("scSpatialReconstructCompare-Paper")

# Supp. Figure 1 

load("RDATA/dataReady_bothData_genesMapped.RData")

detect.gene.halpern <- apply(halpernUMI.norm[useHalpern.genes,], 1, function(x) mean(x!=0))
detect.gene.morten <- apply(data.norm[useMorten.genes,], 1, function(x) mean(x!=0))
detect.gene.droplet <- apply(droplet.norm[useDroplet.genes, ], 1, function(x) mean(x!=0))

mean.noz.halpern <- apply(halpernUMI.norm[useHalpern.genes,], 1, function(x) mean(x[x!=0]))
mean.noz.morten <- apply(data.norm[useMorten.genes,], 1, function(x) mean(x[x!=0]))
mean.noz.droplet <- apply(droplet.norm[useDroplet.genes,], 1, function(x) mean(x[x!=0]))


diff_logMean.mh <- (log2((mean.noz.morten[geneSet[,1]])) - log(mean.noz.halpern[geneSet[,2]]))
diff_detect.mh <- (detect.gene.morten[geneSet[,1]] - detect.gene.halpern[geneSet[,2]]) 

diff_logMean.md <- (log2((mean.noz.morten[geneSet[,1]])) - log2(mean.noz.droplet[geneSet[,1]]))
diff_detect.md <- (detect.gene.morten[geneSet[,1]] - detect.gene.droplet[geneSet[,1]])


## Check GC content and Gene Length differences in the subset of genes in left side of plot.

library(biomaRt)

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
get_gene_props <- getBM(c("mgi_symbol","percentage_gene_gc_content","transcript_length"), mart= mouse)
head(get_gene_props)


diffDetect.neg.mh <- names(which(diff_detect.mh < 0)); length(diffDetect.neg.mh) #halpern bigger than morten
diffDetect.pos.mh <- names(which(diff_detect.mh >= 0)); length(diffDetect.pos.mh) #halpern smaller than morten

diffDetect.neg.md <- names(which(diff_detect.md < 0)); length(diffDetect.neg.md) #droplet bigger than morten
diffDetect.pos.md <- names(which(diff_detect.md >= 0)); length(diffDetect.pos.md) #droplet smaller than morten

length(diffDetect.neg.mh) / length(diff_detect.mh)
length(diffDetect.neg.md) / length(diff_detect.md)

diffMean.neg.mh <- names(which(diff_logMean.mh < 0)); length(diffMean.neg.mh)
diffMean.pos.mh <- names(which(diff_logMean.mh >= 0)); length(diffMean.pos.mh)

diffMean.neg.md <- names(which(diff_logMean.md < 0)); length(diffMean.neg.md)
diffMean.pos.md <- names(which(diff_logMean.md >= 0)); length(diffMean.pos.md)

length(diffMean.neg.mh) / length(diff_logMean.mh)
length(diffMean.neg.md) / length(diff_logMean.md)

###

diffDetect.neg.mh <- unique(subset(get_gene_props, mgi_symbol %in% diffDetect.neg.mh)) #halpern bigger
diffDetect.pos.mh <- unique(subset(get_gene_props, mgi_symbol %in% diffDetect.pos.mh)) #morten bigger

diffDetect.neg.md <- unique(subset(get_gene_props, mgi_symbol %in% diffDetect.neg.md)) #droplet bigger
diffDetect.pos.md <- unique(subset(get_gene_props, mgi_symbol %in% diffDetect.pos.md)) #morten bigger

###

diffMean.neg.mh <- unique(subset(get_gene_props, mgi_symbol %in% diffMean.neg.mh)) #halpern bigger
diffMean.pos.mh <- unique(subset(get_gene_props, mgi_symbol %in% diffMean.pos.mh)) #morten bigger

diffMean.neg.md <- unique(subset(get_gene_props, mgi_symbol %in% diffMean.neg.md)) #droplet bigger
diffMean.pos.md <- unique(subset(get_gene_props, mgi_symbol %in% diffMean.pos.md)) #morten bigger


pdf("PLOTS/detectionRat_vs_Expr_compareDatasets_on_GC_Fig1Supp_HadleyRemake2.pdf", height=6, width=7)
par(mar=c(5,5,1,1), mgp=c(3,1,0), mfrow=c(2,2))
plot(density(diffDetect.neg.mh$percentage_gene_gc_content), ylim=c(0,.1),
     lwd=1.5, xlab= "GC Content", ylab= "Density", main="", xlim=c(20,80),
     cex.lab=1.5, cex.axis=2, pch=20, col="cadetblue3", bty='n', xaxt='n', yaxt='n')
axis(1, at=seq(20,80, by=10), cex.axis=1.3)
axis(2, at=seq(0,.1, by=.05), cex.axis=1.3)
lines(density(diffDetect.pos.mh$percentage_gene_gc_content), col="gray40", lwd=2)
lines(density(diffMean.neg.mh$percentage_gene_gc_content), col="cadetblue3", lwd=2, lty=2)
lines(density(diffMean.pos.mh$percentage_gene_gc_content), col="gray40", lwd=2, lty=2)
legend(20,.099, title="Genes with higher \n detection or mean in: ", 
       c("Smart-seq Data", "MARS-seq Data"), bty='n', xjust=0,
       lwd=2, col=c("gray40", "cadetblue3"), cex=.5)
legend(60,.1, title="Comparison of:", 
       c("Detection fractions", "Mean"), bty='n', xjust=0,
       lwd=2, lty=c(1,3), col=c("black", "black"), cex=.5)
       
       
plot(density(log(diffDetect.neg.mh$transcript_length)), ylim=c(0,.6),
     lwd=1.5, xlab= "log2 Gene Length", ylab= "Density", main="", xlim=c(3, 12),
     cex.lab=1.5, cex.axis=2, pch=20, col="cadetblue3", bty='n', xaxt='n', yaxt='n')
axis(1, at=seq(3,12, by=2), cex.axis=1.3)
axis(2, at=seq(0,.6, by=.2), cex.axis=1.3)
lines(density(log(diffDetect.pos.mh$transcript_length)), col="gray40", lwd=2)
lines(density(log(diffMean.neg.mh$transcript_length)), col="cadetblue3", lwd=2, lty=2)
lines(density(log(diffMean.pos.mh$transcript_length)), col="gray40", lwd=2, lty=2)
legend(3,.6, title="Genes with higher \n detection or mean in: ", 
       c("Smart-seq Data", "MARS-seq Data"), bty='n', xjust=0,
       lwd=2, col=c("gray40", "cadetblue3"), cex=.5)
legend(9,.62, title="Comparison of:", 
       c("Detection fractions", "Mean"), bty='n', xjust=0,
       lwd=2, lty=c(1,3), col=c("black", "black"), cex=.5)
       
       
plot(density(diffDetect.neg.md$percentage_gene_gc_content), ylim=c(0,.1),
     lwd=1.5, xlab= "GC Content", ylab= "Density", main="", xlim=c(20,80),
     cex.lab=1.5, cex.axis=2, pch=20, col="chartreuse3", bty='n', xaxt='n', yaxt='n')
axis(1, at=seq(20,80, by=10), cex.axis=1.3)
axis(2, at=seq(0,.1, by=.05), cex.axis=1.3)
lines(density(diffDetect.neg.md$percentage_gene_gc_content), col="chartreuse3", lwd=2)
lines(density(diffDetect.pos.md$percentage_gene_gc_content), col="gray40", lwd=2)
lines(density(diffMean.neg.md$percentage_gene_gc_content), col="chartreuse3", lwd=2, lty=2)
lines(density(diffMean.pos.md$percentage_gene_gc_content), col="gray40", lwd=2, lty=2)
legend(20,.099, title="Genes with higher \n detection or mean in: ", 
       c("Smart-seq Data", "Droplet Data"), bty='n', xjust=0,
       lwd=2, col=c("gray40",  "chartreuse3"), cex=.5)
legend(60,.1, title="Comparison of:", 
       c("Detection fractions", "Mean"), bty='n', xjust=0,
       lwd=2, lty=c(1,3), col=c("black", "black"), cex=.5)
       
plot(density(log(diffDetect.neg.md$transcript_length)), ylim=c(0,.6),
     lwd=1.5, xlab= "log2 Gene Length", ylab= "Density", main="", xlim=c(3, 12),
     cex.lab=1.5, cex.axis=2, pch=20, col="cadetblue3", bty='n', xaxt='n', yaxt='n')
axis(1, at=seq(3,12, by=2), cex.axis=1.3)
axis(2, at=seq(0,.6, by=.2), cex.axis=1.3)
lines(density(log(diffDetect.neg.md$transcript_length)), col="chartreuse3", lwd=2)
lines(density(log(diffDetect.pos.md$transcript_length)), col="gray40", lwd=2)
lines(density(log(diffMean.neg.md$transcript_length)), col="chartreuse3", lwd=2, lty=2)
lines(density(log(diffMean.pos.md$transcript_length)), col="gray40", lwd=2, lty=2)
legend(3,.6, title="Genes with higher \n detection or mean in: ", 
       c("Smart-seq Data", "Droplet Data"), bty='n', xjust=0,
       lwd=2, col=c("gray40","chartreuse3"), cex=.5)
legend(9,.62, title="Comparison of:", 
       c("Detection fractions", "Mean"), bty='n', xjust=0,
       lwd=2, lty=c(1,3), col=c("black", "black"), cex=.5)

dev.off()