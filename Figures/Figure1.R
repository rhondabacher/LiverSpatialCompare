setwd("LiverSpatialCompare")

# Make Figure 1 

load("RDATA/dataReady_bothData_genesMapped.RData")

########
# Some basic summaries (detection rate cell and gene, mean of gene with and without zeros.)
########

detect.halpern <- apply(halpernUMI.norm[useHalpern.genes,], 2, function(x) mean(x!=0))
detect.morten <- apply(data.norm[useMorten.genes,], 2, function(x) mean(x!=0))
detect.droplet <- apply(droplet.norm[useDroplet.genes,], 2, function(x) mean(x!=0))

detect.gene.halpern <- apply(halpernUMI.norm[useHalpern.genes,], 1, function(x) mean(x!=0))
detect.gene.morten <- apply(data.norm[useMorten.genes,], 1, function(x) mean(x!=0))
detect.gene.droplet <- apply(droplet.norm[useDroplet.genes,], 1, function(x) mean(x!=0))

mean.halpern <- apply(halpernUMI.norm[useHalpern.genes,], 1, function(x) mean(x))
mean.morten <- apply(data.norm[useMorten.genes,], 1, function(x) mean(x))
mean.droplet <- apply(droplet.norm[useDroplet.genes,], 1, function(x) mean(x))

mean.noz.halpern <- apply(halpernUMI.norm[useHalpern.genes,], 1, function(x) mean(x[x!=0]))
mean.noz.morten <- apply(data.norm[useMorten.genes,], 1, function(x) mean(x[x!=0]))
mean.noz.droplet <- apply(droplet.norm[useDroplet.genes,], 1, function(x) mean(x[x!=0]))


summary(detect.halpern)
mean(detect.halpern)
mean(detect.halpern) * length(useHalpern.genes)

summary(detect.morten)
mean(detect.morten)
mean(detect.morten) * length(useMorten.genes)

summary(detect.droplet)
mean(detect.droplet)
mean(detect.droplet) * length(useDroplet.genes)

library(yarrr)

pdf("PLOTS/detectionRate_compareData_Fig1.pdf", height=6, width=7)
toPlot <- data.frame(Detection = c(detect.halpern, detect.morten, detect.droplet), 
      Data = c(rep("MARS-seq",ncol(halpernUMI.norm)), 
               rep("Smart-seq", ncol(data.norm)),
               rep("10X", ncol(droplet.norm))))
toPlot$Data <- factor(toPlot$Data,  levels = c("Smart-seq", "MARS-seq", "10X"))

par(mar=c(3,5,1,1), mgp=c(3,1,0))
yarrr::pirateplot(formula = Detection ~ Data, 
                  data = toPlot, ylim=c(0,.5), point.o = .6, inf.f.o=.5,
                  bean.f.o = .7,
                  main = "", yaxt='n', inf.method='iqr',
                  xlab = "", pal=c("#fc8d59", "#91bfdb", "#73f997"),
                  ylab = "Detection Fraction (per cell)", cex.lab=2, cex.axis=2, cex.names=2)
axis(2, at=seq(0,.6, by=.1), cex.axis=2)
dev.off()

## May put this in the supplement instead of main paper...

# Full-length vs UMIonly
diff_logMean <- (log2((mean.noz.morten[geneSet[,1]])) - log2(mean.noz.halpern[geneSet[,2]]))
diff_detect <- (detect.gene.morten[geneSet[,1]] - detect.gene.halpern[geneSet[,2]]) 

FIT <- lm(diff_logMean ~ diff_detect)
summary(FIT)

# Intercept is 3.4 on log2 scale
2^(coef(summary(FIT))[1,1])
# =~ 10
# C1 data has on average expression 10 times larger than UMI.

alldata <- data.frame(diff_detect=diff_detect, diff_logMean=diff_logMean, stringsAsFactors=F)
alldata <- alldata[which(!is.na(alldata$diff_logMean)),]
library(ggplot2)
qq <- ggplot(alldata, aes(x = diff_detect, y = diff_logMean)) + 
  geom_point(color='gray40', size=.1) +
  stat_smooth(method="lm", size=1, color="black")  +
  labs(title="", y="Log2 Fold Change", x="Difference in detection fraction (per gene)") + 
      theme_bw()  + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.title.x=element_text(size=22),  # X axis title
        axis.title.y=element_text(size=22),  # Y axis title
        axis.text.x=element_text(size=18, color="black"),  # X axis text
        axis.text.y=element_text(size=18, color="black"))+
  scale_x_continuous(breaks=seq(-1,1, by=.5), labels=seq(-1,1, by=.5), limits=c(-1, 1)) +
  scale_y_continuous(breaks=seq(-3,15, by=2), labels=seq(-3,15, by=2), limits=c(-3, 15))+ 
  geom_hline(yintercept=0, linetype="dashed", 
                color = "gray30", size=.5)+ 
  geom_vline(xintercept=0, linetype="dashed", 
                color = "gray30", size=.3)

library(grid)
grob <- grobTree(textGrob("Intercept = 3.1 \n Slope = 3.4", x=0.1, y=0.85, hjust=0,
  gp=gpar(col="black", fontsize=14, fontface="bold")))

qq <- qq + annotation_custom(grob) 

pdf("PLOTS/detectionRat_vs_Expr_compareDatasets_MvsH_Fig1.pdf", height=6, width=6)
ggExtra::ggMarginal(qq, type = "histogram", color='black', size=10,
      xparams = list(fill="#fc8d59"), yparams = list(fill="#91bfdb"))
dev.off()
 







diff_logMean <- (log2((mean.noz.morten[geneSet[,1]])) - log2(mean.noz.droplet[geneSet[,1]]))
diff_detect <- (detect.gene.morten[geneSet[,1]] - detect.gene.droplet[geneSet[,1]]) 

FIT <- lm(diff_logMean ~ diff_detect)
summary(FIT)

# Intercept is 3.4 on log2 scale
2^(coef(summary(FIT))[1,1])
# =~ 12
# C1 data has on average expression 10 times larger than UMI.

alldata <- data.frame(diff_detect=diff_detect, diff_logMean=diff_logMean, stringsAsFactors=F)
alldata <- alldata[which(!is.na(alldata$diff_logMean)),]
library(ggplot2)
qq <- ggplot(alldata, aes(x = diff_detect, y = diff_logMean)) + 
  geom_point(color='gray40', size=.1) +
  stat_smooth(method="lm", size=1, color="black")  +
  labs(title="", y="Log2 Fold Change", x="Difference in detection fraction (per gene)") + 
  theme_bw()  + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_text(size=22),  # X axis title
        axis.title.y=element_text(size=22),  # Y axis title
        axis.text.x=element_text(size=18, color="black"),  # X axis text
        axis.text.y=element_text(size=18, color="black"))+
  scale_x_continuous(breaks=seq(-.25,1, by=.25), labels=seq(-.25,1, by=.25), limits=c(-.25, 1)) +
  scale_y_continuous(breaks=seq(-5,15, by=5), labels=seq(-5,15, by=5), limits=c(-5, 15))+ 
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray30", size=.5)+ 
  geom_vline(xintercept=0, linetype="dashed", 
             color = "gray30", size=.3)

library(grid)
grob <- grobTree(textGrob("Intercept = 2.9 \n Slope = 3.7", x=0.1, y=0.9, hjust=0,
                          gp=gpar(col="black", fontsize=14, fontface="bold")))

qq <- qq + annotation_custom(grob) 

pdf("PLOTS/detectionRat_vs_Expr_compareDatasets_MvD_Fig1.pdf", height=6, width=6)
ggExtra::ggMarginal(qq, type = "histogram", color='black', size=10,
                    xparams = list(fill="#fc8d59"), yparams = list(fill="#73f997"))
dev.off()










diff_logMean <- (log2((mean.noz.halpern[geneSet[,1]])) - log2(mean.noz.droplet[geneSet[,1]]))
diff_detect <- (detect.gene.halpern[geneSet[,1]] - detect.gene.droplet[geneSet[,1]]) 

FIT <- lm(diff_logMean ~ diff_detect)
summary(FIT)

# Intercept is 3.4 on log2 scale
2^(coef(summary(FIT))[1,1])
# =~ 4
# C1 data has on average expression 10 times larger than UMI.

alldata <- data.frame(diff_detect=diff_detect, diff_logMean=diff_logMean, stringsAsFactors=F)
alldata <- alldata[which(!is.na(alldata$diff_logMean)),]
library(ggplot2)
qq <- ggplot(alldata, aes(x = diff_detect, y = diff_logMean)) + 
  geom_point(color='gray40', size=.1) +
  stat_smooth(method="lm", size=1, color="black")  +
  labs(title="", y="Log2 Fold Change", x="Difference in detection fraction (per gene)") + 
  theme_bw()  + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.title.x=element_text(size=22),  # X axis title
        axis.title.y=element_text(size=22),  # Y axis title
        axis.text.x=element_text(size=18, color="black"),  # X axis text
        axis.text.y=element_text(size=18, color="black"))+
  scale_x_continuous(breaks=seq(-1,1, by=.5), labels=seq(-1,1, by=.5), limits=c(-1, 1)) +
  scale_y_continuous(breaks=seq(-5,5, by=2), labels=seq(-5,5, by=2), limits=c(-5, 5))+ 
  geom_hline(yintercept=0, linetype="dashed", 
             color = "gray30", size=.5)+ 
  geom_vline(xintercept=0, linetype="dashed", 
             color = "gray30", size=.3)

library(grid)
grob <- grobTree(textGrob("Intercept = -0.1 \n Slope = 2.1", x=0.1, y=0.85, hjust=0,
                          gp=gpar(col="black", fontsize=14, fontface="bold")))

qq <- qq + annotation_custom(grob) 

pdf("PLOTS/detectionRat_vs_Expr_compareDatasets_HvD_Fig1.pdf", height=6, width=6)
ggExtra::ggMarginal(qq, type = "histogram", color='black', size=10,
                    xparams = list(fill="#91bfdb"), yparams = list(fill="#73f997"))
dev.off()


