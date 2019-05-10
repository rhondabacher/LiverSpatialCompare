setwd("LiverSpatialCompare")

load("RDATA/analysis_WaveCrest_Morten.RData")

library(viridis)
library(gplots)


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

# Scaling for heatmap:
data.norm.scale.pushOL <- PushOL(data.norm.scale[markers,], qt1=.05, qt2=.95)


# heatcols <- viridis(n=100)
library(RColorBrewer)
initilcol<- brewer.pal(11,"Spectral")
colors <- c("darkorchid4", initilcol[1])


heatcols <- cividis(100)
heatcol2 <- heatcols[c(1:5, seq(6, 94, by=5), 95:100)]


pdf("PLOTS/wavecrest_Morten_heatMarkers_Fig1_v2.pdf", height=4, width=6)
heatmap.2(data.norm.scale.pushOL[,wc.order], Colv=F, dendrogram='row', 
          trace="none",RowSideColors=colors[mk.groups],
          col=heatcol2, cexCol =.5, cexRow=1, 
           lhei=c(1.2,3), lwid=c(1,3), key.par = list(cex=0.45, cex.axis=1))
dev.off()


