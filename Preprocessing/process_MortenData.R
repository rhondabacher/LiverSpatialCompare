setwd("scSpatialReconstructCompare-Paper")

#########
## Process Morten's full-length data, needs to be normalized and also QC
#########


library(readxl)
indata <- read_excel("DATA/GSE116140_Data.xlsx", sheet=1) # Genes is sheet 1.
data.use <- as.data.frame(indata, stringsAsFactors=F)

#Formatting:
rownames(data.use) <- data.use[,1]
data.use <- data.use[,-c(1,151)]

# Remove cells that are low quality due to too much ERCC, use threshold of 20%:
spikes <- grep("ERCC-", rownames(data.use), value=T)
propSpikes = colSums(data.use[spikes,]) / colSums(data.use)

# Use 20% as filter:
data.use.filter <- data.use[,names(which(propSpikes <= .2))]

#########
# Normalize using SCnorm:
#########
library(SCnorm)

data.use.filter <- data.matrix(data.use.filter)
Conditions <- rep("Cells", ncol(data.use.filter))

pdf("PLOTS/rawData_count-depth_FilterExpr0_filterCellsP20.pdf", height=5, width=7)
mkCDplot <- plotCountDepth(Data = data.use.filter, Conditions = Conditions,
                             FilterExpression=2)
dev.off()
pdf("PLOTS/fit_SCnorm_FilterExpr0_filterCellsP20.pdf", height=5, width=7)
DataNorm <- SCnorm(Data = data.use.filter,
                  Conditions = Conditions,
                  PrintProgressPlots = TRUE)
dev.off()
# K = 7 converged

data.norm <- SingleCellExperiment::normcounts(DataNorm)

pdf("PLOTS/norm_count-depth_FilterExpr0_filterCellsP20.pdf", height=5, width=7)
mkCDplot <- plotCountDepth(Data = data.use.filter, NormalizedData = data.norm, Conditions = Conditions,
                             FilterExpression=0)
dev.off()


save(data.use.filter, data.norm, file="RDATA/dataReady_Morten.RData")
