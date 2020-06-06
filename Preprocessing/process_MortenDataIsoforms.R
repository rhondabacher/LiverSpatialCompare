setwd("scSpatialReconstructCompare-Paper")

# Read in isoform data:
library(readxl)
indata <- read_excel("DATA/GSE116140_Data.xlsx", sheet=2) # Isoforms is sheet 2.

# Format data:
data.iso.all <- as.data.frame(indata)
rownames(data.iso.all)<- data.iso.all[,1]
data.iso.use <- data.iso.all[,-c(1,2,152)]

spikes <- grep("ERCC-", rownames(data.iso.use), value=T)

# Remove any cell with more than 20% spike-ins
PropSpikes = colSums(data.iso.use[spikes,]) / colSums(data.iso.use)

data.use.iso.filter <- data.iso.use[,names(which(PropSpikes <= .2))]


# Normalize the data:
library(SCnorm)

data.use.iso.filter <- data.matrix(data.use.iso.filter)
Conditions <- rep("Cells", ncol(data.use.iso.filter))
pdf("PLOTS/rawIsoData_count-depth_FilterExpr0_filterCellsP20.pdf", height=5, width=7)
mkCDplot <- plotCountDepth(Data = data.use.iso.filter, Conditions = Conditions,
                           FilterExpression=0)
dev.off()
pdf("PLOTS/fit_SCnorm_FilterExpr0_filterCellsP20_isoData.pdf", height=5, width=7)
DataNorm <- SCnorm(Data = data.use.iso.filter,
                   Conditions = Conditions,
                   PrintProgressPlots = TRUE)
dev.off()

data.iso.norm <- SingleCellExperiment::normcounts(DataNorm)
pdf("PLOTS/norm_count-depth_FilterExpr0_filterCellsP20_isoData.pdf", height=5, width=7)
mkCDplot <- plotCountDepth(Data = data.use.iso.filter, 
                           NormalizedData = data.iso.norm, Conditions = Conditions,
                           FilterExpression=0)
dev.off()


save.image("RDATA/dataReady_MortenIsoforms.RData")
