setwd("LiverSpatialCompare")

#########
## Load and process all Halpern data, downloaded from supplement.
#########

#########
# First the layer means:
#########
library(readxl)
layerStats <- read_excel("DATA/NIHMS70855-supplement-Supplementary_Table_3.xlsx", sheet=1)
head(layerStats); dim(layerStats)

# Formatting:
layerMeans <- data.matrix(layerStats[3:27299, 2:10])
colnames(layerMeans) <- layerStats[2,][2:10]
rownames(layerMeans) <- t(layerStats[,1])[3:27299]
rownames(layerMeans) <- gsub(" ", "", rownames(layerMeans), fixed=TRUE)

layerStatsPvalue <- data.matrix(layerStats[3:27299, 20:21])
colnames(layerStatsPvalue) <- layerStats[2,][20:21]
rownames(layerStatsPvalue) <- t(layerStats[,1])[3:27299]

#########
# Next the UMI table in the supplement:
#########
# I manually fixed an extra space in this file:
halpernUMI.supp <- readr::read_table2("DATA/Table_S1.txt", col_names = T)

library(dplyr)

# Formatting:
halpernUMI <- halpernUMI.supp
hgenes <- halpernUMI.supp %>% pull('GeneSymbol')
halpernUMI <- halpernUMI[,-c(1, 1417)]

halpernUMI <- data.matrix(halpernUMI)
rownames(halpernUMI) <- hgenes[1:nrow(halpernUMI)]

dim(halpernUMI)

#########
## Data needs to be normalized: Use scran for UMI data.
#########

library(scran)

sce <- SingleCellExperiment(list(counts=halpernUMI))
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))

sce <- normalize(sce, return_log=FALSE)


halpernUMI.norm <- normcounts(sce)


####################################

save(halpernUMI, halpernUMI.norm, layerStats, layerMeans, layerStatsPvalue, 
      file="RDATA/dataReady_HalpernPaper.RData")

