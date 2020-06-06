setwd("scSpatialReconstructCompare-Paper")

load("RDATA/dataReady_DropletPaper.RData")

use.data <- droplet_liver_data

# Only keep genes that are not all zero for Monocle analysis:
use.data <- use.data[names(which(rowSums(use.data) > 0)),]  

library(monocle)
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
plot_ordering_genes(use.data.monocle.cds)

############### Reduce Dim #################################
use.data.monocle.cds <- reduceDimension(use.data.monocle.cds,
                                        max_components = 2, norm_method="log")

############### Order Cells ################################
use.data.monocle.cds <- orderCells(use.data.monocle.cds, reverse = T)

## Format dataframe for plotting
pt_data <- use.data.monocle.cds@phenoData@data[,3]
names(pt_data) <- use.data.monocle.cds@phenoData@data[,1]

pt_data.droplet <- pt_data

##############Significant Genes############################

diff_test_res <- differentialGeneTest(use.data.monocle.cds, fullModelFormulaStr = "~sm.ns(Pseudotime)")
head(diff_test_res[,c("gene_short_name", "pval", "qval")])

de.droplet <- diff_test_res

save.image("RDATA/analysis_Monocle_Droplet.RData")


