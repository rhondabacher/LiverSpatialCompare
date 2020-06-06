setwd("scSpatialReconstructCompare-Paper")

## We want to compare datasets, but not all the same genes are annotated. Let's try to match them up here:

load("RDATA/dataReady_HalpernPaper.RData")
load("RDATA/dataReady_Morten.RData")
load("RDATA/dataReady_DropletPaper.RData")


# First, note that the original paper since it is UMI cannot distinguish certain genes that shared an exon.
# The authors seperated these genes using a semicolon ;

####### 
## How many unique genes?
####### 
getGeneLengths <- sapply(rownames(halpernUMI), function(x) {
  length(strsplit(x, ";")[[1]])
}
)
table(getGeneLengths) 

# Unique genes, 16975
uniGenes <- intersect(intersect(names(which(getGeneLengths == 1)), rownames(data.norm)), rownames(droplet.norm))
length(uniGenes)

# Shared exon genes, we still want to be able to compare them...
# There are 2109 number of genes that have share exons:
sum(getGeneLengths > 1)


####### 
# Create a gene map for those shared exon genes:
####### 

multiGenes <- names(which(getGeneLengths > 1))

# For those that are not unique, which of them are in Morten's data?:
getMapping <- sapply(unique(c(rownames(data.norm), rownames(droplet.norm))), function(x) {
  
  getX = grep(x, multiGenes, value=TRUE)
  
  # grep isn't an exact match so this next bit figures out if it is exact:
  if (length(getX) > 0 & any(x == unlist(strsplit(getX, ";")))) {
    
    for (j in 1:length(getX)) {
      if (any(x == strsplit(getX, ";")[[j]])) {
        Y <- getX[j]
      } else {Y <- NA}
    }
  } else { Y <- NA}
  return(Y)
}
)
head(getMapping)

# 2299 genes in our/droplet data are not distinguishable in Halpern's data:
sum(!is.na(getMapping))

# Make into a nice dataframe to use:
getMapping <- data.frame(GENE = names(getMapping[!is.na(getMapping)]), MULTIGENE = unlist(getMapping[!is.na(getMapping)]),
                         stringsAsFactors=FALSE)

geneSet <- data.frame(Gene = c(getMapping$GENE, uniGenes), 
                      HalpernGene = c(getMapping$MULTIGENE, uniGenes),stringsAsFactors=F)
head(geneSet)

# These are the -unique- genes to use in each dataset:
useMorten.genes <- unique(intersect(rownames(data.norm), c(getMapping$GENE, uniGenes)))
useHalpern.genes <- unique(c(getMapping$MULTIGENE, uniGenes))
useDroplet.genes <- unique(intersect(rownames(droplet.norm), c(getMapping$GENE, uniGenes)))


save.image("RDATA/dataReady_bothData_genesMapped.RData")

