# install the package
if(!require(devtools)) install.packages("devtools")
devtools::install_github("AAlhendi1707/countToFPKM", build_vignettes = TRUE)
library(countToFPKM)

# prepare read count matrix
setwd('~/maternal_loading/2.public_data/analysis_embryoExon_tissueExon')

# Import the read count matrix data into R.
gene.counts <- as.matrix(read.delim("gene_count_matrix.csv", row.names=1, quote="", stringsAsFactors=FALSE, sep = ','))
gene.counts <- gene.counts[order(row.names(gene.counts)),]

# Import feature annotations. 
# Assign feature lenght into a numeric vector.
### A numeric vector with feature lengths which can be retrieved using the 'biomaRt' package.
gene.annotations <- read.table("~/source/bySpecies/genomeVersion/ensGene/genomeVersion.ensGene.gene_len.bed", sep="\t", header=FALSE)
gene.annotations <- gene.annotations[order(gene.annotations$V4),]
featureLength <- gene.annotations[[7]]

# Import sample metrics. 
# Assign mean fragment length into a numeric vector.
### a numeric vector with mean fragment length which can be calculated using the 'CollectInsertSizeMetrics(Picard)' tool.
meanFragmentLength <-

# Return FPKM into a numeric matrix.
# fpkm_matrix <- fpkm(counts, featureLength, meanFragmentLength)
fpkm_matrix <- fpkm (gene.counts, featureLength, meanFragmentLength)

# Plot log10(FPKM+1) heatmap of top 30 highly variable features
fpkmheatmap(fpkm_matrix, topvar=30, showfeaturenames=TRUE, return_log = TRUE)