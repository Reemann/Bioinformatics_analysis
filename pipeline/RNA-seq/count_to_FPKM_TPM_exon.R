# prepare read count matrix
setwd('~/maternal_loading/2.public_data/analysis_embryoExon_tissueExon')

# Import the read count matrix data into R.
exon.counts <- as.matrix(read.delim("exon_count_matrix.csv", row.names=1, quote="", stringsAsFactors=FALSE, sep = ','))
exon.counts <- exon.counts[order(row.names(exon.counts)),]

# Import feature annotations. 
# Assign feature lenght into a numeric vector.
### A numeric vector with feature lengths which can be retrieved using the 'biomaRt' package.
exon.annotations <- read.table("~/source/bySpecies/genomeVersion/ensGene/genomeVersion.ensGene.exon_len.bed", sep="\t", header=FALSE)
exon.annotations <- exon.annotations[order(exon.annotations$V4),]
featureLength <- exon.annotations[[7]]

# Process one column at a time for fpkm calculation
fpkm_exon <- function(counts, featureLength){
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  
  # Process one column at a time for fpkm calculation
  fpkm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    N <- sum(counts[,i])
    exp( log(counts[,i]) + log(1e9) - log(featureLength[i]) - log(N) )
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(fpkm) <- colnames(counts)
  rownames(fpkm) <- rownames(counts)

  return(fpkm)
}

exon.fpkm.mat <- fpkm_exon(exon.counts, featureLength)

