setwd("./Amp-AD paper/RNAseq_normalization")
getbiocGet("GenomicRanges")
biocGet("rtracklayer")
biocGet("genoset")

library(GenomicRanges)
library(rtracklayer)
library(genoset)

ROWCNTfile = './Mayo_TCX_all_counts_matrix.txt' #row counts temporal cortex Mayo
COVfile = './MAYO_CBE_TCX_Covariates.tsv' #covariates file - PMI is post mortem int.
mayoTCXRowCountsrowCounts = read.table(ROWCNTfile) 
mayoTCXCOV = read.table(COVfile)

GTFfile = "./gencode.v24.annotation.gtf"
GTF <- import.gff(con=GTFfile , format="gtf", genome="GRCh38.p5", feature.type="exon")

#Load the annotation and reduce it to get gene length
#GTF <- import.gff(GTFfile, format="gtf", genome="GRCm38.71", asRangedData=F, feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))
elementMetadata(reducedGTF)$widths <- width(reducedGTF)
calc_length <- function(x) {
  sum(elementMetadata(x)$widths)
}
output <- sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_length)
colnames(output) <- c("Length")


#### DESeq2 package includes a variance-stabilizing transformation method (DESeq-vst) that transforms bulk RNA-seq data toward homoscedasticit

#Add the GC numbers
#elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length

write.table(output, file="Gene_lengths.tsv", sep="\t", row.names = TRUE)

