setwd("./Amp-AD paper/RNAseq_normalization")

#install & load required packages
biocGet("GenomicRanges")
biocGet("rtracklayer")
biocGet("genoset")
biocGet("MASS")

library(GenomicRanges)
library(rtracklayer)
library(genoset)
library(MASS)


#input files - counts and covariates
RAWCNTfile = './Mayo_TCX_all_counts_matrix.txt' #row counts temporal cortex Mayo, downloaded from: https://www.synapse.org/#!Synapse:syn12104376 
COVfile = './MAYO_CBE_TCX_Covariates.tsv' #covariates file - PMI is post mortem int, downloaded from: https://www.synapse.org/#!Synapse:syn8466814



#load input files
mayoTCXRawCounts = read.table(RAWCNTfile, stringsAsFactors = FALSE, header = TRUE, row.names=1)
mayoTCXRawCounts <- mayoTCXRawCounts[5:nrow(mayoTCXRawCounts),]
#mayoTCXRawCounts = head(mayoTCXRawCounts, n=100)
#feature_depth = colSums(mayoTCXRawCounts)

mayoTCXCOV <- read.table(COVfile, stringsAsFactors = FALSE,  header = TRUE, row.names=1)
#xtabs(~ mayoTCXCOV[,'Sex'] + mayoTCXCOV[,'Tissue.Diagnosis.Sex'])



#function to analytically calculate the output for lqs similarly to lm, from this we will extract the p-values for coefficients
beta.t.tests <- function(covars,residuals,coeff) {
  n <- nrow(covars)
  p <- 1 + ncol(covars)
  std.err.resid <- sqrt(sum(residuals*residuals)/(n-p))
  for (i in c(1:ncol(covars))) {
    if ('numeric' != class(covars[,i])) {
      covars[,i] <- remap.values(covars[,i])
    }
  }
  X <- as.matrix(cbind(Intercept = rep(1,nrow(covars)),covars))
  CovMatrix <- solve(crossprod(X))
  std.err.betas <- std.err.resid * sqrt(diag(CovMatrix))
  t <- coeff/std.err.betas
  list(
    'Estimate'  = coeff,
    'Std. Error' = std.err.betas,
    't value'    = t,
    'Pr(>|t|)'   = 2.0*pt(-abs(t), df=n-p) # lower.tail=TRUE, the default
  )
}




#extract gene types for filtering
GTFfile = "./gencode.v24.annotation.gtf"   # gene annotation file, downloaded from: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz
GTF <- import.gff(con=GTFfile , format="gtf", genome="GRCh38.p5", feature.type="exon")
#GTF <- head(GTF, n=100)

#extract all possible gene types
miRNAgene_types = unique(elementMetadata(GTF)$gene_type)
#create a dataframe with geneid - gene_type to filter genes if necessary
gene_types <- data.frame(list(
  gene_id = elementMetadata(GTF)$gene_id,
  gene_type = elementMetadata(GTF)$gene_type),stringsAsFactors = FALSE)
gene_types <- gene_types[match(unique(gene_types$gene_id),gene_types$gene_id),]

#save gen_types 
write.table(gene_types, file="gene_types_id.tsv", sep="\t", row.names = TRUE)


gt2 <- function(x) { length(which(x > 1)) }

#Filter according to gene type (extracted from genecode), In this example - only protein_coding
PCGs <- gene_types$gene_id[gene_types$gene_type %in% c('protein_coding')]
mayoTCXRawCounts <- mayoTCXRawCounts[rownames(mayoTCXRawCounts) %in% PCGs,]

#Remove genes that have less than 1 cpm counts in at least 50% of samples
df <- t(mayoTCXRawCounts) # Rows are samples, columns are genes
df <- df[,2*apply(df,2,gt2) > nrow(df)]
df <- 1e6*sweep(df,2,colSums(df),'/')

#add covatiate columns to data matrix: Sex, AgeAtDeath, PMI
dfCOV <- mayoTCXCOV[,c('Sex','AgeAtDeath','PMI')]
j <- match(paste('X',rownames(dfCOV),sep=''),rownames(df))
j <- j[! is.na(j)]
df <- as.data.frame(cbind(df,matrix(0,nrow=nrow(df),ncol=ncol(dfCOV))))
colnames(df)[(ncol(df)-ncol(dfCOV)) + c(1:ncol(dfCOV))] <- colnames(dfCOV)
df[j,colnames(dfCOV)] <- dfCOV[j,]

#clean samples without covariate information 
#(improve, currently I'm putting zeros and removing those that don't have sex info)
df <- df[df[,'Sex'] != 0,]
class(df)
dim(df)
df0 <- df

#calculate LQS models for all genes using the covariates
system.time (
LL <- lapply(c(1:(ncol(df)-ncol(dfCOV))),  ## for testing: lapply(c(1:1000),
             function (i) {
               lqs(as.formula(paste(colnames(df)[i],'~',paste(colnames(df)[ncol(df)-2:0],collapse =' + '))),data=df)
             })
)

# Compute significance of each covariate for each gene
BB <- lapply(LL,function(L) {beta.t.tests(L$model[,-1:0],L$residuals,L$coefficients)})

saveRDS(LL,file="listOfCovariateModels.rds") # Saves the models




#### convert categorical variables to 0,1,2.... FEMALE=0, MALE=1, need to improve
remap.values <- function(x) {
  as.numeric(as.factor(x)) - 1
}


#go over all the genes, (i=gene, j=covariates) 
#   for each - calculate the substraction for all significant covariates 
corrected <- df
# You could have Substructions be just the total effects
Substructions <- lapply(c(1:(ncol(df)-ncol(dfCOV))),  ## for testing: lapply(c(1:1000),
             function (i) { #for gene i
               total.effect <- rep(0,nrow(df)) # sums the covariate effects to subtract
               
               # go over all covariates for this gene, and add to the toal effect if significant
               for (j in c(1:length(names(LL[[i]]$coefficients[-1:0])))) {  #for each effect, don't look at the intercept
                cov.name = names(LL[[i]]$coefficients[-1:0])[j]
                 if (BB[[i]][[4]][cov.name] < 0.05) { #if effect is significant
                   print(paste(i,cov.name))
                   if (cov.name == 'SexMALE') {  #remap sex to 0/1
                     cov.effect <- (remap.values(LL[[i]]$model[,names(LL[[i]]$coefficients) == 'SexMALE']) *
                                   LL[[i]]$coefficients[cov.name])
                  } else {
                    cov.effect <- (LL[[i]]$model[,cov.name] *
                                   LL[[i]]$coefficients[cov.name])
                  }
                  total.effect <- total.effect + cov.effect #total effect for all covatiates in this model
                 }
               }
               total.effect <- total.effect - mean(total.effect)
               total.effect
             })


original.data <- do.call(rbind, df[,1:(ncol(df)-ncol(dfCOV))])
corrections = do.call(rbind, Substructions)
corrected_data <- original.data - corrections # We can easily compare
corrected_data[corrected_data < 0] <- 0
rownames(corrected_data) = colnames(df)[1:(ncol(df)-ncol(dfCOV))]
colnames(corrected_data) = rownames(df)

write.table(corrected_data, file="corrected_data.tsv", sep="\t", row.names = TRUE)



##Code to extract the gene length from the GTF file, in case we want to compute FPKM

## *****not currently used************
# #Load the annotation and reduce it to get gene length
# GTF <- import.gff(GTFfile, format="gtf", genome="GRCm38.71", asRangedData=F, feature.type="exon")
# grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
# reducedGTF <- unlist(grl, use.names=T)
# elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))
# elementMetadata(reducedGTF)$widths <- width(reducedGTF)
# calc_length <- function(x) {
#   sum(elementMetadata(x)$widths)
# }
# output <- sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_length)
# colnames(output) <- c("Length")
# 
# write.table(output, file="Gene_lengths.tsv", sep="\t", row.names = TRUE)


#notes:

#filtering:
#Remove genes that have less than 1 cpm counts in at least 50% of samples per Tissue x Diagnosis
#Remove genes with missing gene length and percentage GC content - not implemented yet

#gene types proportions:
#antisense	0.0464993
#lincRNA	0.0449709
#processed_pseudogene	0.0109929
#protein_coding	0.8436306
#TEC	0.0100523
#transcribed_unprocessed_pseudogene	0.0121098

#nonuseful code:
# system.time (
#   LLM <- lapply(c(1:1000),#(ncol(df)-ncol(dfCOV))),
#                 function (i) {
#                   lm(as.formula(paste(colnames(df)[i],'~',paste(colnames(df)[ncol(df)-2:0],collapse =' + '))),data=df)
#                 })
# )



# BB <- lapply(LL,function(L) {beta.t.tests(L$model[,-1:0],L$residuals,L$coefficients)})
# zz <- sapply(c(1:1000),function(k){BB[[k]][[4]]['PMI']})
# length(which(zz < 0.05/1000))
# p.Sex <- sapply(c(1:1000),function(k){BB[[k]][[4]]['SexMALE']})
# 
# Substructions <- lapply(c(1:(ncol(df)-ncol(dfCOV))),  ## for testing: lapply(c(1:1000),
#                         function (i) {
#                           total.effect <- rep(0,samples) # sums the effects to substruct from all covariates
#                           # go over all covariates for this gene, and add to the toal effect if significant
#                           for (j in ncol(dfCOV)) {  #for each effect
#                             if (effect is significant) {
#                               
#                               
#                               sex.effect <- remap.values(LL[[i]]$model[,names(LL[[i]]$coefficients) == 'SexMALE']) *
#                                 LL[[i]]$coefficients['SexMALE']
#                               total.effect <- total.effect + sex.effect
#                             }
#                           }
#                           total.effect <- total.effect - mean(total.effect)
#                           corrected <- df
#                           corrected[,i] <- corrected[,i] - total.effect
#                           
#                         })
