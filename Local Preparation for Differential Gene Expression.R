###############################################################  
#
# Code written by: 
#  Haim Bar
#
# With the collaboration of:
#  Nathan J. Clement
#  Douglas Meyer
#  Anton A. Komar
#  Ryan C. Hunt
#  Michael DiCuccio
#  Chava Kimchi-Sarfaty
#
# For the manuscript:
#  Codon usage and gene expression patterns associated with five-year survival rates in cancer patients
#
############################################################### 

# Prepare the data for differential expression analysis using the DVX
# package (available for download from https://haim-bar.uconn.edu/software/dvx/)

# install Biobase:
# https://bioconductor.org/packages/release/bioc/html/Biobase.html
library("Biobase")
# BiocManager::install("qvalue")
# BiocManager::install("limma")

library("betaMix")

# choose a cancer type:
cancertype <- "LGG" 
#cancertype <- "KIRC"
#cancertype <- "LUSC"
#cancertype <- "HNSC"
datadir <- "~/Data/Nate/AllCodonSurvData/"
expdat <- paste0(datadir , cancertype,
                 "_gencode_patient_gene_expression.csv")
dat <- read.csv(expdat, header=TRUE)

# keep just the expression data in matrix format:
datnum <- as.matrix(dat[,-1], nrow=nrow(dat), ncol=ncol(dat))
rownames(datnum) <- dat[,1] # the row names are the gene names
# keep only genes for which the median expression (before taking log) is > 2^6
keep <- which(apply(dat[,-1], 1, min) > 64)
dat <- datnum[keep,]
cat(dim(dat))

# read the clinical data for the patients:
clinfile <- paste0(datadir, cancertype,
                   "gencodePatientCodonUsageandSurvivalData.csv")
clin <- read.csv(clinfile, header=TRUE)
cat(dim(clin), "\n")
clin$file_name <- gsub("-",".", clin$file_name)
clin$file_name <- gsub("^(\\d)","X\\1", clin$file_name)
# exclude patients with missing age:
exc <- which(clin$age_at_diagnosis == "'--")
dat <- dat[ ,-which(colnames(dat) %in% clin$file_name[exc])]
clin <- clin[-exc,]
# convert age to years:
clin$age_at_diagnosis <- as.numeric(clin$age_at_diagnosis)/365
clin$End <- as.numeric(clin$end)/365
# the vital status:
clin$vital_status <- as.factor(clin$vital_status)
clin$FYalive <- clin$FYcensor <- clin$FYdead <- rep(FALSE, nrow(clin))
clin$FYalive[which(clin$End >= 5)] <- TRUE
clin$FYcensor[which((clin$End < 5) & (clin$vital_status == "Alive"))] <- TRUE
clin$FYdead <- (!clin$FYalive) & (!clin$FYcensor)

tumor <- which(clin$tissue_type == "Primary Tumor")
dat <- dat[ ,which(colnames(dat) %in% clin$file_name[tumor])]
clin <- clin[tumor, ]
cat(dim(dat), dim(clin), "\n")

print(all(colnames(dat) %in% clin$file_name)) # Must be TRUE
# normalize the expression data:
exprs <- log2(1+dat)
minimalSet <- ExpressionSet(assayData=exprs)

rownames(clin) <- clin$file_name
all(rownames(clin) == colnames(exprs)) # must be TRUE
pData <- clin[,c("tissue_type", "age_at_diagnosis", "vital_status", 
                 "FYalive", "FYcensor", "FYdead")]
pData$tissue_type <- as.factor(pData$tissue_type)
pData$vital_status <- as.factor(pData$vital_status)
pData$FYalive <- as.factor(pData$FYalive)
pData$FYcensor <- as.factor(pData$FYcensor)
pData$FYdead <- as.factor(pData$FYdead)
#print(table(pData$tissue_type))
# select a subset of patients for fitting the model:
set.seed(231114)
set1 <- sort(sample(which(pData$FYdead == TRUE), 50))
set2 <- sort(sample(which(pData$FYdead == FALSE), 50))
# group 1: dead within 5 years, group 1: alive or censored
testpData <- pData[-c(set1, set2), ]
testexprs <- exprs[ ,-c(set1, set2)]
testphenoData <- new("AnnotatedDataFrame", data=testpData)
pData <- pData[c(set1, set2), ]
exprs <- exprs[ ,c(set1, set2)]
phenoData <- new("AnnotatedDataFrame", data=pData)
#pData(phenoData)
fitSet <- ExpressionSet(assayData=exprs,
                        phenoData=phenoData,
                        annotation=paste("TCGA", cancertype, "fit"))
testSet <- ExpressionSet(assayData=testexprs,
                         phenoData=testphenoData,
                         annotation=paste("TCGA", cancertype, "test"))
save(fitSet, file=paste0(datadir, "fit", cancertype, ".RData"))
save(testSet, file=paste0(datadir, "test", cancertype, ".RData"))

#############

# Create the matrix that betaMix uses:
M <- t(rbind(as.numeric(pData$FYdead) - 1.5, exprs))
res <- betaMix(M, maxalpha = 1/choose(ncol(M),2), ppr = 1e-3)
shortSummary(res)
plotFittedBetaMix(res)
A <- getAdjMat(res)
cat(length(which(A[1,])))

#####################
# combine results with the output from DVX (the DVX analysis is done
# via a graphical use interface.)
resvar <- read.csv(paste0(datadir,"var", cancertype,".csv"))
resmean <- read.csv(paste0(datadir,"mean", cancertype,".csv"))
#resdvx <- read.csv("LGGminlt5Mean.csv")
cors <- (sin(res$angleMat)^2)[,1]
inc <- which(names(cors) %in% resvar[,1])
cors1 <- cors[inc]
cors2 <- as.numeric(cors1 < res$ppthr)
ressave <- cbind(resmean, resvar[, -1], cors1, cors2)
colnames(ressave)[c(1, 14, 15)] <- c("Gene", "sin2angle", "Correlated")
write.csv(ressave, file=paste0(datadir,cancertype,"minlt5Cor.csv"), row.names = FALSE)
