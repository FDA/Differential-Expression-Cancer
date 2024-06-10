###############################################################  
#
# Code written by: 
#  Nathan J. Clement
#
# With the collaboration of:
#  Haim Bar
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

#Load libraries for analysis
library(foreach)
library(doParallel)

#Cancer types modeled in this study.
cancer_types <- c("KIRC", "HNSC", "LUSC", "LGG")

setwd("/home/nathan.clement/Cancer Survival/data/final_gene_table/")

#Load ENSEMBL's fasta of human genetic coding sequences.
#https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/
ensembl_fasta <- read.fasta("Homo_sapiens.GRCh38.cds.all.fa")

#Extract ENST identifiers from the fasta.
ensembl_names <- names(ensembl_fasta)

#Reload the fasta with no annotations for simpler indexing.
ensembl_fasta <- read.fasta("Homo_sapiens.GRCh38.cds.all.fa", seqonly = T)

#Initialize parallel back-end.
#Un-comment the next line, and add the number of CPUs to parallelize over.
#co <- 
cl <- makeCluster(co)
registerDoParallel(cl)

############################################################### 

for (c in cancer_types) {
  
  ############################################################### 
  
  #Strings for paths to reference later in the script.
  #Change "..." to a path relevant to your file system.
  production_wd <- paste0(".../Production_Analysis/")
  cancer_wd <- paste0(production_wd, c, "/")
  mutation_wd <- paste0(cancer_wd, "Mutations/")
  
  ############################################################### 
  
  #Load mutation metadata and sample sheet.
  setwd(cancer_wd)
  mutation_metadata <- jsonlite::read_json("mutation_metadata.json")
  mutation_sample_type <- read.csv(file = "mutation_sample_sheet.csv")
  
  #Load reference codon usage.
  load("original_gencode_gene_codon.counts")
  
  #Read significant gene results. 
  sig_gene_results <- read.csv(paste0(c, "_significant_genes.csv"))
  
  #Extract the differentially expressed genes. 
  dge_inds <- which(sig_gene_results[,which(colnames(sig_gene_results)=="dE.q.val")]<=0.05)
  dge <- sig_gene_results[dge_inds,]
  
  #Extract the differentially variable genes. 
  dgv_inds <- which(sig_gene_results[,which(colnames(sig_gene_results)=="dV.q.val")]<=0.05)
  dgv <- sig_gene_results[dgv_inds,]
  
  #Extract the significant genes. 
  sig_gene_inds <- which(sig_gene_results[,which(colnames(sig_gene_results)=="Correlated")]==1)
  sig_gene <- sig_gene_results[sig_gene_inds,]
  
  #Extract transcript ID's and exclude their version number. 
  tran_search <- foreach(e = unlist(gene_codon_counts[,3]), .combine = 'c') %do% {
    
    strsplit(e, split = ".", fixed = T)[[1]][1]
    
  }
  
  ############################################################### 
  
  setwd(mutation_wd)
  
  #The directory names are collected so that they can be used to connect a file's name to a patient's case ID. 
  mutation_patients <- list.dirs(recursive = F, full.names = F)
  
  #This loop also pairs file ID's and case ID's, but it also removes mutations that are not in the sample sheet.
  mutation_case_id <- foreach(e = 1:length(mutation_metadata), .combine = 'rbind') %do% {
    
    pfid <- mutation_patients[which(mutation_patients==mutation_metadata[[e]]$file_id)]
    cid <- mutation_metadata[[e]][["associated_entities"]][[1]][["case_id"]]
    
    sample_check <- pfid %in% mutation_sample_type[,1]
    
    if (sample_check==T) {
      c(pfid, cid)
    }
    
  }
  
  colnames(mutation_case_id) <- c("File ID", "Case ID")
  
  ############################################################### 
  
  #Extracts the relevant mutation information from all of the mutation files, and binds it into a data table.
  mutation_table <- foreach(p = mutation_patients, 
                             .combine = 'rbind', 
                             .multicombine = T,
                             .inorder = F,
                             .packages = c("doParallel", "foreach")) %dopar% {
                               
                               setwd(paste0(mutation_wd, p))
                               
                               #TCGA mutation files have the .maf file type.
                               maf_name <- list.files(pattern = "*.maf")
                               
                               #Load the .maf into a data.table. 
                               maf <- data.table::fread(file = maf_name, 
                                                        sep = "\t",
                                                        stringsAsFactors = F,
                                                        verbose = F,
                                                        data.table = T,
                                                        showProgress = T,
                                                        header = T,
                                                        fill = T,
                                                        skip = "Hugo_Symbol",
                                                        quote = "")
                               
                               #Filter non-canonical, mutations that were not in protein-coding genes, non-SNP mutations, and entries with no listed codon change.
                               canonical_inds <- which(maf$CANONICAL=="YES")
                               maf <- maf[canonical_inds,]
                               maf <- maf[which(maf$BIOTYPE=="protein_coding"),]
                               maf <- maf[which(maf$Variant_Type=="SNP"),]
                               maf <- maf[which(nchar(maf$Codons)==7),]
                               
                               #If the filtered .maf has no data, then no changes are made to loop counts.
                               #Some .maf files were also natively empty.
                               if (nrow(maf)>0) {
                                 
                                 #Only extract mutations from transcripts used in the study.
                                 good_transcripts <- maf$Transcript_ID %in% tran_search
                                 good_transcripts <- which(good_transcripts==T)
                                 maf <- maf[good_transcripts,]
                                 
                                 #Count the length of the transcript so that it can be included in the table.
                                 transcript_length <- foreach(t = maf$Transcript_ID, .combine = 'c') %do% {
                                   
                                   ensembl_ind <- grep(t, ensembl_names)
                                   nchar(ensembl_fasta[[ensembl_ind]])
                                     
                                 }
                                 
                                 #Identify differentially expressed genes.
                                 differentially_expressed <- maf$Gene %in% dge[,1]
                                 
                                 #Identify differentially variable genes.
                                 differentially_variable <- maf$Gene %in% dgv[,1]
                                 
                                 #Identify survival-associated genes.
                                 significant_for_survival <- maf$Gene %in% sig_gene[,1]
                                 
                                 #Extract relevant columns from the .maf.
                                 cbind(
                                   # "GeneSymbol",
                                   maf$SYMBOL,
                                   # "GeneID",
                                   maf$Gene,
                                   # "TranscriptID",
                                   maf$Transcript_ID,
                                   # "TranscriptLength",
                                   transcript_length,
                                   # "DifferentiallyExpressed",
                                   differentially_expressed,
                                   # "DifferentiallyVariable",
                                   differentially_variable,
                                   #"SignificantForSurvival",
                                   significant_for_survival,
                                   # "GenomicDNAChange",
                                   maf$HGVSc,
                                   #Amino Acid Change
                                   maf$HGVSp_Short,
                                   #Codon change
                                   maf$Codons,
                                   # "MutationSubtype",
                                   maf$Variant_Classification,
                                   # "ConsequenceType",
                                   maf$Consequence,
                                   # "PolyphenImpact",
                                   maf$PolyPhen,
                                   # "SiftImpact",
                                   maf$SIFT,
                                   # "VepImpact",
                                   maf$IMPACT,
                                   # Tumor UUID
                                   maf$Tumor_Sample_UUID,
                                   # Matched normal UUID
                                   maf$Matched_Norm_Sample_UUID
                                 )
                               } #End if (nrow(maf)>0)
                               
                             } #End mutation_table foreach
  
  colnames(mutation_table) <- c(
    "GeneSymbol",
    "GeneID",
    "TranscriptID",
    "TranscriptLength",
    "DifferentiallyExpressed",
    "DifferentiallyVariable",
    "SignificantForSurvival",
    "GenomicDNAChange",
    "AminoAcidChange",
    "CodonChange",
    "MutationSubtype",
    "ConsequenceType",
    "PolyphenImpact",
    "SiftImpact",
    "VepImpact",
    "TumorUUID",
    "MatchedNormalUUID")
  
  #Write the cancer mutation table to a .csv.
  setwd(cancer_wd)
  write.csv(mutation_table, file = paste0(c, "_mutation_table.csv"))    
  
} #End cancer types loop

#Stop the parallel back-end.
stopCluster(cl)
