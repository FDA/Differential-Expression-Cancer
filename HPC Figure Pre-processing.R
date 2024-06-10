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

#Gene sets to compute transcriptomic-weighted codon usage. 
canon_names <- c("original", "dge", "exp_threshold", "significant_genes")

###############################################################  

for (c in cancer_types) {
  
  for (i in canon_names) {
    
    ###############################################################  
    
    #Strings for paths to reference later in the script.
    #Change "..." to a path relevant to your file system.
    production_wd <- paste0(".../Production_Analysis/")
    cancer_wd <- paste0(production_wd, c, "/")
    expression_wd <- paste0(cancer_wd, "Expression/")
    
    ###############################################################  
    
    setwd(cancer_wd)
    
    #Load reference codon usage.
    load(file = paste0(i, "_gencode_gene_codon.counts"))
    
    #Format reference codon usage as a matrix.
    gene_cu <- foreach(w = 4:ncol(gene_codon_counts), .combine = 'cbind') %do% {
      
      unlist(gene_codon_counts[,w])
      
    }
    
    colnames(gene_cu) <- colnames(gene_codon_counts)[4:ncol(gene_codon_counts)]
    
    #Extract transcript ID's and exclude their version number. 
    canonical_transcripts <- foreach(w = 1:nrow(gene_codon_counts), .combine = 'c') %do% {
      
      strsplit(unlist(gene_codon_counts[w,3]), split = ".", fixed = T)[[1]][1]
      
      
    }
    
    rownames(gene_cu) <- canonical_transcripts
    
    #Load expression metadata and sample sheets.
    expression_metadata <- jsonlite::read_json("expression_metadata.json")
    expression_sample_type <- read.csv(file = "expression_sample_sheet.csv")

    #Calculate the number of codons in each gene, and save them for local analysis.
    gene_lengths <- rowSums(gene_cu)
    names(gene_lengths) <- unlist(gene_codon_counts[,2])
    save(gene_lengths, file = paste0(i,"_gene.lengths"))
    
    ###############################################################  
    
    setwd(expression_wd)
    
    #TCGA file ID's are unique to each sample, but their name is misleading.
    #TCGA file ID's correspond to the directory that a file is in.
    #The name of the file inside of the directory is not derived from the file ID. 
    #The directory names are collected so that they can be used to connect a file's name to a patient's case ID. 
    patients <- list.dirs(recursive = F, full.names = F)
    
    #Extracting file ID and case ID pairs so that file names can be mapped to a case ID downstream. 
    patient_case_id <- foreach(e = 1:length(expression_metadata), .combine = 'rbind') %do% {
      
      pfid <- patients[which(patients==expression_metadata[[e]]$file_id)]
      cid <- expression_metadata[[e]][["associated_entities"]][[1]][["case_id"]]
      
      c(pfid, cid)
    }
    
    colnames(patient_case_id) <- c("File ID", "Case ID")
    
    ###############################################################  
    
    #This section identifies the expression indices for the genes in this loop's gene set.
    #It is done before looping through all of the expression files so that the indices can be used in the loop. 
    setwd(paste0(expression_wd, patients[1]))
    
    #There is only one .tsv saved in each expression directory.
    #Since the file name and file ID (directory) do not correspond with one another, 
    #This method allows for STAR files to be loaded automatically without knowing their name.
    file_name <- list.files(pattern = "*.tsv")
    
    #Load the STAR file.
    star <- read.table(file = file_name, skip = 6)
    
    #Extract the protein-coding genes
    star <- star[which(star[,3]=="protein_coding"),]
    
    #Filter genes with PAR_Y in their name.
    par_y <- grep("PAR_Y", star[,1])
    star <- star[-par_y,]
    
    #Identify the STAR indices that correspond to the gene set.
    star_inds <- foreach(s = 1:nrow(gene_codon_counts), 
                         .combine = 'rbind', 
                         .multicombine = T) %do% {
                           
                           tran_search <- strsplit(unlist(gene_codon_counts[s,2]), split = ".", fixed = T)[[1]][1]
                           
                           tran_ind <- grep(tran_search, star[,1])
                           
                           c(tran_search, tran_ind)
                           
                         }
    
    if (class(star_inds)=="character") {
      
      #Implement the filter, order the STAR genes by the order of the gene set, and add ENSG ID's with no version number.
      star <- star[as.numeric(star_inds[2]),]
      star <- cbind(star_inds[1], star)
    } else {
      
      #Implement the filter, order the STAR genes by the order of the gene set, and add ENSG ID's with no version number.
      star <- star[as.numeric(star_inds[,2]),]
      star <- cbind(star_inds[,1], star)
    }
    
    #Sanity check to ensure the gene set and STAR genes are in the same order.
    #If they are not, stop the script and return an error.
    star_check <- foreach(s = 1:nrow(gene_codon_counts),
                          .combine = 'rbind',
                          .multicombine = T) %do% {
                            
                            tran_search <- strsplit(unlist(gene_codon_counts[s,2]), split = ".", fixed = T)[[1]][1]
                            
                            grep(tran_search, star[,1])
                            
                          }
    
    star_check <- star_check == 1:nrow(gene_codon_counts)
    if (F %in% star_check) {
      stop("STAR names do not match gene_codon_counts.")
    }
    
    ###############################################################  
    
    #Calculate RNA sequence relative frequencies for each sample. 
    patient_gene_frequencies <- foreach(p = patients, 
                          .combine = 'cbind', 
                          .multicombine = T,
                          .packages = c("doParallel", "foreach")) %do% {
                            
                            ###############################################################  
                            
                            #Same file loading logic as described in STAR filtration.
                            setwd(paste0(expression_wd, p))
                            
                            file_name <- list.files(pattern = "*.tsv")
                            
                            star <- read.table(file = file_name, skip = 6)
                            star <- star[which(star[,3]=="protein_coding"),]
                            par_y <- grep("PAR_Y", star[,1])
                            star <- star[-par_y,]
                            
                            if (class(star_inds)=="character") {
                              
                              #Implement the filter, order the STAR genes by the order of the gene set, and add ENSG ID's with no version number.
                              star <- star[as.numeric(star_inds[2]),]
                              star <- cbind(star_inds[1], star)
                            } else {
                              
                              #Implement the filter, order the STAR genes by the order of the gene set, and add ENSG ID's with no version number.
                              star <- star[as.numeric(star_inds[,2]),]
                              star <- cbind(star_inds[,1], star)
                            }
                            
                            ###############################################################  
                            
                            #Compute RNA sequence relative frequencies from filtered STAR data.
                            gene_frequencies <- unlist(star[,5])/sum(unlist(star[,5]))
                            names(gene_frequencies) <- star[,1]
                            
                            gene_frequencies
                            
                            ###############################################################  
                            
                          }
    
    colnames(patient_gene_frequencies) <- patients
    
    #Save each sample's RNA sequence relative frequencies for local analysis.
    setwd(cancer_wd)
    save(patient_gene_frequencies, file = paste0(i,"_patient_gene.freqs"))
    
    ###############################################################  
    
  } #End gene sets loop.
  
} #End cancer types loop.