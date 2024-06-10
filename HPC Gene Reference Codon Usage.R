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
library(seqinr)
library(foreach)
library(doParallel)

#Subsets of genes to compute transcriptomic-weighted codon usage. 
#Run this script with this canon names before the differential gene expression script is ran.
#canon_names <- c("original")

#Run this script with this canon names after the differential gene expression script is ran. 
canon_names <- c("original", "dge", "exp_threshold", "significant_genes")

#Cancer types modeled in this study.
cancer_types <- c("KIRC", "HNSC", "LUSC", "LGG")

#Initialize parallel back-end.
#Un-comment the next line, and add the number of CPUs to parallelize over.
#co <- 
cl <- makeCluster(co)
registerDoParallel(cl)

###############################################################  

for (g in cancer_types) {

  ###############################################################  
  
  #Change "..." to a path relevant to your file system.
  setwd("...")
  
  #Load an expression file for reference.
  star <- read.table(file = "d1f7678a-8106-46c9-95b7-5ec5d940f08c.rna_seq.augmented_star_gene_counts.tsv", skip = 6)
  #Extract the protein-coding genes
  star <- star[which(star[,3]=="protein_coding"),]
  #Filter genes with PAR_Y in their name.
  par_y <- grep("PAR_Y", star[,1])
  star <- star[-par_y,]
  
  ###############################################################  
  
  #Load the list of canonical ENSG and ENST IDs, and save them to reference later. 
  canonical_transcripts <- read.csv("canonical ENSG-ENST.csv")
  original_canon <- canonical_transcripts
  
  #Load ENSEMBL's fasta of human genetic coding sequences.
  #https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cds/
  ensembl_fasta <- read.fasta("Homo_sapiens.GRCh38.cds.all.fa")
  
  #Extract ENSG identifiers from the annotations in the fasta.
  ensemble_ensg <- foreach(e = 1:length(ensembl_fasta), .combine = 'c', .multicombine = T) %do% {
    
    att_test <- attributes(ensembl_fasta[[e]])
    
    enst_num <- stringr::str_extract(att_test$Annot, "(?<=ENSG)\\d+")
    paste0("ENSG", enst_num)
    
  }
  
  #Extract ENST identifiers from the fasta.
  ensembl_names <- names(ensembl_fasta)
  #Reload the fasta with no annotations for simpler indexing.
  ensembl_fasta <- read.fasta("Homo_sapiens.GRCh38.cds.all.fa", seqonly = T)
  
  #Remove the version number from the transcript IDs.
  ens_name_search <- foreach(e = ensembl_names, .combine = 'c') %do% {
    
    strsplit(e, split = ".", fixed = T)[[1]][1]
    
  }
  
  ###############################################################  
  
  #Change "..." to a path relevant to your file system.
  setwd("...")
  
  #Read significant gene results. 
  sig_gene_results <- read.csv(paste0(g, "_significant_genes.csv"))
  
  #Extract the differentially expressed genes. 
  dge_inds <- which(sig_gene_results[,which(colnames(sig_gene_results)=="dE.q.val")]<=0.05)
  dge <- sig_gene_results[dge_inds,]
  
  #Extract the differentially variable genes. 
  dgv_inds <- which(sig_gene_results[,which(colnames(sig_gene_results)=="dV.q.val")]<=0.05)
  dgv <- sig_gene_results[dgv_inds,]
  
  #Extract the significant genes. 
  sig_gene_inds <- which(sig_gene_results[,which(colnames(sig_gene_results)=="Correlated")]==1)
  sig_gene <- sig_gene_results[sig_gene_inds,]
  
  #Identify retired ENSG ID's in the STAR file.
  star_ensembl_overlap <- foreach(c = 1:nrow(star), .combine = 'c', .multicombine = T, .packages = c("foreach", "doParallel")) %dopar% {
    
    tran_search <- strsplit(star[c,1], split = ".", fixed = T)[[1]][1]
    
    tran_test <- tran_search %in% ensemble_ensg
    
    if (tran_test==T) {
      tran_search
    }
    
  }
  
  #Identify which STAR genes have a canonical transcript. 
  canonical_overlap <- foreach(c = 1:nrow(canonical_transcripts), .combine = 'c', .multicombine = T, .packages = c("foreach", "doParallel")) %dopar% {
    
    tran_search <- strsplit(canonical_transcripts[c,2], split = ".", fixed = T)[[1]][1]
    
    tran_ind <- grepl(tran_search, star[,1])
    
    T %in% tran_ind
    
  }
  
  #Filter the canonical transcripts so that only matches remain.
  canonical_transcripts <- canonical_transcripts[which(canonical_overlap==T),]
  
  #Remove the version number from the transcripts for more convenient matching.
  canonical_ensg_no_version <- foreach(c = 1:nrow(canonical_transcripts), .combine = 'c', .multicombine = T, .packages = c("foreach", "doParallel")) %dopar% {
    
    strsplit(canonical_transcripts[c,2], split = ".", fixed = T)[[1]][1]
    
  }
  
  #Make a list of each gene set's canonical transcripts ENST ID's to iterate over.
  gene_set_ensg <- list(original = canonical_ensg_no_version,
                        dge = unique(c(dge[,1], dgv[,1])),
                        exp_threshold = sig_gene_results[,1],
                        significant_genes = sig_gene[,1])
  
  #Ensure each ENST ID is associated with a single gene.
  canon_list <- foreach(i = canon_names) %:% foreach(s = 1:length(gene_set_ensg[[which(names(gene_set_ensg)==i)]]), .combine = 'rbind', .multicombine = T) %do% {
    
    gene_ensg <- gene_set_ensg[[which(names(gene_set_ensg)==i)]][s]
    gene_ind <- grep(gene_ensg, canonical_transcripts[,2])
    gene_ind_length <- length(gene_ind)
    
    #If there are multiple or zero matches return NA so these cases can be investigated further. 
    if (gene_ind_length == 1) {
      canonical_transcripts[gene_ind,]
    } else if (gene_ind_length > 1) {
      c(NA, NA, gene_ind_length)
    } else {
      c(NA, NA, NA)
    }
    
  }
  
  #If a transcript has matches with multiple genes stop the script and return an error.
  quality_check <- foreach(n = 1:length(canon_list), .combine = 'c') %do% {
    
    check <- is.na(canon_list[[n]][,2])
    T %in% check
    
  }
  
  if (any(quality_check)) {
    stop("Check list of canonical transcripts for NA values.")
  }
  
  ###############################################################  
  
  for (i in 1:length(canon_names)) {
    
    #Canonical transcripts for loop gene set. 
    canonical_transcripts <- canon_list[[i]]
    
    save(canonical_transcripts, file = paste0(canon_names[i], "_canonical.transcripts"))
    
    #Count each gene's codon usage and bind into a matrix. 
    gene_codon_counts <- foreach(s = 1:nrow(canonical_transcripts), .combine = 'rbind', .packages = c("foreach", "doParallel")) %dopar% {
      
      #Search the ENSEMBL fasta for the canonical transcript.
      tran_search <- strsplit(canonical_transcripts[s,3], split = ".", fixed = T)[[1]][1]
      tran_ind <- grepl(tran_search, ens_name_search)
      
      ens_inds <- which(tran_ind==T)
      
      #If a transcript is not found, then return NA for all codon usage values.
      if (length(ens_inds)==0) {
        
        c(canonical_transcripts[s,], rep(NA, times = 64))
        
      } else if (length(ens_inds)>1) {
        
        #Ensure errors caused by multiple matches are prevented.
        stop("More than one transcript match.")
        
      } else {
        
        #If a transcript is found, extract its sequence.
        seq <- ensembl_fasta[[ens_inds]]
        
        #Change the sequence from a string into a character vector of nucleotides.
        seq <- unlist(strsplit(seq, split = character(0)))
        
        #Change the nucleotide vector into a vector of codons. 
        seq <- foreach(w = 1:(length(seq)/3), .combine = 'c') %do% {
          
          paste0(seq[(1:3)+3*(w-1)], collapse = "")
          
        }
        
        #Count the codon usage of the sequence.
        #plyr::count() only counts codons that are observed in the sequence.
        #plyr::count() also assumes the sequence is in the correct frame.
        #These assumptions greatly reduce the computational time of codon counts.
        seq_count <- plyr::count(seq)
        
        #List all possible codons, and create an empty, named vector to add the counts to. 
        all_codons_count <- rep(0, times = 64)
        names(all_codons_count) <- toupper(seqinr::words())
        
        #Add observed counts to the empty vector. 
        for (h in 1:nrow(seq_count)) {
          ind <- which(names(all_codons_count)==seq_count[h,1])
          
          all_codons_count[ind] <- seq_count[h,2]
          
        }
        
        #Creates a vector where the first element is the transcript name, and the rest of the elements are codon counts.
        c(canonical_transcripts[s,], all_codons_count)
      }
    }
    
    #Identify which transcripts did not have a sequence.
    
    if (is.null(dim(gene_codon_counts))) {
      
      gene_codon_counts <- as.data.frame(t(as.matrix(unlist(gene_codon_counts))))
      
      foreach(b = 4:67) %do% {
        
        gene_codon_counts[,b] <- as.numeric(gene_codon_counts[,b])
        
      }
      
    } 
    
    missing_genes <- is.na(gene_codon_counts[,4])
    missing_gene_inds <- which(missing_genes==T)
    
    #If there were transcripts without a sequence, report them and write a list for manual curation.
    if (length(missing_gene_inds)>0) {
      
      #If any transcripts were missing, print them to screen and write a .csv containing their names. 
      print(paste(canon_names[i], gene_codon_counts[missing_gene_inds,1]))
      missing_canon <- canonical_transcripts[missing_gene_inds,]
      write.csv(missing_canon, file = paste0(canon_names[i], "_transcript_not_in_ensembl_human.csv"))
      
      #Write the canonical transcripts that were observed to a .csv.
      canonical_transcripts <- canonical_transcripts[which(missing_genes==F),]
      save(canonical_transcripts, file = paste0(canon_names[i], "_canonical.transcripts"))
      
    }
    
    #Write the codon usage of the observed transcripts to a .csv
    gene_codon_counts <- gene_codon_counts[which(missing_genes==F),]
    save(gene_codon_counts, file = paste0(canon_names[i], "_gencode_gene_codon.counts"))
    write.csv(gene_codon_counts, file = paste0(canon_names[i], "_gencode_gene_codon_counts.csv"))
    
  }
  
} #End cancer types loop.

#Stop the parallel back-end.
stopCluster(cl)
