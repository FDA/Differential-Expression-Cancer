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

for (c in cancer_types) {
  
  for (i in cannon_names) {
    
    ###############################################################  
    
    #Strings for paths to reference later in the script.
    #Change "..." to a path relevant to your file system.
    production_wd <- paste0(".../Production_Analysis/")
    cancer_wd <- paste0(production_wd, c, "/")
    expression_wd <- paste0(cancer_wd, "Expression/")
    mutation_wd <- paste0(cancer_wd, "Mutations/")
    
    ###############################################################  
    
    setwd(cancer_wd)
    
    #Load reference codon usage.
    load(file = paste0(i, "_gencode_gene_codon.counts"))
    
    #Format reference codon usage as a matrix.
    gene_counts_only <- foreach(w = 4:ncol(gene_codon_counts), .combine = 'cbind') %do% {
      
      unlist(gene_codon_counts[,w])
      
    }
    
    colnames(gene_counts_only) <- colnames(gene_codon_counts)[4:ncol(gene_codon_counts)]
    
    #Extract transcript ID's and exclude their version number. 
    canonical_transcripts <- foreach(w = 1:nrow(gene_codon_counts), .combine = 'c') %do% {
      
      strsplit(unlist(gene_codon_counts[w,3]), split = ".", fixed = T)[[1]][1]
      
      
    }
    
    rownames(gene_counts_only) <- canonical_transcripts
    
    #Load mutation/expression metadata and sample sheets.
    expression_metadata <- jsonlite::read_json("expression_metadata.json")
    mutation_metadata <- jsonlite::read_json("mutation_metadata.json")
    
    expression_sample_type <- read.csv(file = "expression_sample_sheet.csv")
    mutation_sample_type <- read.csv(file = "mutation_sample_sheet.csv")
    
    #The mutation case ID's are all duplicated. 
    #This loop extracts the first copy so that it can be used downstream. 
    mutation_cases <- foreach(l = 1:length(mutation_sample_type[,6]), .combine = 'c') %do% {
      
      unlist(strsplit(mutation_sample_type[l,6], split = "\\, "))[1]
      
    }
    
    #This loop filters mutation files with ambiguous identifiers.
    #Mutations passed through this loop have unique identifiers for each entry.
    mutation_sample_type <- foreach(u = unique(expression_sample_type[,6]), .combine = 'rbind') %do% {
      
      mut_matches <- which(mutation_cases==u)
      
      if (length(mut_matches)>0) {
        loop_sample <- mutation_sample_type[which(mutation_cases==u),]
        l_test <- length(loop_sample[,7])==length(unique(loop_sample[,7]))

        if (l_test == T) {
           loop_sample
        }
      } 
      
    }
    
    #Load clinical data.
    clinical <- read.csv(paste0(cancer_wd, "Clinical/clinical.csv"), header = T)
    
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
    
    setwd(mutation_wd)
    
    #The logic in this section for mutation file ID's and case ID's the same as what was described above. 
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
    
    #This section identifies the expression indices for the genes in this loop's gene set.
    #It is done before looping through all of the expression files so that the indices can be used in the loop. 
    setwd(paste0(expression_wd, patients[1]))
    
    #There is only one .tsv saved in each expression directory.
    #Since the file name and file ID (directory) do not correspond with one another, 
    #this method allows for STAR files to be loaded automatically without knowing their name.
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
    #If they are not, stop the script and send an error.
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
    
    #Calculate transcriptomic-weighted codon usage with mutational codon changes for each patient. 
    patient_cu <- foreach(p = patients, 
                          .combine = 'rbind', 
                          .multicombine = T,
                          .inorder = F,
                          .packages = c("doParallel", "foreach")) %dopar% {
                            
                            #Same file loading logic as described in STAR filtration.
                            setwd(paste0(expression_wd, p))
                            
                            file_name <- list.files(pattern = "*.tsv")
                            
                            #Pre-filter the STAR files in the same way.
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
                            
                            #Save the gene set specific expression file so it can be combined with the others downstream.
                            save(star, file = paste0(i, "_star.rows"))
                            
                            #Calculate unstranded frequencies for each gene.
                            gene_frequencies <- unlist(star[,5])/sum(unlist(star[,5]))
                            names(gene_frequencies) <- star[,1]
                            
                            #Copy of reference codon usage. 
                            #The copy will be modified by any mutations associated with this expression file.
                            loop_counts <- gene_counts_only
                            
                            #Identify the sample type and sample ID of the expression file.
                            sample_type <- expression_sample_type[which(expression_sample_type[,2]==file_name),8]
                            sample_id <- expression_sample_type[which(expression_sample_type[,2]==file_name),7]
                            
                            #Search for a mutation file associated with the sample ID.
                            ID_match <- grep(sample_id, mutation_sample_type[,7])
                        
                            #If there is not an associated mutation file, then no changes are made to loop_counts.
                            if (length(ID_match)>0) {
                              
                              if (sample_type == "Solid Tissue Normal") {
                              #TCGA open access mutation files do not provide germline mutation data.
                              #Consequently, no changes are made to loop counts. 
                            } else {
                              
                              #Navigate to the associated mutation's directory.
                              mut_sample_dir <- mutation_sample_type[ID_match,1]
                              
                              setwd(paste0(mutation_wd, mut_sample_dir))
                              
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
                              maf <- maf[which(nchar(maf$Codons)>0),]
                              
                              #If the filtered .maf has no data, then no changes are made to loop counts.
                              #Some .maf files were also natively empty.
                              if (nrow(maf)>0) {
                                
                                #If .maf has data, then make the corresponding codon changes in loop counts.
                                foreach(r = 1:nrow(maf), .combine = 'rbind') %do% {
                                  
                                  #Codon changes are listed in a WildTypeCodon//MutatedCodon format in the .maf.
                                  #This code separates WildTypeCodon and MutatedCodon.
                                  codon_change <- unlist(strsplit(maf$Codons[r], split = "\\/"))
                                  transcript <- maf$Transcript_ID[r]
                                  
                                  #Find the transcript in loop counts.
                                  transcript_row <- which(rownames(loop_counts)==transcript)
                                  
                                  #If the transcript is not in loop counts, then it is ignored.
                                  if (length(transcript_row)>0) {
                                    
                                    #Identify the codon columns.
                                    codon_lost <- which(colnames(loop_counts)==toupper(codon_change[1]))
                                    codon_gained <- which(colnames(loop_counts)==toupper(codon_change[2]))
                                    
                                    #Remove a copy of the lost codon, and add a copy of the gained codon.
                                    loop_counts[transcript_row, codon_lost] <- loop_counts[transcript_row, codon_lost] - 1
                                    loop_counts[transcript_row, codon_gained] <- loop_counts[transcript_row, codon_gained] + 1
                                    
                                  } #End length(transcript_row)>0
                                } #End foreach
                              } #End nrow(maf)>0
                            } #End abnormal tissue
                          } #End length(ID_match)>0

                            #Multiple each gene's usage frequency by its codon usage.
                            #This step is why there was a sanity check to ensure the star genes and gene set were in the same order.
                            #This method is significantly faster than matching each gene's frequency with its counts in a loop.
                            frequency_weighed_gene_counts <- gene_frequencies*loop_counts

                            #Each codon's transcriptomic-weighted average is the sum of its usage in the gene set.
                            colSums(frequency_weighed_gene_counts)
                            
                        } #End patient_cu foreach
    
    #Combine all the gene set unstranded expression results for each patient into a matrix.
    patient_gene_expression <- foreach(p = patients, 
                                       .combine = 'cbind', 
                                       .multicombine = T,
                                       .inorder = F) %dopar% {
                                         
                                         setwd(paste0(expression_wd, p))
                                         
                                         file_name <- list.files(pattern = "*.tsv")
                                         load(paste0(i, "_star.rows"))
                                         
                                         express <- star[,5]
                                         names(express) <- star[,1]
                                         unlist(express)
                                         
                                       } #End patient_gene_expression foreach
    
    
    rownames(patient_cu) <- patients
    colnames(patient_gene_expression) <- patients
    
    #Save the results into the project directory.
    setwd(cancer_wd)
    save(patient_cu, file = paste0(i, "_patient.cu"))
    write.csv(patient_cu, file = paste0(i, "_patient_cu.csv"))
    
    #Identify which clinical samples correspond to an expression file.
    clinical_overlap <- foreach(s = clinical[,1], .combine = 'c') %do% {
      
      s %in% patient_case_id[,2]
      
    }
    
    #Remove the clinical samples that did not have an expression file.
    clinical <- clinical[which(clinical_overlap==T),]
    
    #Each patient has two entries.
    #The entries document differences in treatment that are not relevant to this project.
    #One entry for each patient is retained for downstream analysis.
    clinical <- foreach(v = unique(clinical[,1]), .combine = 'rbind') %do% {
      
      clinical[which(clinical[,1]==v)[1],]
      
    }
    
    #Patients with prior treatment are removed from the analysis. 
    prior_treatment <- clinical[,which(colnames(clinical)=="prior_treatment")]
    clinical <- clinical[which(prior_treatment=="No"),]
    
    #Combine each patient's transcriptomic-weighted codon usage with their survival data into a data table.
    patient_cuSurvival <- foreach(y = 1:nrow(clinical), 
                                  .combine = 'rbind', 
                                  .packages = c("doParallel", "foreach"), 
                                  .multicombine = T,
                                  .inorder = F) %dopar% {
                                    
                                    #Cancer type
                                    cancer_type <- c
                                    
                                    #Case ID
                                    case_id <- clinical[y,1]
                                    
                                    #Number
                                    number <- NA
                                    
                                    #Sample(s) associated with this patient.
                                    file_id <- patient_case_id[which(patient_case_id[,2]==case_id),1]
                                    
                                    #This loop extracts the results for each sample. 
                                    foreach(f = file_id, .combine = 'rbind', .multicombine = T) %do% {
                                      
                                      #Tissue type
                                      tissue_type <- expression_sample_type[which(expression_sample_type[,1]==f),8]
                                      
                                      #Codon usage of the sample.
                                      codon_usage <- patient_cu[which(rownames(patient_cu)==f),]
                                      
                                      #Age at diagnosis
                                      age_at_diagnosis <- clinical[y,which(colnames(clinical)=="age_at_diagnosis")]
                                      
                                      #Vital status
                                      vital_status <- clinical[y,which(colnames(clinical)=="vital_status")]
                                      
                                      #Days to last follow up.
                                      follow_up <- clinical[y,which(colnames(clinical)=="days_to_last_follow_up")]
                                      
                                      #Days to death.
                                      death <- clinical[y,which(colnames(clinical)=="days_to_death")]
                                      
                                      #When the study stopped.
                                      if (follow_up == "'--") {
                                        end <- death
                                      } else if (death>follow_up) {
                                        end <- death
                                      } else {
                                        end <- follow_up
                                      }
                                      
                                      c(cancer_type, 
                                        case_id, 
                                        tissue_type, 
                                        number, 
                                        codon_usage, 
                                        age_at_diagnosis, 
                                        vital_status, 
                                        follow_up, 
                                        death, 
                                        end,
                                        f)
                                      
                                    } #End file ID foreach
                                    
                                    
                                  } #End patient_cuSurvival foreach
    
    
    cusurvival_colnames <-      c("cancer_type", 
                                  "case_id", 
                                  "tissue_type", 
                                  "number", 
                                  colnames(patient_cu), 
                                  "age_at_diagnosis", 
                                  "vital_status", 
                                  "follow_up", 
                                  "death", 
                                  "end", 
                                  "file_name")
    
    colnames(patient_cuSurvival) <- cusurvival_colnames
    
    #Write codon usage and clinical information file.
    write.csv(patient_cuSurvival, file = paste(i, "gencode Patient Codon Usage and Survival Data.csv"))
    
    #Identify which patients in the expression matrix were filtered out for the codon usage and survival computation.
    filter_file_names <- foreach(h = 1:nrow(patient_cuSurvival), .combine = 'c') %do% {
      
      which(patient_cuSurvival[h,ncol(patient_cuSurvival)]==colnames(patient_gene_expression))
      
    }
    
    #Remove the patients who were not in the codon usage and survival results.
    patient_gene_expression <- patient_gene_expression[,filter_file_names]
    
    #Write the expression matrix file. 
    save(patient_gene_expression, file = paste0(i, "_patient_gene.expression"))
    write.csv(patient_gene_expression, file = paste0(i, "_", "gencode_patient_gene_expression.csv"))
    
  } #End gene sets loop.
  
} #End cancer types loop.

#Stop the parallel back-end.
stopCluster(cl)
