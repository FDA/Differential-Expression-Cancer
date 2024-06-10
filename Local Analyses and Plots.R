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
library(tidyverse)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggh4x)

###############################################################  

#Cancer types modeled in this study.
cancer_types <- c("KIRC", "HNSC", "LGG", "LUSC")

#Different types of codon quartile analysis.
types <- c("AllOncogenes", "MutatedOncogenes", "AllTSG", "MutatedTSG", "AllMutations")

############################################################### 

for (c in cancer_types) {
  
  ###############################################################  
  
  if (c == "LUSC") {
    
    #LUSC did not have any genes in its survival-associated gene set.
    gene_sets <- c("original", "exp_threshold", "dge")
    
  } else {
    
    #Gene sets to loop over for the rest of the cancer types
    gene_sets <- c("original", "exp_threshold", "dge", "significant_genes")
    
  } #End gene set if
  
  ###############################################################  
  
  #Strings for paths to reference later in the script.
  #Change "..." to a path relevant to your file system.
  oncogene_analysis_path <- ".../Oncogene Analysis/"
  significant_codon_path <- paste0(".../Survival Codon Analysis/", c, "/")
  tcga_path <- paste0(".../TCGA-Projects/", c, "/")
  cancer_path <- paste0(oncogene_analysis_path, c, "/")
  cancer_input_path <- paste0(cancer_path,"Input Data/")
  cancer_output_path <- paste0(cancer_path,"Output Data/")
  codon_path <- paste0(cancer_path,"Codon Plots/")
  
  ############################################################### 
  
  setwd(oncogene_analysis_path)
  
  #Load ENSG identifiers from Ongene, and remove their version numbers for simpler matching.
  #bioinfo-minzhao.org
  oncogenes <- read.csv(file = "ongene_human.csv")
  oncogene_ensg <- stringr::str_extract(oncogenes[,4], "(?<=ENSG)\\d+")
  ensg_na <- is.na(oncogene_ensg)
  oncogene_ensg <- oncogene_ensg[which(ensg_na==F)]
  oncogene_ensg <- paste0("ENSG", oncogene_ensg)
  
  #Load ENSG identifiers from TSgene, and remove their version numbers for simpler matching.
  #https://bioinfo.uth.edu/TSGene/
  tsg <- read.csv(file = "Human_TSGs.csv")
  tsg_ensg <- stringr::str_extract(tsg[,4], "(?<=ENSG)\\d+")
  ensg_na <- is.na(tsg_ensg)
  tsg_ensg <- tsg_ensg[which(ensg_na==F)]
  tsg_ensg <- paste0("ENSG", tsg_ensg)
  
  ###############################################################  
  
  #Load the cancer type's mutation table.
  setwd(cancer_input_path)
  
  tcga_mutation <- read.csv(file = paste0(c, "_mutation_table.csv"))
  
  ############################################################### 
  
  #Prepare each cancer type's data for analysis and plotting.
  setwd(cancer_input_path)
  
  gene_set_data <- foreach(set = gene_sets) %do% {
    
    #Load each gene set's reference codon usage table and extract the ENSG identifiers.
    load(file = paste0(set, "_gencode_gene_codon.counts"))
    ensg <- unlist(gene_codon_counts[,2])
    
    #Extract the sense codon counts. 
    gene_codon_counts <- gene_codon_counts[,4:ncol(gene_codon_counts)]
    gene_codon_counts_codon_names <- colnames(gene_codon_counts)
    stop_codon_columns <- c(which(gene_codon_counts_codon_names=="TAG"),
                   which(gene_codon_counts_codon_names=="TGA"),
                   which(gene_codon_counts_codon_names=="TAA"))
    gene_codon_counts <- gene_codon_counts[,-stop_codon_columns]
    gene_codon_counts_codon_names <- gene_codon_counts_codon_names[-stop_codon_columns]
    
    #Convert the counts to a matrix.
    gene_codon_counts <- foreach(n = 1:ncol(gene_codon_counts), .combine = 'cbind', .multicombine = T) %do% {
      
      unlist(gene_codon_counts[,n])
      
    } 
    
    colnames(gene_codon_counts) <- gene_codon_counts_codon_names
    rownames(gene_codon_counts) <- ensg
    
    #Translate the codons into their amino acids. 
    gene_counts_aa <- foreach(n = colnames(gene_codon_counts), .combine = 'c') %do% {
      
      seqinr::translate(unlist(strsplit(n, split = character(0))))
      
    }
    
    #Order the codons and amino acids in alphabetical amino acid order.
    aa_order <- order(gene_counts_aa)
    gene_codon_counts <- gene_codon_counts[,aa_order]
    gene_counts_aa <- gene_counts_aa[aa_order]
    
    ############################################################### 

    #Calculate the RSCU of each gene in the reference codon usage table.
    dge_rscu <- foreach(b = 1:ncol(gene_codon_counts), .combine = 'cbind', .multicombine = T) %do% {
        
        #What amino acid is coded by this codon?
        aa <- gene_counts_aa[b]

        #How many synonymous codons encode the same amino acid?
        syn_inds <- which(gene_counts_aa==aa)
        degen <- length(syn_inds)
        
        #How many of this codon were observed in each gene?
        cod_count <- gene_codon_counts[,b]
        
        #If there is only one synonymous codon, syn_count is a vector.
        #If there are more than one it is a matrix.
        if (degen==1) {
          syn_count <- gene_codon_counts[, syn_inds]
        } else {
          syn_count <- rowSums(gene_codon_counts[, syn_inds])
        }
        
        #Calculate the codon's RFac for each gene.
        cod_rscu <- cod_count/syn_count

        return(cod_rscu)

      } #End dge_rscu
    
    colnames(dge_rscu) <- colnames(gene_codon_counts)
    rownames(dge_rscu) <- rownames(gene_codon_counts)
    
    #Save the reference RSCU's and RNA sequence relative frequencies for reference later. 
    save(dge_rscu, file = paste0(set, "_RSCU.mat"))
    load(paste0(set, "_patient_gene.freqs"))
    
    #Calculate the RSCU relative to the entire genome.
    genome_rscu <- foreach(b = 1:ncol(gene_codon_counts), .combine = 'c') %do% {
      
      #What amino acid is coded by this codon?
      aa <- gene_counts_aa[b]
      
      #How many synonymous codons encode the same amino acid?
      syn_inds <- which(gene_counts_aa==aa)
      degen <- length(syn_inds)
      
      #How many of this codon were observed in each gene?
      cod_count <- sum(gene_codon_counts[,b])
      
      syn_count <- sum(gene_codon_counts[, syn_inds])
 
      
      #Calculate the codon's RFac for each gene.
      cod_rscu <- cod_count/syn_count
      
      return(cod_rscu)
      
    } #End genome_rscu
    
    #Save genome_rscu.
    save(genome_rscu, file = paste0(set, "_genome_RSCU.mat"))
    
    ###############################################################  
    
    #This function takes a reference RSCU matrix and each sample's RNA sequence relative frequencies as input,
    metric_results <- function(rscu_mat, expression_frequencies) {
        
      exp_rscu_mats <- foreach(v = 1:ncol(rscu_mat)) %do% {
        
        return(expression_frequencies*rscu_mat[,v])
        
      }
      
      #and outputs expression-weighed RSCU for each gene in the reference RSCU matrix.
      exp_rscu_mats <- foreach(u = 1:length(exp_rscu_mats), .combine = 'cbind', .multicombine = T) %do% {
          
          return(rowSums(exp_rscu_mats[[u]], na.rm = T))
          
      }
        
      gc()
      
      colnames(exp_rscu_mats) <- colnames(rscu_mat)
      rownames(exp_rscu_mats) <- rownames(rscu_mat)
      
      #Output each codon's expression-weighed RSCU, and quantify codon-specific expression-weighed RSCU quartiles.
      all_metric <- foreach(n = 1:ncol(exp_rscu_mats), .combine = 'rbind', .multicombine = T) %do% {
        
        codon_vec <- exp_rscu_mats[,n]
        codon <- rep(colnames(exp_rscu_mats)[n], times = length(codon_vec))
        codon_vec_quartiles <- infotheo::discretize(codon_vec, disc = "equalfreq", nbins = 4)
        
        cbind(codon,
              codon_vec,
              codon_vec_quartiles)
        
      } #End all_metric
      
      colnames(all_metric) <- c("Codon", 
                                "ExprRSCU",
                                "CodonSpecificQuartile")
      return(all_metric)
      
    } #End metric_results.
      
    set_rscu <- metric_results(dge_rscu, patient_gene_frequencies)

    return(set_rscu)
    
  } #End gene_set_data
  
  names(gene_set_data) <- gene_sets
  
  ###############################################################  
  
  #This loop prepares each gene set's data for plotting. 
  setwd(codon_path)
  
  foreach(g = gene_sets) %do% {
    
    ###############################################################   

      #Extract the gene_set dataframe from the list.
      gene_set_ind <- which(names(gene_set_data)==g)
      plot_data <- gene_set_data[[gene_set_ind]]
      original_rownames <- rownames(plot_data)

      #Add codon labels to the data.
      plot_data$Codon <- factor(plot_data$Codon, levels = unique(plot_data[,1]))

      #Add ENSG labels to the data.
      ensg_no_version <- foreach(r = which(plot_data$Codon=="AAA"), .combine = 'c') %do% {

        original_name <- original_rownames[r]
        unlist(strsplit(original_name, split = "\\."))[1]

      }

      ensg_no_version <- rep(ensg_no_version, times = 61)
      plot_data$ENSG <- ensg_no_version

      ###############################################################
      
      #Add labels for oncogenes and their mutations to the data.
      oncogene_overlap <- intersect(plot_data$ENSG, oncogene_ensg)
      oncogene_tcga_overlap <- intersect(unique(tcga_mutation[,"GeneID"]), oncogene_overlap)
      oncogene_overlap <- plot_data$ENSG %in% oncogene_overlap
      plot_data$Oncogene <- oncogene_overlap
      df_oncogene_tcga_overlap <- rep(F, times = nrow(plot_data))

      foreach(o = oncogene_tcga_overlap) %do% {

        mutation_rows <- which(tcga_mutation[,"GeneID"]==o)
        plot_data_rows <- which(plot_data$ENSG==o)

        codons <- toupper(unlist(strsplit(tcga_mutation$CodonChange[mutation_rows], split = "\\/")))
        codon_lost <- codons[c(TRUE, FALSE)]
        codon_gained <- codons[c(FALSE, TRUE)]
        plot_data_codon_rows <- foreach(cod = codon_lost, .combine = 'c') %do% {
          intersect(which(plot_data$Codon==cod), plot_data_rows)
        }

        df_oncogene_tcga_overlap[plot_data_codon_rows] <- T

      }

      plot_data$MutatedOncogene <- df_oncogene_tcga_overlap

      ###############################################################
      
      #Add labels for tumor suppressor genes and their mutations to the data.
      tsg_overlap <- intersect(plot_data$ENSG, tsg_ensg)
      tsg_tcga_overlap <- intersect(unique(tcga_mutation[,"GeneID"]), tsg_overlap)
      tsg_overlap <- plot_data$ENSG %in% tsg_overlap
      plot_data$TSG <- tsg_overlap

      df_tsg_tcga_overlap <- rep(F, times = nrow(plot_data))

      foreach(o = tsg_tcga_overlap) %do% {

        mutation_rows <- which(tcga_mutation[,"GeneID"]==o)
        plot_data_rows <- which(plot_data$ENSG==o)

        codons <- toupper(unlist(strsplit(tcga_mutation$CodonChange[mutation_rows], split = "\\/")))
        codon_lost <- codons[c(TRUE, FALSE)]
        codon_gained <- codons[c(FALSE, TRUE)]
        plot_data_codon_rows <- foreach(cod = codon_lost, .combine = 'c') %do% {
          intersect(which(plot_data$Codon==cod), plot_data_rows)
        }

        df_tsg_tcga_overlap[plot_data_codon_rows] <- T

      }

      plot_data$MutatedTSG <- df_tsg_tcga_overlap

      ###############################################################

      #Add labels for all observed mutations to the data.
      all_gene_mutations <- intersect(plot_data$ENSG, unique(tcga_mutation[,"GeneID"]))

      df_all_gene_mutations <- rep(F, times = nrow(plot_data))

      foreach(o = all_gene_mutations) %do% {

        mutation_rows <- which(tcga_mutation[,"GeneID"]==o)
        plot_data_rows <- which(plot_data$ENSG==o)

        codons <- toupper(unlist(strsplit(tcga_mutation$CodonChange[mutation_rows], split = "\\/")))
        codon_lost <- codons[c(TRUE, FALSE)]
        codon_gained <- codons[c(FALSE, TRUE)]
        plot_data_codon_rows <- foreach(cod = codon_lost, .combine = 'c') %do% {
          intersect(which(plot_data$Codon==cod), plot_data_rows)
        }

        df_all_gene_mutations[plot_data_codon_rows] <- T

      }

      plot_data$AllMutations <- df_all_gene_mutations
      
      ###############################################################
      
      #Summarize the quartile results for oncogenes, tumor suppressor genes, and mutations by codon.
      ind_codon <- foreach(u = 1:4, .combine = 'rbind') %do% {

        ###############################################################
        
        #Extract the quartile's results. 
        quartile_inds <- which(plot_data$CodonSpecificQuartile==u)
        quartile_data <- plot_data[quartile_inds,]

        #Compute each codon's percentage for that quartile. 
        codon_quartile_percentages <- foreach(q = unique(quartile_data$Codon), .combine = 'rbind') %do% {

          codon_quartile_inds <- which(quartile_data$Codon==q)
          codon_quartile_data <- quartile_data[codon_quartile_inds,]

          all_codon_quartile_inds <- which(plot_data$Codon==q)
          all_codon_quartile_data <- plot_data[all_codon_quartile_inds,]

          codon <- as.character(q)
          quartile <- u
          onco_per <- length(which(codon_quartile_data$Oncogene==T))/length(which(all_codon_quartile_data$Oncogene==T))
          mutated_onco_per <- length(which(codon_quartile_data$MutatedOncogene==T))/length(which(all_codon_quartile_data$MutatedOncogene==T))
          tsg_per <- length(which(codon_quartile_data$TSG==T))/length(which(all_codon_quartile_data$TSG==T))
          mutated_tsg_per <- length(which(codon_quartile_data$MutatedTSG==T))/length(which(all_codon_quartile_data$MutatedTSG==T))
          all_mut_per <- length(which(codon_quartile_data$AllMutations==T))/length(which(all_codon_quartile_data$AllMutations==T))

          return(c(codon, quartile, onco_per, mutated_onco_per, tsg_per, mutated_tsg_per, all_mut_per))

        }

        return(codon_quartile_percentages)

      }

      colnames(ind_codon) <- c("Codon", "Quartile", types)

      #Add the results for each individual codon to the data. 
      ind_codon <- as.data.frame(ind_codon)

      ###############################################################

      #Compute each codon's percentage across all quartiles.
      all_type_weighted_quartiles <- foreach(t = types, .combine = "rbind") %do% {

        ind_type_inds <- which(colnames(ind_codon)==t)
        ind_type_data <- as.numeric(ind_codon[,ind_type_inds])
        codons <- ind_codon$Codon
        quartiles <- as.numeric(ind_codon$Quartile)

        weighted_quartiles <- quartiles*ind_type_data

        weighted_quartiles <- tapply(weighted_quartiles, codons, sum)
        return(weighted_quartiles)

      }

      rownames(all_type_weighted_quartiles) <- types
      save(all_type_weighted_quartiles, file = paste(g, c, "cod.quartiles"))

      #Add the results for all_type_weighted_quartiles to the data. 
      gene_set_data[[gene_set_ind]] <- plot_data

  } #end significant_codons_vs_all_codons plot loop

  #Save the results.
  setwd(cancer_output_path)
  save(gene_set_data, file = "gene_set.data")
  
} #end cancer loop

############################################################### 

#Reset gene_sets
gene_sets <- c("original", "exp_threshold", "dge", "significant_genes")

#Load oncogene and tumor suppressor gene data.
oncogene_analysis_path <- "..."
setwd(oncogene_analysis_path)

oncogenes <- read.csv(file = "ongene_human.csv")
tsg <- read.csv(file = "Human_TSGs.csv")

#Extract the ENSG of oncogenes.
oncogene_ensg <- stringr::str_extract(oncogenes[,4], "(?<=ENSG)\\d+")
ensg_na <- is.na(oncogene_ensg)
oncogene_ensg <- oncogene_ensg[which(ensg_na==F)]
oncogene_ensg <- paste0("ENSG", oncogene_ensg)

#Extract the ENSG of tumor suppressor genes.
tsg_ensg <- stringr::str_extract(tsg[,4], "(?<=ENSG)\\d+")
ensg_na <- is.na(tsg_ensg)
tsg_ensg <- tsg_ensg[which(ensg_na==F)]
tsg_ensg <- paste0("ENSG", tsg_ensg)

############################################################### 

#Generate a matrix for plotting significant age and codon usage codons from Kaplan-Meier curves. 
all_significant_codon_data <- foreach(c = cancer_types, .combine = 'rbind', .multicombine = T) %do% {
    
  ############################################################### 
  
    #Load significant age and codon usage data.
    significant_codon_path <- paste0(".../Survival Codon Analysis/", c, "/")
    setwd(significant_codon_path)

    load(file = paste0(c, "_Significant_Codon_Supplemental.File"))
    
    #Ensure that the caner type has significant codon usage results and mark them with "1" for plotting.
    cu_result_summary <- plyr::count(supp[,c(1,2,6)])
    cu_result_summary <- cu_result_summary[which(cu_result_summary[,3]==T),]
    
    if (nrow(cu_result_summary)>0) {
      cu_result_summary[,3] <- 1
    }
    
    colnames(cu_result_summary) <- c("CancerType", "Codons", "CodonSignificance", "Frequency")
    
    #Ensure that the caner type has significant age results and mark them with "1" for plotting.
    age_result_summary <- plyr::count(supp[,c(1,2,10)])
    age_result_summary <- age_result_summary[which(age_result_summary[,3]==T),]
    
    if (nrow(age_result_summary)>0) {
      age_result_summary[,3] <- 2
    }
    
    colnames(age_result_summary) <- c("CancerType", "Codons", "CodonSignificance", "Frequency")
    
    result_summary <- rbind(cu_result_summary, age_result_summary)
    
    #If both coefficients have significant codons, combine their results by codon. 
    if (nrow(result_summary)>0) {
      result_summary <- result_summary %>%
        group_by(CancerType, Codons) %>%
        summarize(across(everything(), sum))
      
      result_summary <- ungroup(result_summary)
      
    }
    
    ############################################################### 
    
    #Generate a matrix containing a cancers results for each gene set.
    significant_codons <- foreach(g = gene_sets, .combine = 'rbind', .multicombine = T) %do% {
      
      ############################################################### 
      
      codons <- toupper(seqinr::words())
      
      #Name the results for each gene set.
      if (g=="original") {
        heatmap_rownames <- rep("AG", times = length(codons)-3)
      } else if (g=="exp_threshold") {
        heatmap_rownames <- rep("ETG", times = length(codons)-3)
      } else if (g=="dge") {
        heatmap_rownames <- rep("DEVG", times = length(codons)-3)
      } else if (g=="significant_genes") {
        heatmap_rownames <- rep("SAG", times = length(codons)-3)
      } else {
        stop("Undefined gene set.")
      }
      
      #Extract the cancer's results. 
      codon_inds <- grep(g, result_summary$CancerType)
      codon_fill <- result_summary$CodonSignificance[codon_inds]
      names(codon_fill) <- result_summary$Codons[codon_inds]

      #Make a vector of zeroes for each codon.
      #Then change significant codons to 1. 
      all_codons_significance <- rep(0, times = length(codons))
      names(all_codons_significance) <- codons
      
      if (length(codon_fill)>0) {
        all_codons_significance <- c(all_codons_significance, codon_fill)
        all_codons_significance <- tapply(all_codons_significance, names(all_codons_significance), sum)
      }

      #Sort codons by amino acid order.
      aa <- foreach(m = names(all_codons_significance), .combine = 'c') %do% {
        
        seqinr::translate(unlist(strsplit(m, split = character(0))))
        
      }
      
      aa_order <- order(aa)
      aa <- aa[aa_order]
      aa <- aa[4:length(aa)]
      all_codons_significance <- all_codons_significance[aa_order]
      all_codons_significance <- all_codons_significance[4:length(all_codons_significance)]
      codons <- names(all_codons_significance)
      
      heatmap_results <- cbind(c, heatmap_rownames, codons, aa, all_codons_significance)
      heatmap_results
      
    } #End significant codons loop.
    
    colnames(significant_codons) <- c("Cancer",
                                      "HeatmapResults",
                                      "Codons",
                                      "AminoAcids",
                                      "Significance")
    significant_codons
    
  }

#Coerce the matrix into a data.frame, and prepare its variables for plotting.
all_significant_codon_data <- as.data.frame(all_significant_codon_data)
all_significant_codon_data$Codons <- factor(all_significant_codon_data$Codons, levels = unique(all_significant_codon_data$Codons))
all_significant_codon_data$AminoAcids <- factor(all_significant_codon_data$AminoAcids, levels = unique(all_significant_codon_data$AminoAcids))
all_significant_codon_data$HeatmapResults <- factor(all_significant_codon_data$HeatmapResults, levels = unique(all_significant_codon_data$HeatmapResults))

#Generate Figure 3 heatmaps.
ggplot(all_significant_codon_data, aes(Codons, HeatmapResults, fill = Significance)) + 
  geom_tile(colour = "black") + 
  facet_nested("Cancer Types" + Cancer ~ "Amino Acids" + AminoAcids, scales = "free", space = "free") +
  theme(strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = paste("Surivival Associated Codons (Bonferroni corrected p<0.01)"),
       x = "Codons",
       y = "Gene Sets") +
  scale_fill_manual(values = c("lightblue","yellow", "red", "orange"), name = "Legend", labels = c("All p>0.01", "Codon usage p<0.01", "Age p<0.01", "All p<0.01")) 

 
setwd("...")

#Save Figure 3 heatmaps.
ggsave(filename = paste("heatmap of significant codons in each cancer and gene set.tiff"), width = 10, height = 7, units = "in")

############################################################### 

#Generate a matrix for plotting differential codon abundance.
all_differential_codon_data <- foreach(c = cancer_types, .combine = 'rbind', .multicombine = T) %do% {
  
  ############################################################### 
  
  #Save paths necessary for matrix.
  oncogene_analysis_path <- "..."
  tcga_path <- "..."
  cancer_path <- paste0(oncogene_analysis_path, c, "/")
  cancer_input_path <- paste0(cancer_path,"Input Data/")
  setwd(cancer_input_path)
  
  #Load differentially expressed and variable gene RSCU values.
  load(file = paste0("dge_RSCU.mat"))
  dge_rscu <- dge_rscu
  dge_rscu_sums <- colSums(dge_rscu, na.rm = T)
  dge_rscu_avg <- dge_rscu_sums/nrow(dge_rscu)
  
  #Load expression threshold gene RSCU values.
  load(file = paste0("exp_threshold_RSCU.mat"))
  exp_threshold_rscu <- dge_rscu
  exp_threshold_rscu_sums <- colSums(dge_rscu, na.rm = T)
  exp_threshold_rscu_avg <- exp_threshold_rscu_sums/nrow(exp_threshold_rscu)
  
  #Load genome RSCU.
  load(file = paste0("dge_genome_RSCU.mat"))

  ############################################################### 
  
  #Load significant genes.
  setwd(tcga_path)
  fn <- paste0("...")
  dat <- read.csv(fn, header=T, row.names = 1) 
  
  #Identify significantly differentially expressed and differentially variable genes.
  dv_inds <- which(dat[,which(colnames(dat)=="dV.q.val")]<=0.05)
  de_inds <- which(dat[,which(colnames(dat)=="dE.q.val")]<=0.05)
  
  #Extract their results to matrices.
  dv <- dat[dv_inds ,which(colnames(dat)=="dV")]
  names(dv) <- rownames(dat)[dv_inds]
  de <- dat[de_inds ,which(colnames(dat)=="dE")]
  names(de) <- rownames(dat)[de_inds]
  
  expr_de <- dat[,which(colnames(dat)=="dE")]
  names(expr_de) <- rownames(dat)

  #If there are differentially variable genes, then compute their differential abundance.
  if (length(dv)>0) {
    diff_dv <- foreach(d = 1:length(dv), .combine = 'rbind', .multicombine = T) %do% {
      
      rscu_row <- grep(names(dv[d]), rownames(dge_rscu))
      rscu <- dge_rscu[rscu_row,]
      dv[d] * rscu
      
    }
    
    if (!is.null(nrow(diff_dv))) {
      diff_dv <- colSums(diff_dv, na.rm = T)
      diff_dv <- diff_dv/dge_rscu_sums

    }
    
  }
  
  #If there are differentially expressed genes, then compute their differential abundance.
  if (length(de)>0) {
    diff_de <- foreach(d = 1:length(de), .combine = 'rbind', .multicombine = T) %do% {
      
      rscu_row <- grep(names(de[d]), rownames(dge_rscu))
      rscu <- dge_rscu[rscu_row,]
      de[d] * rscu
      
    }
    
    if (!is.null(nrow(diff_de))) {
      diff_de <- colSums(diff_de, na.rm = T)/max(abs(colSums(diff_de, na.rm = T)))

    }
    
  }
  
  #If there are differentially expressed or differentially variable genes, then combine their results. 
  if ((length(de)>0 && length(dv)>0)) {

    dat <- c(diff_de, diff_dv)
    dat_labels <- sort(rep(c("DE", "DV"), times = length(dat)/2))
    print(paste(c, "both"))
  } else if (length(de)==0) { #If there are no differentially expressed genes, then only save the differentially variable results.
    dat <- diff_dv
    print(paste(c, "de = 0"))
    dat_labels <- rep(c("DV"), times = length(dat))
  }else if (length(dv)==0) {#If there are no differentially variable genes, then only save the differentially expressed results.

    dat <- diff_de
    
    print(paste(c, "dv = 0"))
    dat_labels <- rep(c("DE"), times = length(dat))
  }
  
  
  aa <- foreach(m = names(dat), .combine = 'c') %do% {
    
    seqinr::translate(unlist(strsplit(m, split = character(0))))
    
  }
  
  
  #Column bind the results into a matrix and add amino acid labels.
  cbind(rep(c, times = length(dat)), dat_labels, aa, names(dat), dat)
  
}

colnames(all_differential_codon_data) <- c("Cancer",
                                           "Differential",
                                           "AminoAcids",
                                           "Codons",
                                           "Value")

#Coerce the matrix into a data.frame, and prepare its variables for plotting.
all_differential_codon_data <- as.data.frame(all_differential_codon_data)
all_differential_codon_data$Codons <- factor(all_differential_codon_data$Codons, levels = unique(all_differential_codon_data$Codons))
all_differential_codon_data$AminoAcids <- factor(all_differential_codon_data$AminoAcids, levels = unique(all_differential_codon_data$AminoAcids))
all_differential_codon_data$Differential <- factor(all_differential_codon_data$Differential, levels = unique(all_differential_codon_data$Differential))
all_differential_codon_data$Value <- as.numeric(all_differential_codon_data$Value)
all_differential_codon_data <- all_differential_codon_data[which(all_differential_codon_data$Differential=="DE"),]

#Plot figure 5.
ggplot(all_differential_codon_data, aes(Codons, Differential, fill = Value)) + 
  geom_tile(colour = "black") + 
  facet_nested("Cancer Types" + Cancer ~ "Amino Acids" + AminoAcids, scales = "free", space = "free", switch = "y") +
  theme(strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(title = paste("Differential Codon Analysis Results"),
       x = "Codons") +
  scale_fill_gradient2(low="#7fbf7b", mid="white", high="#af8dc3", name = "RF weighted Dg")


setwd("...")

#Save figure 5.
ggsave(filename = paste("heatmap of differential codon changes.tiff"), width = 10, height = 7, units = "in")

############################################################### 

#Generate a matrix for plotting mutational deltas.
all_mutation_codon_data <- foreach(c = cancer_types, .combine = 'rbind', .multicombine = T) %do% {

  #Save paths for reference.
  oncogene_analysis_path <- "..."
  cancer_path <- paste0(oncogene_analysis_path, c, "/")
  cancer_input_path <- paste0(cancer_path,"Input Data/")
  setwd(cancer_input_path)
  
  #Read mutation data.
  tcga_mutation <- read.csv(file = paste0(c, "..."))
  
  #Codon changes are listed in a WildTypeCodon//MutatedCodon format in the .maf.
  #This code separates WildTypeCodon and MutatedCodon.
  codon_changes <- foreach(r = 1:nrow(tcga_mutation), .combine = 'rbind') %do% {
    
    codons <- toupper(unlist(strsplit(tcga_mutation$CodonChange[r], split = "\\/")))
    gene_symbol <- tcga_mutation$GeneSymbol[r]
    ensg <- tcga_mutation$GeneID[r]
    codon_lost <- codons[1]
    codon_gained <- codons[2]
    
    c(gene_symbol, ensg, codon_lost, codon_gained)
    
  } #end codon_changes foreach
  
  colnames(codon_changes) <- c("GeneSymbol",
                               "ENSG", 
                               "CodonLost",
                               "CodonGained")
  
  codon_changes <- as.data.frame(codon_changes)
  
  #Compute mutational deltas for each codon and format into a matrix.
  all_codon_delta <- foreach(g = gene_sets, .combine = 'rbind', .multicombine = T) %do% {
    
    #Label each gene set's results.
    if (g=="original") {
      heatmap_rownames <- rep("AG", times = length(codons))
    } else if (g=="exp_threshold") {
      heatmap_rownames <- rep("ETG", times = length(codons))
    } else if (g=="dge") {
      heatmap_rownames <- rep("DEVG", times = length(codons))
    } else if (g=="significant_genes") {
      heatmap_rownames <- rep("SAG", times = length(codons))
    } else {
      stop("Undefined gene set.")
    }
    
    load(paste0(g, "_gene.lengths"))
    
    gene_set_names <- foreach(n = names(gene_lengths), .combine = 'c', .multicombine = T) %do% {
      unlist(strsplit(n, split = "\\."))[1]
    }
    
    #Identify which codons changed occurred in the set.
    change_in_set <- intersect(codon_changes[,2], gene_set_names)
    
    set_codon_changes <- codon_changes[codon_changes$ENSG %in% change_in_set,]
    
    #Compute codon gains, losses, and delta. 
    overall_codon_loss <- plyr::count(set_codon_changes[,3])
    codon_loss <- overall_codon_loss[,2]*-1
    names(codon_loss) <- overall_codon_loss[,1]
    overall_codon_gain <- plyr::count(set_codon_changes[,4])
    codon_gain <- overall_codon_gain[,2]
    names(codon_gain) <- overall_codon_gain[,1]
    
    #Make a vector of zeros and label it by codon.
    empty_codons <- rep(0, times = 64)
    names(empty_codons) <- toupper(seqinr::words())
    
    #Fill the vector with codon deltas.
    codon_delta <- c(codon_loss, codon_gain, empty_codons)
    codon_delta <- tapply(codon_delta, names(codon_delta), sum)
    codon_delta <- codon_delta/length(change_in_set)

    #Label and organize codon deltas by amino acid.
    aa <- foreach(m = names(codon_delta), .combine = 'c') %do% {
      
      seqinr::translate(unlist(strsplit(m, split = character(0))))
      
    }
    
    aa_order <- order(aa)
    codon_delta <- codon_delta[aa_order]
    aa <- aa[aa_order]
    
    heatmap_results <- cbind(rep(c, times = length(codons)), heatmap_rownames, names(codon_delta), aa, codon_delta)
    heatmap_results
    
  }
  
  colnames(all_codon_delta) <- c("Cancer",
                                 "HeatmapResults",
                                 "Codons",
                                 "AminoAcids",
                                 "CodonDelta")
  
  all_codon_delta
   
}

#Coerce the matrix into a data.frame, and prepare its variables for plotting.
all_mutation_codon_data <- as.data.frame(all_mutation_codon_data)
all_mutation_codon_data$Codons <- factor(all_mutation_codon_data$Codons, levels = unique(all_mutation_codon_data$Codons))
all_mutation_codon_data$AminoAcids <- factor(all_mutation_codon_data$AminoAcids, levels = unique(all_mutation_codon_data$AminoAcids))
all_mutation_codon_data$HeatmapResults <- factor(all_mutation_codon_data$HeatmapResults, levels = unique(all_mutation_codon_data$HeatmapResults))
all_mutation_codon_data$CodonDelta <- as.numeric(all_mutation_codon_data$CodonDelta)

#Plot figure six.
ggplot(all_mutation_codon_data, aes(Codons, HeatmapResults, fill = CodonDelta)) + 
  geom_tile(colour = "black") + 
  facet_nested("Cancer Types" + Cancer ~ "Amino Acids" + AminoAcids, scales = "free", space = "free") +
  theme(strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title = paste("Mutational Codon Changes Per Mutated Genes in Set"),
       x = "Codons",
       y = "Gene Sets") +
  scale_fill_gradient2(low="#006400", mid="white", high="#8B0000",
                       midpoint=0, name = "Legend") 


setwd("...")

#Save figure 6.
ggsave(filename = paste("heatmap of mutational codon deltas.tiff"), width = 10, height = 7, units = "in")

############################################################### 

#Generate a matrix for plotting composite plot 7. 
all_weighted_quartiles_data <- foreach(c = cancer_types, .combine = 'rbind', .multicombine = T) %do% {
  
  ############################################################### 
  
  #Only iterate over large enough gene sets for meaningful quartiles.
  if (c == "LUSC") {
    gene_sets <- c("original", "exp_threshold")
  } else   if (c == "HNSC") {
    gene_sets <- c("original", "exp_threshold")
  } else   if (c == "KIRC") {
    gene_sets <- c("original", "exp_threshold", "dge")
  } else {
    
    gene_sets <- c("original", "exp_threshold", "dge", "significant_genes")
    
  }
  
  #Bind gene sets into a matrix.
  #The matrices are then bound into an all-cancer matrix by the first loop.
  foreach(g = gene_sets, .combine = 'rbind', .multicombine = T) %do% {
    
    ############################################################### 
    
    #Save paths for reference.
    oncogene_analysis_path <- "..."
    cancer_path <- paste0(oncogene_analysis_path, c, "/")
    codon_path <- paste0(cancer_path,"Codon Plots/")
    
    setwd(codon_path)
    
    #Load quartile data for the gene set and cancer.
    load(file = paste(g, c, "cod.quartiles"))
    
    #Generate labels for the data.
    cancer <- rep(c, times = ncol(all_type_weighted_quartiles))
    if (g=="original") {
      heatmap_rownames <- rep("AG", times = ncol(all_type_weighted_quartiles))
    } else if (g=="exp_threshold") {
      heatmap_rownames <- rep("ETG", times = ncol(all_type_weighted_quartiles))
    } else if (g=="dge") {
      heatmap_rownames <- rep("DEVG", times = ncol(all_type_weighted_quartiles))
    } else if (g=="significant_genes") {
      heatmap_rownames <- rep("SAG", times = ncol(all_type_weighted_quartiles))
    } else {
      stop("Undefined gene set.")
    }
    
    codons <- colnames(all_type_weighted_quartiles)
    
    #Label and order the data by amino acid.
    aa <- foreach(m = codons, .combine = 'c') %do% {
      
      seqinr::translate(unlist(strsplit(m, split = character(0))))
      
    }
    
    aa_order <- order(aa)
    codons <- codons[aa_order]
    all_type_weighted_quartiles <- all_type_weighted_quartiles[,aa_order]
    aa <- aa[aa_order]
    
    heatmap_results <- cbind(cancer, heatmap_rownames, codons, aa, t(all_type_weighted_quartiles))
    
  }
  
} 

colnames(all_weighted_quartiles_data) <- c("Cancer",
                                           "GeneSet",
                                           "Codons",
                                           "AminoAcids",
                                           "AllOncogenes",
                                           "MutatedOncogenes",
                                           "AllTSG",
                                           "MutatedTSG",
                                           "AllMutations")

#Coerce the matrix into a data.frame, and prepare its variables for plotting.
all_weighted_quartiles_data <- as.data.frame(all_weighted_quartiles_data)
all_weighted_quartiles_data$Codons <- factor(all_weighted_quartiles_data$Codons, levels = unique(all_weighted_quartiles_data$Codons))
all_weighted_quartiles_data$AminoAcids <- factor(all_weighted_quartiles_data$AminoAcids, levels = unique(all_weighted_quartiles_data$AminoAcids))
all_weighted_quartiles_data$GeneSet <- factor(all_weighted_quartiles_data$GeneSet, levels = unique(all_weighted_quartiles_data$GeneSet))

#Plot and save the individual plots that make up composite plot 7.
foreach(t=types) %do% {
  
  all_weighted_quartiles_data[,t] <- as.numeric(all_weighted_quartiles_data[,t])
  
  ggplot(all_weighted_quartiles_data, aes(Codons, GeneSet, fill = all_weighted_quartiles_data[,t])) + 
    geom_tile(colour = "black") + 
    facet_nested("Cancer Types" + Cancer ~ "Amino Acids" + AminoAcids, scales = "free", space = "free") +
    theme(strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(title = paste(t, "Frequency Weighted Quartiles"),
         x = "Codons",
         y = "Gene Sets") +
    scale_fill_gradient2(low="navy", mid="white", high="red",
                         midpoint=2.5, limits=c(1,4), name = "Legend") 
  
  setwd("...")

  ggsave(filename = paste(t, "frequency-weighted codon quartiles.tiff"), width = 10, height = 7, units = "in")
  
} #End foreach(t=types)

