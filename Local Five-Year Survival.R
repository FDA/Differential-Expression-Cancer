###############################################################  
#
# Code written by: 
#   Haim Bar and Nathan J. Clement
#
# With the collaboration of:
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
library(dplyr)
library(SEMMS)
library(betaMix)
library(survival)
library(igraph)
library(purrr)
library(ggplot2)
library(foreach)
library(doParallel)

#Cancer types modeled in this study.
cancer_types <- c("KIRC", "HNSC", "LUSC", "LGG")

#Gene sets to compute transcriptomic-weighted codon usage. 
canon_names <- c("original", "dge", "exp_threshold", "significant_genes")

###############################################################  

# The following function is part of the betaMix package 
# https://github.com/haimbar/betaMix
plotCluster <- function(AdjMat, clustNo, clusterInfo=NULL, labels=FALSE, nodecol="blue", labelsize=1, figtitle=NULL, edgecols="grey88") {
  if(is.null(clusterInfo))
    clusterInfo <- graphComponents(AdjMat, minCtr=2, type=0)
  if(length(nodecol) < nrow(AdjMat))
    nodecol <- rep(nodecol[1],length=nrow(AdjMat))
  ids <- which(clusterInfo$clustNo == clustNo)
  if (length(ids) > 0) {
    tmpA <- AdjMat[ids,ids]
    tmpclusterInfo <- clusterInfo[ids,]
    rads <- round(10*tmpclusterInfo$distCenter/
                    max(0.1,max(tmpclusterInfo$distCenter)))
    thetas <- rep(0,length(rads))
    intvls <- findInterval(rads,seq(1,10))
    for (intvl in unique(sort(intvls))) {
      pts <- which(intvls == intvl)
      thetas[pts] <- 3*intvl*pi/max(0.01,max(intvls))+seq(0,1.9*pi,length=length(pts))
    }
    sizes <- pmax(0.3,tmpclusterInfo$degree/max(tmpclusterInfo$degree))
    opacity <- 0.25+tmpclusterInfo$intEdges/tmpclusterInfo$degree
    opacity <- 0.7*opacity/max(opacity)
    nodecol <- rgb(t(col2rgb(nodecol)/255),alpha=opacity)[ids]
    plot(rads*cos(thetas), rads*sin(thetas),cex=sizes*1.2, pch=19,axes=FALSE,
         xlab="",ylab="",col=nodecol, main=figtitle,
         ylim=c(min(rads*sin(thetas)), 1.1*max(rads*sin(thetas))))
    for (i in 1:ncol(tmpA)) {
      nbrs <- setdiff(which(abs(tmpA[i,]) == 1), 1:i)
      if(length(nbrs) > 0) {
        edgecol <- rep(edgecols[1], ncol(tmpA))
        if (edgecols[2] %in% colours()) {
          edgecol[which(tmpA[i,nbrs] == -1)] <- edgecols[2]
        }
        for (j in nbrs) {
          lines(c(rads[i]*cos(thetas[i]), rads[j]*cos(thetas[j])),
                c(rads[i]*sin(thetas[i]), rads[j]*sin(thetas[j])),
                col=edgecol[j], lwd=0.5)
        }
      }
    }
    points(rads*cos(thetas), rads*sin(thetas),cex=sizes*1.2, pch=19, col=nodecol)
    if (labels)
      text(rads*cos(thetas), rads*sin(thetas), tmpclusterInfo$labels, pos=3, cex=labelsize)
    ctr <- which(tmpclusterInfo$iscenter==1)
    points(rads[ctr]*cos(thetas[ctr]), rads[ctr]*sin(thetas[ctr]),pch=21,
           cex=sizes[ctr]*2, col="black",lwd=2)
  } else {
    cat("Invalid cluster number\n")
  }
}

###############################################################  

# Loop over four cancer types: "KIRC", "HNSC", "LUSC", "LGG"
for (c in cancer_types) {

  ###############################################################  
  
  #Strings for paths to reference later in the script.
  dir_path <- paste0("C:/Users/Nathan.Clement/OneDrive - FDA/Documents/letter to the editor/Genome Medicine/Survival Codon Analysis/", c)
  figure_path <- paste0(dir_path, "/Figures/")
  setwd(dir_path)
  
  ###############################################################  
  
  # Load the data for each gene set into a matrix.
  dat <- foreach(g = canon_names, .combine = 'rbind') %do% {
  
    fn <- paste(g, "gencode Patient Codon Usage and Survival Data.csv")
    dat <- read.csv(fn, header=T) 
    dat <- dat[,2:ncol(dat)]
    
    dat[,1] <- paste0(g, dat[,1])
    return(dat)
    
  }
  
  ###############################################################  
  
  #dathdr <- read.csv(fn, skip=1, nrows = 1, header = FALSE)
  #colnames(dat) <- dathdr
  colnames(dat)[1] <- "Cancer"
  colnames(dat)[3] <- "TissueType"
  
  #remove patients with unknown age and make the age column integers
  dat <- dat[which(dat$age_at_diagnosis != "'--"),]
  dat$age_at_diagnosis <- as.numeric(dat$age_at_diagnosis)
  dat$End <- as.numeric(dat$end)
  dat <- dat[which(dat$vital_status %in% c("Alive", "Dead")),]
  
  # convert age from days to years, and set the 5-year surival status
  dat$age_at_diagnosis <- as.numeric(dat$age_at_diagnosis)/365
  dat$End <- as.numeric(dat$end)/365
  dat$vital_status <- as.factor(dat$vital_status)
  dat$FYalive <- dat$FYcensor <- dat$FYdead <- rep(FALSE, nrow(dat))
  dat$FYalive[which(dat$End >= 5)] <- TRUE
  dat$FYcensor[which((dat$End < 5) & (dat$vital_status == "Alive"))] <- TRUE
  dat$FYdead <- (!dat$FYalive) & (!dat$FYcensor)
  
  #hist(dat$age_at_diagnosis)
  #table(dat$Cancer)
  #table(dat$vital_status)
  
  # Significance level
  pthres <- 0.05
  codons <- colnames(dat[,5:68])
  stopcodons <- c('TAG','TGA','TAA')
  usecodons <- setdiff(codons, stopcodons)
  minP <- 80 # the minimum required sample size
  
  ###############################################################  
  
  for (cancertype in unique(dat$Cancer)) {
    
    ###############################################################  
    # create a subset by cancer type (separate analysis for each cancer)
    
    cnc <- dat[which((dat$Cancer == cancertype) &
                       (dat$TissueType == "Primary Tumor")), ]
    #cnc <- dat[which((dat$Cancer == cancertype)), ]
    cat(cancertype, "\n")
    # One-codon-at-a-time survival analysis, contrasting the
    # high vs. low codon usage groups
    for (colnm in usecodons) {
      
      ###############################################################  
      
      coduse <- unlist(cnc[colnm])
      grpLo <- which(coduse < quantile(coduse, 0.25))
      grpHi <- which(coduse > quantile(coduse, 0.75))
      tmpdat <- cnc[c(grpLo, grpHi), c('End', 'vital_status', 'age_at_diagnosis')]
      tmpdat$group <- as.factor(c(rep("L", length(grpLo)), 
                                  rep("H", length(grpHi))))
      tmpdat$vital_status <- as.numeric(tmpdat$vital_status) - 1
      #survfit0 <- survdiff(Surv(End, vital_status) ~ group, tmpdat)
      #survfit1 <- logrank_test(Surv(End, vital_status) ~ group, data = tmpdat)
      #pval <- pvalue(survfit1) 
      wb_fit <- survreg(Surv(End, vital_status) ~ group + age_at_diagnosis, 
                        data = tmpdat, dist = "logistic")
      agerange <- seq(min(cnc$age_at_diagnosis), max(cnc$age_at_diagnosis))
      predL <- predict(wb_fit, type = "quantile", p = 1 - 0.5, 
                       newdata = data.frame(group= "L",
                                            age_at_diagnosis=agerange))
      predH <- predict(wb_fit, type = "quantile", p = 1 - 0.5, 
                       newdata = data.frame(group= "H",
                                            age_at_diagnosis=agerange))
      #cat(agerange[which(predL <= 5)][1], "\n")
      #cat(agerange[which(predH <= 5)][1], "\n")
      pval <- summary(wb_fit)$table[2,4]
      # 1-pchisq(survfit0$chisq, 1)
      if (pval < pthres/length(usecodons)) {
        # plot(survfit(Surv(End, vital_status) ~ group, tmpdat),
        #      col=tmpdat$group, main=cancertype)
        cat("\t", colnm, pval*length(usecodons),"\n")
        
        
      } #End if (pval < pthres/length(usecodons))
    } # End for (colnm in usecodons)
    
    ###############################################################  
    
    # All codons at once - find their network structure, and fit
    # a GLM model with SEMMS
    # https://github.com/haimbar/SEMMS
    
    alive <- which(cnc$vital_status == "Alive")
    M <- as.matrix(cnc[alive, c(69, which(colnames(cnc)%in%usecodons))])
    if(nrow(M) < minP) {
      cat("Sample size of Alive is too small\n")
      # next
    }
    MA <- M
    
    dead <- which(cnc$vital_status == "Dead")
    M <- as.matrix(cnc[dead, c(69, which(colnames(cnc)%in%usecodons))])
    if(nrow(M) < minP) {
      cat("Sample size of Dead is too small\n")
      #next
    }
    MB <- M
    # Save the data in a format that SEMMS expects:
    dat4semms <- data.frame(SURV=c(rep(0, nrow(MA)), rep(1, nrow(MB))),
                            rbind(MA, MB))
    
    save(dat4semms, file="dat4semms.RData")
    dataYXZ <- readInputFile("dat4semms.RData", ycol=1, Xcols = 2, 
                             Zcols=3:63, addIntercept = TRUE)
    # The initialization of SEMMS can be done via the edgefinder
    # package, or with a simple threshold method. The former is
    # preferred, but with a couple of the LUSC subset, all the codons
    # are correlated, which causes edgefinder to crash, so in this case
    # we choose the latter initialization method.
    edgefinderinit <- ifelse(cancertype == "LUSC", FALSE, TRUE)
    fittedSEMMS <- fitSEMMS(dataYXZ, mincor=0.3, nn=6, minchange= 1,
                            distribution="P", verbose=T, rnd=F, 
                            initWithEdgeFinder = edgefinderinit)
    fittedGLM <- runLinearModel(dataYXZ, fittedSEMMS$gam.out$nn, "P")
    cat(fittedSEMMS$gam.out$nn, dataYXZ$originalZnames[fittedSEMMS$gam.out$nn],"\n")
    cat(which(fittedSEMMS$gam.out$lockedOut != 0),"\n")

    # Network analysis using the betaMix package
    bmfit <- betaMix(dat4semms)
    # connect pairs of variable (codons, survival, and age) with an edge,
    # if the correlation between them is sufficiently large
    # (according to betaMix)
    A <- getAdjMat(bmfit, signed = TRUE)
    # Plot the network:
    Cairo::Cairo(
      15.24, #length
      15.24, #width
      file = paste0(figure_path, cancertype, "1.tiff"),
      type = "tiff", #tiff
      bg = "transparent", #white or transparent depending on your requirement 
      dpi = 300,
      units = "cm" #you can change to pixels etc 
    )
    
    plot(graph.adjacency(A!=0, mode="undirected"), 
         vertex.label=rownames(A),
         vertex.label.cex=0.5, vertex.size=0, 
         vertex.frame.color="white",
         vertex.label.color='blue',
         edge.color="grey80",  
         asp=1)
    title(cancertype, cex.main=0.7)
    g <- graph.adjacency(A!=0, mode="undirected")
    edg <- get.edgelist(g)
    edgecol <- rep("grey80", nrow(edg))
    for (i in 1:nrow(edg)) {
      if (A[which(rownames(A) == edg[i,1]),
            which(rownames(A) == edg[i,2])] < 0)
        edgecol[i] <- "orange"
    }
    dev.off()
    3
    Cairo::Cairo(
      15.24, #length
      15.24, #width
      file = paste0(figure_path, cancertype, "2.tiff"),
      type = "tiff", #tiff
      bg = "transparent", #white or transparent depending on your requirement 
      dpi = 300,
      units = "cm" #you can change to pixels etc 
    )
    
    plot(g, vertex.label=rownames(A),
         vertex.label.cex=0.5, vertex.size=0, 
         vertex.frame.color="white",
         vertex.label.color='blue', 
         edge.color=edgecol,  asp=1)
    title(cancertype, cex.main=0.7)
    cat(colnames(A)[which(A[1,] != 0)],"\n\n\n")
    dev.off()
  }
  
  diffpred <- est <- pval <- age_est <- age_pval <- rep(0, 18)
  i <- 1
  
  for (cancertype in unique(dat$Cancer)) {
    cnc <- dat[which((dat$Cancer == cancertype) &
                       (dat$TissueType == "Primary Tumor")), ]
    #cnc <- dat[which((dat$Cancer == cancertype)), ]
    cat(cancertype, "\n")
    # One-codon-at-a-time survival analysis, contrasting the
    # high vs. low codon usage groups
    #usecodons
    for (colnm in usecodons) {
      coduse <- unlist(cnc[colnm])
      grpLo <- which(coduse < quantile(coduse, 0.25))
      grpHi <- which(coduse > quantile(coduse, 0.75))
      tmpdat <- cnc[c(grpLo, grpHi), c('End', 'vital_status', 'age_at_diagnosis')]
      tmpdat$group <- as.factor(c(rep("L", length(grpLo)), 
                                  rep("H", length(grpHi))))
      tmpdat$vital_status <- as.numeric(tmpdat$vital_status) - 1
      wb_fit <- survreg(Surv(End, vital_status) ~ group + age_at_diagnosis, 
                        data = tmpdat, dist = "logistic")
      predL <- predict(wb_fit, type = "quantile", p = 1 - 0.5, 
                       newdata = data.frame(group= "L",
                                            age_at_diagnosis=55))
      predH <- predict(wb_fit, type = "quantile", p = 1 - 0.5, 
                       newdata = data.frame(group= "H",
                                            age_at_diagnosis=55))
      est[i] <- summary(wb_fit)$table[2,1]
      pval[i] <- summary(wb_fit)$table[2,4]#*61
      diffpred[i] <- predH - predL
      age_est[i] <- summary(wb_fit)$table[3,1]
      age_pval[i] <- summary(wb_fit)$table[3,4]#*61
      #print(summary(wb_fit)$table)
      i <- i+1
      
      new_dat <- expand.grid(
        group = levels(tmpdat$group), 
        survival = seq(.01, .99, by = .01)
      ) %>%
        mutate(
          pred = map2(group, survival, 
                      ~predict(wb_fit, type = "quantile", p = 1 - .y, se = TRUE, 
                               newdata = data.frame(group = .x, age_at_diagnosis=55))),
          t_plotvar = map_dbl(pred, ~pluck(.x, "fit")),
          se_plotvar = map_dbl(pred, ~pluck(.x, "se.fit")),
          ucl = t_plotvar + 1.96 * se_plotvar,
          lcl = t_plotvar - 1.96 * se_plotvar
        )
      
      palette_group <- c("#E7B800", "#2E9FDF")
      names(palette_group) <- c("H", "L")
      
      new_dat %>%
        ggplot(aes(y = survival)) +
        geom_line(aes(x = t_plotvar, color = group)) +
        geom_ribbon(aes(xmin = lcl, xmax = ucl, fill = group), alpha = 0.2) +
        scale_color_manual(values = palette_group) +
        scale_fill_manual(values = palette_group) +
        xlab("Time (years)") +
        ylab("Survival Probability (%)") + 
        theme_light()
      
      #Uncomment this to make kaplan meier curves for each cancer type, gene set, and sense codon.
      #ggsave(filename = paste0(figure_path, cancertype, "_", colnm, ".tiff"), width = 6, height = 6)
      
    }
  }
  
  #Create cancer name labels for plots.
  cancer_names <- foreach(u = unique(dat$Cancer), .combine = 'c') %do% {
    
    rep(u, times = 61)
    
  }
  
  cnames <- rep(usecodons, times = length(unique(dat$Cancer)))

  #Compute corrected p-values and identify which are significant.
  bon_pval <- pval*61
  bon_age_pval <- age_pval*61
  est_sig <- bon_pval<0.01
  age_sig <- bon_age_pval<0.01
  
  #Combine information into a data table. 
  supp <- cbind(cancer_names, cnames, est, pval, bon_pval, est_sig, age_est, age_pval, bon_age_pval, age_sig)
  colnames(supp) <- c("CancerType",
                      "Codons",
                      
                      "CodonCoefficient",
                      "CodonPValue",
                      "CodonBonferroniPValue",
                      "CodonSignificance",
                      
                      "AgeCoefficient",
                      "AgePValue",
                      "AgeBonferroniPValue",
                      "AgeSignificance")
  
  #Write the results.
  save(supp, file = paste0(c, "_Significant_Codon_Supplemental.File"))
  write.csv(supp, file = paste0(c, "_Significant_Codon_Supplemental_File.csv"))
  
  result_summary <- plyr::count(cbind(cancer_names, cnames, est_sig))
  result_summary <- result_summary[which(result_summary[,3]==T),]
  colnames(result_summary) <- c("CancerType", "Codons", "CodonSignificance", "Frequency")
  write.csv(result_summary, file = paste0(c, "_Significant_Codon_Result_Summary.csv"))
    
  result_summary <- read.csv(file = paste0(c, "_Significant_Codon_Result_Summary.csv"))
  test <- plyr::count(result_summary[,3])
  test[which(test[,2]>1),]
  plyr::count(result_summary[,2])

} #End for (c in cancer_types).

