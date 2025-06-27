# how should we define our sorted proteomics signature for Monocyte?
# assuming trends with WV deconvolution will hold true with other deconvolution methods
# 1- what is the best filter for accurately scoring monocytes (based on t-test, correct category per sample)
# 2- what is the best filter for predicting drug sensitivity (based on Ven, Aza+Ven AUC)

library(synapser)
library(DMEA)
library(tidyr)
library(tibble)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
synapser::synLogin()
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/circBar.R")
source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/helperScripts/panSEA_helper_20240913.R")

evalOneMonoSig <- function(temp.sig, BeatAML, type="rna", gmt) {
  # calculate WV scores
  temp.sig <- na.omit(temp.sig)
  global.df100 <- BeatAML$global[,colSums(is.na(BeatAML$global)) == 0]
  if (nrow(temp.sig) > 0) {
    # perform weighted voting on input omics
    temp.sig <- temp.sig[,c("Gene", colnames(temp.sig)[1])]
    temp.sig2 <- temp.sig[temp.sig$Gene %in% colnames(global.df100)[2:ncol(global.df100)],]
    if (nrow(temp.sig2) > 0) {
      global.df100 <- global.df100[,c("Barcode.ID",temp.sig2$Gene)]
      if (nrow(global.df100) > 0) {
        sorted.wv <- panSEA::WV(global.df100, temp.sig2)
        temp.wv <- sorted.wv$scores }
      } else {
        temp.wv <- data.frame()
      }
    } else {
      temp.wv <- data.frame()
    }
    
  AUC.df <- merge(temp.wv, BeatAML$drug, by="Barcode.ID")
  Drug <- colnames(AUC.df)[3:ncol(AUC.df)]
  Barcode.ID <- AUC.df$Barcode.ID
  pt.corr <- data.frame(Barcode.ID, SSE = NA, Pearson.est = NA, Pearson.p = NA,
                        Spearman.est = NA, Spearman.p = NA, N = NA)
  all.corr <- data.frame()
  for (i in 1:nrow(AUC.df)) {
    # leave out one patient at a time
    temp.AUC <- AUC.df[-i,]
    corr <- data.frame(Drug, Slope = NA, Intercept = NA, R.squared = NA, N = NA)
    # for each drug:
    for (j in 3:ncol(temp.AUC)) {
      # determine line of best fit 
      x <- as.numeric(temp.AUC$WV)
      y <- as.numeric(temp.AUC[,j])
      Regression <- stats::lm(y ~ x)
      corr$Slope[j-2] <- Regression$coeff[[2]]
      corr$Intercept[j-2] <- Regression$coeff[[1]]
      corr$R.squared[j-2] <- summary(Regression)$r.squared
      corr$N[j-2] <- length(x)
    }
    corr$Barcode.ID <- AUC.df$Barcode.ID[i]
    
    # predict drug sensitivity using line of best fit
    corr$WV <- temp.wv[temp.wv$Barcode.ID == AUC.df$Barcode.ID[i],]$WV
    corr$AUC_predicted <- corr$Intercept + corr$Slope * corr$WV
    
    # calculate accuracy
    corr$AUC_measured <- as.numeric(as.vector(AUC.df[i,Drug]))
    corr$delta_AUC_squared <- ifelse(is.na(corr$AUC_predicted) | is.na(corr$AUC_measured),
                                     NA, (corr$AUC_predicted - corr$AUC_measured)^2)
    all.corr <- rbind(all.corr, corr)
    pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$SSE <- sum(corr$delta_AUC_squared)
    
    corr.df <- na.omit(corr)
    pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$N <- nrow(corr.df)
    if (nrow(corr.df) > 2) {
      Pearson <- stats::cor.test(corr.df$AUC_measured, corr.df$AUC_predicted, method = "pearson")
      pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$Pearson.est <- Pearson$estimate
      pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$Pearson.p <- Pearson$p.value
      
      Spearman <- stats::cor.test(corr.df$AUC_measured, corr.df$AUC_predicted, method = "spearman")
      pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$Spearman.est <- Spearman$estimate
      pt.corr[pt.corr$Barcode.ID == AUC.df$Barcode.ID[i],]$Spearman.p <- Spearman$p.value 
    }
  }
  
  # for each drug, evaluate accuracy
  drug.corr <- data.frame(Drug, SSE = NA, Pearson.est = NA, Pearson.p = NA,
                        Spearman.est = NA, Spearman.p = NA, N = NA)
  for (j in Drug) {
    drug.df <- na.omit(all.corr[all.corr$Drug == j,])
    drug.corr[drug.corr$Drug == j,]$N <- nrow(drug.df)
    drug.corr$SSE <- sum(drug.df$delta_AUC_squared)
    
    if (nrow(drug.df) > 2) {
      Pearson <- stats::cor.test(drug.df$AUC_measured, drug.df$AUC_predicted, method = "pearson")
      drug.corr[drug.corr$Drug == j,]$Pearson.est <- Pearson$estimate
      drug.corr[drug.corr$Drug == j,]$Pearson.p <- Pearson$p.value
      
      Spearman <- stats::cor.test(drug.df$AUC_measured, drug.df$AUC_predicted, method = "spearman")
      drug.corr[drug.corr$Drug == j,]$Spearman.est <- Spearman$estimate
      drug.corr[drug.corr$Drug == j,]$Spearman.p <- Spearman$p.value 
    }
  }
  
  # run DMEA to see which MOAs have the best accuracy based on corr estimate
  dmeaInput <- na.omit(drug.corr[,c("Drug", "Pearson.est", "Spearman.est")])
  if (nrow(dmeaInput) > 6) {
    dmeaPearson <- panSEA::drugSEA_ties(drug.corr, gmt)
    dmeaSpearman <- panSEA::drugSEA_ties(drug.corr, gmt, rank.metric="Spearman.est") 
  } else {
    dmeaPearson <- list()
    dmeaSpearman <- list()
  }
  
  return(list(wv = temp.wv, drug = all.corr, 
              pt.corr = pt.corr, drug.corr = drug.corr, 
              DMEA = dmeaPearson, DMEA.Spearman = dmeaSpearman))
}

evalMonoSig <- function(sig.matrix, BeatAML, types=rep("rna",ncol(sig.matrix)), gmt) {
  # evaluate each signature
  wv.df <- data.frame()
  drug.df <- data.frame()
  pt.corr.df <- data.frame()
  drug.corr.df <- data.frame()
  sig.matrix$Gene <- rownames(sig.matrix)
  DMEA.results <- list()
  DMEA.results.Spearman <- list()
  for (i in 1:(ncol(sig.matrix)-1)) {
    cat("evaluating",names(sig.matrix)[i],"as",types[i],"\n")
    temp.result <- evalOneMonoSig(sig.matrix[,c(i,ncol(sig.matrix))],
                                  BeatAML, types[i], gmt)
    DMEA.results[[names(sig.matrix)[i]]] <- temp.result$DMEA
    DMEA.results.Spearman[[names(sig.matrix)[i]]] <- temp.result$DMEA.Spearman
    temp.wv.df <- temp.result$wv
    temp.drug.df <- temp.result$drug
    temp.pt.corr.df <- temp.result$pt.corr
    temp.drug.corr.df <- temp.result$drug.corr
    temp.wv.df$Signature <- names(sig.matrix)[i]
    temp.drug.df$Signature <- names(sig.matrix)[i]
    temp.pt.corr.df$Signature <- names(sig.matrix)[i]
    temp.drug.corr.df$Signature <- names(sig.matrix)[i]
    
    wv.df <- rbind(wv.df, temp.wv.df)
    drug.df <- rbind(drug.df, temp.drug.df)
    pt.corr.df <- rbind(pt.corr.df, temp.pt.corr.df)
    drug.corr.df <- rbind(drug.corr.df, temp.drug.corr.df)
  }
  
  return(list(wv = wv.df, drug = drug.df, pt.corr = pt.corr.df, 
              drug.corr = drug.corr.df, DMEA = DMEA.results,
              DMEA.Spearman = DMEA.results.Spearman))
}

compareSigs <- function(sigs, value.var = "Log2FC", BeatAML, 
                        types=rep("rna",length(sigs)), gmt, 
                        fillVals = RColorBrewer::brewer.pal(length(sigs), "Set2")) {
  # combine signatures into matrix
  filtered.sigs.df <- data.table::rbindlist(sigs, use.names = TRUE, idcol = "Signature")
  sig.matrix <- reshape2::dcast(filtered.sigs.df, Gene ~ Signature, mean,
                                value.var = value.var)
  rownames(sig.matrix) <- sig.matrix$Gene
  sig.matrix$Gene <- NULL
  sig.matrix <- sig.matrix[,names(sigs)]
  
  # test signature matrix
  sigResults <- evalMonoSig(sig.matrix, BeatAML,types,gmt)
  wv.df <- sigResults$wv
  drug.df <- sigResults$drug
  pt.corr.df <- sigResults$pt.corr
  drug.corr.df <- sigResults$drug.corr
  
  # compare accuracy for Aza, Ven, Aza+Ven across signatures
  rank.metrics <- c("Pearson.est", "Spearman.est", "SSE")
  for (i in rank.metrics) {
    descr <- stringr::str_split_1(i, "[.]")[1]
    if ("Drug" %in% colnames(drug.corr.df)) {
      doi <- c("Azacytidine", "Venetoclax", "Azacytidine - Venetoclax")
      drug.corr.df$`Drug Treatment` <- NA
      drug.corr.df[drug.corr.df$Drug == "Azacytidine",]$`Drug Treatment` <- "Aza"
      drug.corr.df[drug.corr.df$Drug == "Azacytidine - Venetoclax",]$`Drug Treatment` <- "Aza + Ven"
      drug.corr.df[drug.corr.df$Drug == "Venetoclax",]$`Drug Treatment` <- "Ven"
      doi.names <- c("Aza", "Aza + Ven", "Ven")
      
      # aza, aza + ven, ven correlations
      plot.df <- drug.corr.df[drug.corr.df$Drug %in% doi,]
      plot.df$rank <- plot.df[,i]
      sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Signature))
      ylab <- paste0(descr," Correlation Estimate")
      if (grepl("SSE", i)) {
        ylab <- "Sum Squared Error"
      }
      ggplot(na.omit(plot.df), aes(x=Signature, y=rank, fill = Signature)) + 
        geom_col(alpha=0.5) + theme_minimal(base_size = 12) + ylab(ylab) + 
        facet_wrap(~ `Drug Treatment`) +
        ggplot2::scale_x_discrete(limits = sigOrder) +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
        scale_fill_manual(values=fillVals, 
                          breaks=c("Sorted","Lasry","Triana","van Galen"))+
        ggtitle("Monocytic signatures predict drug sensitivity")
      ggsave(paste0("AzaVen_DIA_WV_signatureFill_",descr,".pdf"), width = 5, height = 5)
      
      for (j in doi.names) {
        plot.df <- drug.corr.df[drug.corr.df$`Drug Treatment` == j,]
        plot.df$rank <- plot.df[,i]
        sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Signature))
        plot.df$alpha <- 0.5
        circBar(na.omit(plot.df), x="Signature", y = "rank", fill = "Signature", 
                alpha = "alpha", ymin = 0, ymax = 1, alpha_range=0.5, 
                ytick_yScale = 2/3, ytick_yShift = 0, fillVals=fillVals,
                title=paste("Monocytic signatures predict", j, "sensitivity"),
                fname=paste0(j,"_", descr, "_DIA_WV_signatureFill_circBarPlot.pdf"))
        ggplot(plot.df, aes(x=Signature, y=rank, fill = Signature)) + 
          geom_col(alpha=0.5) + theme_classic(base_size = 12) + 
          ylab(ylab) + 
          ggplot2::scale_x_discrete(limits = sigOrder) +
          scale_fill_manual(values=fillVals, 
                            breaks=c("Sorted","Lasry","Triana","van Galen"))+
          ggtitle(paste("Monocytic signatures predict", j, "sensitivity"))
        ggsave(paste0(j,"_", descr, "_DIA_WV_signatureFill_barPlot.pdf"), width = 5, height = 5)
      }
    }
  }
  
  # compare accuracy across all drugs for each signature
  rank.metrics <- c("Pearson.est", "Spearman.est", "SSE")
  for (i in rank.metrics) {
    descr <- stringr::str_split_1(i, "[.]")[1]
    if ("Drug" %in% colnames(drug.corr.df)) {
      plot.df <- drug.corr.df
      plot.df$rank <- plot.df[,i]
      sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Signature))
      ylab <- paste0(descr," Correlation Estimate")
      if (grepl("SSE", i)) {
        ylab <- "Sum Squared Error"
      }
      ggplot(na.omit(plot.df), aes(x=Signature, y=rank, fill = Signature)) + 
        geom_col(alpha=0.5) + theme_minimal(base_size = 12) + ylab(ylab) + 
        facet_wrap(~ Drug) +
        ggplot2::scale_x_discrete(limits = sigOrder) +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
        scale_fill_manual(values=fillVals, 
                          breaks=c("Sorted","Lasry","Triana","van Galen"))+
        ggtitle("Monocytic signatures predict drug sensitivity")
      ggsave(paste0("Drug_DIA_WV_signatureFill_",descr,".pdf"), width = 5, height = 5)
    }
  }
  
  return(list(wv = wv.df, drug = drug.df, pt.corr = pt.corr.df,
              drug.corr = drug.corr.df, DMEA = sigResults$DMEA,
              DMEA.Spearman = sigResults$DMEA.Spearman))
}

load_not_norm_BeatAML_for_DMEA3 <- function(BeatAML.path = "BeatAML_DMEA_inputs_not_normalized",
                                            exclude.samples = c()) {
  message("Loading Beat AML data for DMEA")
  BeatAML_synapse_id <- list("drug_response.csv" = "syn51674470", 
                             "Ex10_metadata.txt" = "syn25807733",
                             "ptrc_ex10_crosstab_global_gene_original.txt" = "syn25714254",
                             "ptrc_ex10_crosstab_phospho_siteID_original.txt" = "syn25714936")
  
  ### download files if any not already downloaded
  if (!file.exists(BeatAML.path)) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  } else if (!any(FALSE %in% lapply(names(BeatAML_synapse_id), file.exists))) {
    lapply(BeatAML_synapse_id, synapser::synGet, downloadLocation = BeatAML.path)
  }
  
  ### load files
  drug.BeatAML <- read.csv(file.path(BeatAML.path, names(BeatAML_synapse_id)[1]))
  meta.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[2]), 
                             sep = "\t", header = TRUE)
  global.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[3]),
                               sep = "\t", header = TRUE)
  phospho.BeatAML <- read.table(file.path(BeatAML.path, names(BeatAML_synapse_id)[4]),
                                sep = "\t", header = TRUE)
  rna.BeatAML <- synapser::synTableQuery("select * from syn26545877")$asDataFrame()
  
  ### format BeatAML data for DMEA
  sample.names <- "Barcode.ID"
  
  ## format drug sensitivity data frame
  # format drug.BeatAML wide (samples in first column, drug names for rest of columns)
  drug.BeatAML <- reshape2::dcast(drug.BeatAML, sample_id ~ inhibitor, 
                                  value.var = "auc", fill = NA)
  
  # change sample column name to match expression data
  names(drug.BeatAML)[1] <- sample.names
  
  ## format global proteomics data frame
  # change global.BeatAML column names from SampleID.abbrev to 
  # Barcode.ID to match drug.BeatAML
  global.ids <- names(global.BeatAML)
  
  # remove X and any 0's from start of each column name and then
  # replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
  for(i in seq_len(length(global.ids))){
    global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
    
    if(substring(global.ids[i], 1, 1) == 0){
      global.ids[i] <- substr(global.ids[i], 2, nchar(global.ids[i]))
    }
    
    if(global.ids[i] %in% meta.BeatAML$SampleID.abbrev){
      global.ids[i] <- meta.BeatAML[meta.BeatAML$SampleID.abbrev == global.ids[i], ]$Barcode.ID
    }
  }
  
  # replace global.BeatAML column names 
  names(global.BeatAML) <- global.ids
  
  # subtract sample medians
  sample.names <- colnames(dplyr::select_if(global.BeatAML, is.numeric))
  #global.BeatAML[,sample.names] <- log(global.BeatAML[,sample.names], 2)
  global_sample_coef <- apply(global.BeatAML[,sample.names], 2, median, na.rm = T)
  global.BeatAML[,sample.names] <- sweep(global.BeatAML[,sample.names], 2, global_sample_coef, FUN = '-')
  
  # transpose global.BeatAML so that first column is Barcode.ID and 
  # rest of columns are gene symbols
  global.BeatAML <- as.data.frame(t(global.BeatAML))
  
  # make first column Barcode.ID
  global.BeatAML[,"Barcode.ID"] <- rownames(global.BeatAML)
  global.BeatAML <- 
    global.BeatAML[ , c("Barcode.ID", 
                        names(global.BeatAML[ , 1:(ncol(global.BeatAML)-1)]))]
  
  ## format phospho-proteomics data frame
  # change global.BeatAML column names from SampleID.abbrev to Barcode.ID to match drug.BeatAML
  phospho.ids <- names(phospho.BeatAML)
  
  # remove X and any 0's from start of each column name and then
  # replace SampleID.abbrev with Barcode.ID to match drug.BeatAML
  for(i in seq_len(length(phospho.ids))){
    phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
    
    if(substring(phospho.ids[i], 1, 1) == 0){
      phospho.ids[i] <- substr(phospho.ids[i], 2, nchar(phospho.ids[i]))
    }
    
    if(phospho.ids[i] %in% meta.BeatAML$SampleID.abbrev){
      phospho.ids[i] <- meta.BeatAML[
        meta.BeatAML$SampleID.abbrev == phospho.ids[i], ]$Barcode.ID
    }
  }
  
  # replace phospho.BeatAML column names
  names(phospho.BeatAML) <- phospho.ids
  
  # subtract sample medians
  sample.names <- colnames(dplyr::select_if(phospho.BeatAML, is.numeric))
  #phospho.BeatAML[,sample.names] <- log(phospho.BeatAML[,sample.names], 2)
  phospho_sample_coef <- apply(phospho.BeatAML[,sample.names], 2, median, na.rm = T)
  phospho.BeatAML[,sample.names] <- sweep(phospho.BeatAML[,sample.names], 2, phospho_sample_coef, FUN = '-')
  
  # transpose phospho.BeatAML so that first column is Barcode.ID and rest of columns are gene symbols
  phospho.BeatAML <- as.data.frame(t(phospho.BeatAML))
  
  # make first column Barcode.ID
  phospho.BeatAML[, "Barcode.ID"] <- rownames(phospho.BeatAML)
  phospho.BeatAML <- phospho.BeatAML[ , c("Barcode.ID", names(phospho.BeatAML[ , 1:(ncol(phospho.BeatAML)-1)]))]
  
  ## format rnaSeq data frame
  rna.BeatAML$Barcode.ID <- rna.BeatAML$labId
  rna.BeatAML <- reshape2::dcast(rna.BeatAML, Barcode.ID ~ display_label, mean,
                                 value.var = "RNA counts")
  rownames(rna.BeatAML) <- rna.BeatAML$Barcode.ID
  
  return(list(meta = meta.BeatAML[!(meta.BeatAML$Barcode.ID %in% exclude.samples),], 
              drug = drug.BeatAML[!(drug.BeatAML$Barcode.ID %in% exclude.samples),], 
              rna = rna.BeatAML[!(rna.BeatAML$Barcode.ID %in% exclude.samples),],
              global = global.BeatAML[!(global.BeatAML$Barcode.ID %in% exclude.samples),],
              phospho = phospho.BeatAML[!(phospho.BeatAML$Barcode.ID %in% exclude.samples),]))
}

#### predict using full signatures ####
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
setwd("data")

# drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa_2025-01-20.csv",
#                       stringsAsFactors=FALSE, fileEncoding="latin1")
# gmt.drug <- DMEA::as_gmt(drug.info, sep=", ")
# saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA_2025-01-20.rds")
gmt.drug <- readRDS("gmt_BeatAML_drug_MOA_2025-01-20.rds")

# load sorted proteomics signature
sig.paths <- list("Sorted" = "analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_noNA.csv",
                  "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
                  "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HSC_protein-coding.csv")

# import signatures and filter
sigs <- list()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
for (i in names(sig.paths)) {
  sigs[[i]] <- read.csv(sig.paths[[i]])
  sigs[[i]] <- na.omit(sigs[[i]][sigs[[i]]$adj.P.Val <= 0.05,c("Gene","Log2FC")])
}

setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create("Monocyte_vs_progenitor_signatures_beadOnly_LOO_2025-05-30")
setwd("Monocyte_vs_progenitor_signatures_beadOnly_LOO_2025-05-30")

dia.wo.out <- readRDS("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/DIA_2batches_noOutliers.rds")
# sorted.patients <- c("18-00105", "21-00839", "22-00571", "22-00117", "16-01184",
#                      "19-00074", "18-00103", "21-00432", "17-01060", "22-00251")
sorted.patients <- unique(dia.wo.out$meta$patient)
BeatAML <- load_not_norm_BeatAML_for_DMEA3(exclude.samples=sorted.patients)

# evaluate signature
evalResults <- compareSigs(sigs, BeatAML = BeatAML, types=c("global", "rna", "rna", "rna"), gmt = gmt.drug) 
write.csv(evalResults$wv, "wv.csv", row.names = FALSE)
write.csv(evalResults$drug, "predictions.csv", row.names = FALSE)
write.csv(evalResults$drug.corr, "drugAccuracy.csv", row.names = FALSE)
write.csv(evalResults$pt.corr, "patientAccuracy.csv", row.names = FALSE)
saveRDS(evalResults$DMEA, "DMEA.rds")
saveRDS(evalResults$DMEA.Spearman, "DMEA_Spearman.rds")
all.DMEA.files <- list()
for (i in names(sigs)) {
  DMEA.files <- list("DMEA_results.csv" =
                       evalResults$DMEA[[i]]$result,
                     "DMEA_results_Spearman.csv" =
                       evalResults$DMEA.Spearman[[i]]$result,
                     "DMEA_volcano_plot.pdf" =
                       evalResults$DMEA[[i]]$volcano.plot,
                     "DMEA_volcano_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$volcano.plot,
                     "DMEA_bar_plot.pdf" =
                       evalResults$DMEA[[i]]$bar.plot,
                     "DMEA_bar_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$bar.plot,
                     "DMEA_dot_plot.pdf" =
                       evalResults$DMEA[[i]]$dot.plot,
                     "DMEA_dot_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$dot.plot) 
  all.DMEA.files[[i]] <- DMEA.files
}
save_to_synapse_v2(all.DMEA.files)

# redo plots
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create("Monocyte_vs_progenitor_signatures_beadOnly_LOO_2025-05-30")
setwd("Monocyte_vs_progenitor_signatures_beadOnly_LOO_2025-05-30")
drug.corr.df <- read.csv("drugAccuracy.csv")
pt.corr.df <- read.csv("patientAccuracy.csv")
p.df <- pt.corr.df
sigOrder <- p.df[order(p.df$Pearson.est),]$Signature
p.df$Signature <- factor(p.df$Signature, levels=unique(sigOrder))
ggplot2::ggplot(p.df, aes(x=Signature, y=Pearson.est)) + geom_violin(alpha=0) +
  geom_point(#aes(color=Drug)
    ) + 
  geom_boxplot(width=0.2, alpha = 0) + 
  labs(y="Pearson Correlation Estimate") + theme_classic(base_size = 12) +
  ggtitle(paste("Monocyte signatures predict drug sensitivity"))
ggsave(paste0("patientAccuracy","_bySignature.pdf"), width = 5, height = 5)

p.df <- drug.corr.df
sigOrder <- p.df[order(p.df$Pearson.est),]$Signature
p.df$Signature <- factor(p.df$Signature, levels=unique(sigOrder))
ggplot2::ggplot(p.df, aes(x=Signature, y=Pearson.est)) + geom_violin(alpha=0) +
  geom_point(#aes(color=Drug)
  ) + 
  geom_boxplot(width=0.2, alpha = 0) + 
  labs(y="Pearson Correlation Estimate") + theme_classic(base_size = 12) +
  ggtitle(paste("Monocyte signatures predict drug sensitivity"))
ggsave(paste0("drugAccuracy","_bySignature.pdf"), width = 5, height = 5)


#frac.corr.df <- read.csv("cellFractionCorrelations.csv")

drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
                      stringsAsFactors = FALSE, fileEncoding = "latin1")
drug.info <- drug.info[,c("Drug","moa")]
drug.info[drug.info$Drug == "Ralimetinib (LY2228820)",]$moa <- "p38 MAPK inhibitor"
drug.info[drug.info$Drug == "Nilotinib",]$moa <- "Abl kinase inhibitor"
drug.info[drug.info$Drug == "AT-101",]$moa <- "BCL inhibitor"
drug.info[is.na(drug.info$moa),]$moa <- "Other"
library(patchwork); library(ggplot2)
pearson.plots <- NULL
spearman.plots <- NULL
#MOAsInTop50 <- names(gmt.drug$genesets)
#moaColors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(MOAsInTop50))
pearson.venn <- list()
spearman.venn <- list()
for (i in unique(drug.corr.df$Signature)) {
  p.df <- drug.corr.df[drug.corr.df$Signature == i,]
  p.df$Pearson.q <- qvalue::qvalue(p = p.df$Pearson.p, pi0 = 1)$qvalues
  p.df$Spearman.q <- qvalue::qvalue(p = p.df$Spearman.p, pi0 = 1)$qvalues
  plot.df <- merge(p.df, drug.info, by="Drug", all.x = TRUE)
  plot.df$Mechanism <- "Other"
  #plot.df[plot.df$moa %in% MOAsInTop50,]$Mechanism <- plot.df[plot.df$moa %in% MOAsInTop50,]$moa
  plot.df[grepl("Venetoclax",plot.df$Drug),]$Mechanism <- "BCL inhibitor"
  plot.df$Drug <- sub(" [(].*", "", plot.df$Drug) # shorten drug names for plot
  plot.df[plot.df$Drug == "NF-kB Activation Inhibitor",]$Drug <- "NFkB Inhibitor"
  
  rank.metrics <- c("Pearson.est", "Spearman.est")
  for (j in rank.metrics) {
    descr <- stringr::str_split_1(j, "[.]")[1]
    if ("Drug" %in% colnames(plot.df)) {
      if (j == "Pearson.est") {
        plot.df <- plot.df[plot.df$Pearson.est > 0 & plot.df$Pearson.q <= 0.05,]
        ylab <- paste0(descr," r")
        pearson.venn[[i]] <- unique(plot.df$Drug)
      } else {
        plot.df <- plot.df[plot.df$Spearman.est > 0 & plot.df$Spearman.q <= 0.05,]
        ylab <- paste0(descr," rho")
        spearman.venn[[i]] <- unique(plot.df$Drug)
      }
      plot.df$rank <- plot.df[,j]
      sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Drug))
      plot.annot <- paste0(i, "\n(", nrow(plot.df), " / ", nrow(p.df), " Drugs Positively Correlated)")
      corr.plot <- ggplot(plot.df, aes(x=Drug, y=rank, fill = Mechanism)) + 
        geom_col() + theme_minimal(base_size = 12) + ylab(ylab) + 
        ggplot2::scale_x_discrete(limits = sigOrder) +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
              axis.title.x=element_blank()) +
        #scale_fill_manual(breaks=MOAsInTop50, values = moaColors) +
        ggtitle(plot.annot) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16), legend.position="bottom")
      ggsave(paste0("Drug_DIA_WV_moaFill_",descr,"_", i, ".pdf"), corr.plot, width = 10, height = 5)
      if (is.null(pearson.plots) & j == "Pearson.est") {
        pearson.plots <- (corr.plot + theme(legend.position = "none"))
      } else if (j == "Pearson.est") {
        pearson.plots <- pearson.plots / (corr.plot + theme(legend.position = "none"))
      } else if (is.null(spearman.plots) & j == "Spearman.est") {
        spearman.plots <- (corr.plot + theme(legend.position = "none"))
      } else if (j == "Spearman.est") {
        spearman.plots <- spearman.plots / (corr.plot + theme(legend.position = "none"))
      }
    }
  }
}
#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
#pearson.plots <- (pearson.plots / plot_spacer()) + plot_layout(guides='collect')
#pearson.plots <- pearson.plots + theme(legend.position = "none")
ggplot2::ggsave("Drug_DIA_WV_moaFill_Pearson_allSigs.pdf", pearson.plots, width=12, height=12)
ggplot2::ggsave("Drug_DIA_WV_moaFill_Spearman_allSigs.pdf", spearman.plots, width=12, height=12)
ggvenn::ggvenn(pearson.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Pearson_sigOverlap.pdf", width=5, height=5)
ggvenn::ggvenn(spearman.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Spearman_sigOverlap.pdf", width=5, height=5)

#### narrow down sorted signature ####
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
setwd("data")

# drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa_2025-01-20.csv",
#                       stringsAsFactors=FALSE, fileEncoding="latin1")
# gmt.drug <- DMEA::as_gmt(drug.info, sep=", ")
# saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA_2025-01-20.rds")
gmt.drug <- readRDS("gmt_BeatAML_drug_MOA_2025-01-20.rds")

# load sorted proteomics signature
sig.paths <- list("Sorted" = "analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_noNA.csv",
                  "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
                  "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HSC_protein-coding.csv")

# import signatures and filter
sigs <- list()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
for (i in names(sig.paths)) {
  sigs[[i]] <- read.csv(sig.paths[[i]])
  sigs[[i]] <- na.omit(sigs[[i]][sigs[[i]]$adj.P.Val <= 0.05,c("Gene","Log2FC")])
}

setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create("Monocyte_vs_progenitor_sigOpt_beadOnly")
setwd("Monocyte_vs_progenitor_sigOpt_beadOnly")

dia.wo.out <- readRDS("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/DIA_2batches_noOutliers.rds")
# sorted.patients <- c("18-00105", "21-00839", "22-00571", "22-00117", "16-01184",
#                      "19-00074", "18-00103", "21-00432", "17-01060", "22-00251")
sorted.patients <- unique(dia.wo.out$meta$patient)
BeatAML <- load_not_norm_BeatAML_for_DMEA3(exclude.samples=sorted.patients)

# evaluate signature
evalResults <- compareSigs(sigs, BeatAML = BeatAML, types=c("global", "rna", "rna", "rna"), gmt = gmt.drug) 
write.csv(evalResults$wv, "wv.csv", row.names = FALSE)
write.csv(evalResults$drug, "predictions.csv", row.names = FALSE)
write.csv(evalResults$drug.corr, "drugAccuracy.csv", row.names = FALSE)
write.csv(evalResults$pt.corr, "patientAccuracy.csv", row.names = FALSE)
saveRDS(evalResults$DMEA, "DMEA.rds")
saveRDS(evalResults$DMEA.Spearman, "DMEA_Spearman.rds")
all.DMEA.files <- list()
for (i in names(sigs)) {
  DMEA.files <- list("DMEA_results.csv" =
                       evalResults$DMEA[[i]]$result,
                     "DMEA_results_Spearman.csv" =
                       evalResults$DMEA.Spearman[[i]]$result,
                     "DMEA_volcano_plot.pdf" =
                       evalResults$DMEA[[i]]$volcano.plot,
                     "DMEA_volcano_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$volcano.plot,
                     "DMEA_bar_plot.pdf" =
                       evalResults$DMEA[[i]]$bar.plot,
                     "DMEA_bar_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$bar.plot,
                     "DMEA_dot_plot.pdf" =
                       evalResults$DMEA[[i]]$dot.plot,
                     "DMEA_dot_plot_Spearman.pdf" =
                       evalResults$DMEA.Spearman[[i]]$dot.plot) 
  all.DMEA.files[[i]] <- DMEA.files
}
save_to_synapse_v2(all.DMEA.files)

# redo plots
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create("Monocyte_vs_progenitor_signatures_beadOnly_LOO_2025-05-30")
setwd("Monocyte_vs_progenitor_signatures_beadOnly_LOO_2025-05-30")
drug.corr.df <- read.csv("drugAccuracy.csv")
pt.corr.df <- read.csv("patientAccuracy.csv")
p.df <- pt.corr.df
sigOrder <- p.df[order(p.df$Pearson.est),]$Signature
p.df$Signature <- factor(p.df$Signature, levels=unique(sigOrder))
ggplot2::ggplot(p.df, aes(x=Signature, y=Pearson.est)) + geom_violin(alpha=0) +
  geom_point(#aes(color=Drug)
  ) + 
  geom_boxplot(width=0.2, alpha = 0) + 
  labs(y="Pearson Correlation Estimate") + theme_classic(base_size = 12) +
  ggtitle(paste("Monocyte signatures predict drug sensitivity"))
ggsave(paste0("patientAccuracy","_bySignature.pdf"), width = 5, height = 5)

p.df <- drug.corr.df
sigOrder <- p.df[order(p.df$Pearson.est),]$Signature
p.df$Signature <- factor(p.df$Signature, levels=unique(sigOrder))
ggplot2::ggplot(p.df, aes(x=Signature, y=Pearson.est)) + geom_violin(alpha=0) +
  geom_point(#aes(color=Drug)
  ) + 
  geom_boxplot(width=0.2, alpha = 0) + 
  labs(y="Pearson Correlation Estimate") + theme_classic(base_size = 12) +
  ggtitle(paste("Monocyte signatures predict drug sensitivity"))
ggsave(paste0("drugAccuracy","_bySignature.pdf"), width = 5, height = 5)


#frac.corr.df <- read.csv("cellFractionCorrelations.csv")

drug.info <- read.csv("~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
                      stringsAsFactors = FALSE, fileEncoding = "latin1")
drug.info <- drug.info[,c("Drug","moa")]
drug.info[drug.info$Drug == "Ralimetinib (LY2228820)",]$moa <- "p38 MAPK inhibitor"
drug.info[drug.info$Drug == "Nilotinib",]$moa <- "Abl kinase inhibitor"
drug.info[drug.info$Drug == "AT-101",]$moa <- "BCL inhibitor"
drug.info[is.na(drug.info$moa),]$moa <- "Other"
library(patchwork); library(ggplot2)
pearson.plots <- NULL
spearman.plots <- NULL
#MOAsInTop50 <- names(gmt.drug$genesets)
#moaColors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(MOAsInTop50))
pearson.venn <- list()
spearman.venn <- list()
for (i in unique(drug.corr.df$Signature)) {
  p.df <- drug.corr.df[drug.corr.df$Signature == i,]
  p.df$Pearson.q <- qvalue::qvalue(p = p.df$Pearson.p, pi0 = 1)$qvalues
  p.df$Spearman.q <- qvalue::qvalue(p = p.df$Spearman.p, pi0 = 1)$qvalues
  plot.df <- merge(p.df, drug.info, by="Drug", all.x = TRUE)
  plot.df$Mechanism <- "Other"
  #plot.df[plot.df$moa %in% MOAsInTop50,]$Mechanism <- plot.df[plot.df$moa %in% MOAsInTop50,]$moa
  plot.df[grepl("Venetoclax",plot.df$Drug),]$Mechanism <- "BCL inhibitor"
  plot.df$Drug <- sub(" [(].*", "", plot.df$Drug) # shorten drug names for plot
  plot.df[plot.df$Drug == "NF-kB Activation Inhibitor",]$Drug <- "NFkB Inhibitor"
  
  rank.metrics <- c("Pearson.est", "Spearman.est")
  for (j in rank.metrics) {
    descr <- stringr::str_split_1(j, "[.]")[1]
    if ("Drug" %in% colnames(plot.df)) {
      if (j == "Pearson.est") {
        plot.df <- plot.df[plot.df$Pearson.est > 0 & plot.df$Pearson.q <= 0.05,]
        ylab <- paste0(descr," r")
        pearson.venn[[i]] <- unique(plot.df$Drug)
      } else {
        plot.df <- plot.df[plot.df$Spearman.est > 0 & plot.df$Spearman.q <= 0.05,]
        ylab <- paste0(descr," rho")
        spearman.venn[[i]] <- unique(plot.df$Drug)
      }
      plot.df$rank <- plot.df[,j]
      sigOrder <- na.omit(unique(plot.df[order(plot.df$rank, decreasing=TRUE),]$Drug))
      plot.annot <- paste0(i, "\n(", nrow(plot.df), " / ", nrow(p.df), " Drugs Positively Correlated)")
      corr.plot <- ggplot(plot.df, aes(x=Drug, y=rank, fill = Mechanism)) + 
        geom_col() + theme_minimal(base_size = 12) + ylab(ylab) + 
        ggplot2::scale_x_discrete(limits = sigOrder) +
        theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1),
              axis.title.x=element_blank()) +
        #scale_fill_manual(breaks=MOAsInTop50, values = moaColors) +
        ggtitle(plot.annot) + 
        theme(plot.title = element_text(hjust = 0.5, face="bold", size=16), legend.position="bottom")
      ggsave(paste0("Drug_DIA_WV_moaFill_",descr,"_", i, ".pdf"), corr.plot, width = 10, height = 5)
      if (is.null(pearson.plots) & j == "Pearson.est") {
        pearson.plots <- (corr.plot + theme(legend.position = "none"))
      } else if (j == "Pearson.est") {
        pearson.plots <- pearson.plots / (corr.plot + theme(legend.position = "none"))
      } else if (is.null(spearman.plots) & j == "Spearman.est") {
        spearman.plots <- (corr.plot + theme(legend.position = "none"))
      } else if (j == "Spearman.est") {
        spearman.plots <- spearman.plots / (corr.plot + theme(legend.position = "none"))
      }
    }
  }
}
#source("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/MPNST/Chr8/MPNST_Chr8_manuscript/Figure_3_Kinase/guides_build_mod.R")
#pearson.plots <- (pearson.plots / plot_spacer()) + plot_layout(guides='collect')
#pearson.plots <- pearson.plots + theme(legend.position = "none")
ggplot2::ggsave("Drug_DIA_WV_moaFill_Pearson_allSigs.pdf", pearson.plots, width=12, height=12)
ggplot2::ggsave("Drug_DIA_WV_moaFill_Spearman_allSigs.pdf", spearman.plots, width=12, height=12)
ggvenn::ggvenn(pearson.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Pearson_sigOverlap.pdf", width=5, height=5)
ggvenn::ggvenn(spearman.venn, show_percentage=FALSE, set_name_size=5, text_size=5)
ggsave("Drug_DIA_WV_Spearman_sigOverlap.pdf", width=5, height=5)
