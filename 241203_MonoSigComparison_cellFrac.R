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
synapser::synLogin()

evalOneMonoSig <- function(global.df100, frac.df, temp.sig, BeatAML, type="rna", gmt) {
  temp.sig <- na.omit(temp.sig)
  if (nrow(temp.sig) > 0) {
    ## sorted proteomics
    # perform weighted voting on sorted proteomics
    temp.sig <- temp.sig[,c("Gene", colnames(temp.sig)[1])]
    temp.sig2 <- temp.sig[temp.sig$Gene %in% colnames(global.df100)[2:ncol(global.df100)],]
    if (nrow(temp.sig2) > 0) {
      global.df100 <- global.df100[,c("Sample",temp.sig2$Gene)]
      if (nrow(global.df100) > 0) {
        sorted.wv <- panSEA::WV(global.df100, temp.sig2)
        temp.wv <- sorted.wv$scores 
        
        # evaluate accuracy based on t-test
        temp.test <- stats::t.test(temp.wv[grepl("CD14", temp.wv$Sample),]$WV,
                                   temp.wv[!grepl("CD14", temp.wv$Sample),]$WV, 
                                   "greater")
        temp.test.df <- as.data.frame(unlist(temp.test))
        colnames(temp.test.df)[1] <- "value"
        temp.test.df$variable <- rownames(temp.test.df)
        
        # compare to known cell fractions
        temp.wv2 <- temp.wv[grepl("CD14",temp.wv$Sample),]
        temp.wv2$labId <- sub("_.*","",temp.wv2$Sample)
        temp.wv2 <- temp.wv2[,c("labId","WV")]
        temp.wv.frac <- merge(frac.df, temp.wv2, by="labId")
        if (nrow(temp.wv.frac) > 2) {
          frac.corr <- cor.test(temp.wv.frac$WV, temp.wv.frac[,2], 
                                method = "pearson")
          N <- nrow(temp.wv.frac)
          Pearson.est <- frac.corr$estimate
          Pearson.p <-frac.corr$p.value
          frac.corr <- cor.test(temp.wv.frac$WV, temp.wv.frac[,2], 
                                method = "spearman")
          Spearman.est <- frac.corr$estimate
          Spearman.p <-frac.corr$p.value
          frac.corr.result <- data.frame(Pearson.est, Pearson.p, 
                                         Spearman.est, Spearman.p, N)
        } else {
          frac.corr.result <- data.frame()
        }
      } else {
        temp.wv <- data.frame()
        temp.test.df <- data.frame()
        frac.corr.result <- data.frame()
      }
    } else {
      temp.wv <- data.frame()
      temp.test.df <- data.frame()
      frac.corr.result <- data.frame()
    }
    
    ## Beat AML
    # perform DMEA on Beat AML global proteomics
    expr <- BeatAML[[type]]
    temp.sig2 <- temp.sig[temp.sig$Gene %in% colnames(expr)[2:ncol(expr)],]
    nSamples <- length(expr$Barcode.ID[expr$Barcode.ID %in% BeatAML$drug$Barcode.ID])
    if (nrow(temp.sig2) > 2 & nSamples > 2) {
      expr <- expr[,c("Barcode.ID", colnames(expr)[colnames(expr) %in% temp.sig2$Gene])]
      if (ncol(expr) > 3 & nrow(expr) > 2) {
        DMEA.result <- panSEA::mDMEA(BeatAML$drug, gmt, list(expr), 
                                     list(temp.sig2), types=type,
                                     sample.names="Barcode.ID",
                                     weight.values=colnames(temp.sig2)[2], 
                                     scatter.plots = FALSE)
        corr.df <- DMEA.result$all.results[[1]]$corr.result
      } else {
        corr.df <- data.frame()
      }
    } else {
      corr.df <- data.frame()
    }
  } else {
    temp.wv <- data.frame()
    temp.test.df <- data.frame()
    corr.df <- data.frame()
    frac.corr.result <- data.frame()
  }
  
  return(list(wv = temp.wv, test = temp.test.df, drug.corr = corr.df, 
              frac.corr = frac.corr.result))
}

evalMonoSig <- function(global.df100, frac.df, sig.matrix, BeatAML, 
                        types=rep("rna",ncol(sig.matrix)), gmt) {
  # evaluate each signature
  wv.df <- data.frame()
  test.df <- data.frame()
  drug.corr.df <- data.frame()
  frac.corr.df <- data.frame()
  sig.matrix$Gene <- rownames(sig.matrix)
  for (i in 1:(ncol(sig.matrix)-1)) {
    temp.result <- evalOneMonoSig(global.df100, frac.df, 
                                  sig.matrix[,c(i,ncol(sig.matrix))],
                                  BeatAML, types[i], gmt)
    temp.wv.df <- temp.result$wv
    temp.test.df <- temp.result$test
    temp.drug.corr.df <- temp.result$drug.corr
    temp.frac.corr.df <- temp.result$frac.corr
    temp.wv.df$Signature <- names(sig.matrix)[i]
    temp.test.df$Signature <- names(sig.matrix)[i]
    temp.drug.corr.df$Signature <- names(sig.matrix)[i]
    temp.frac.corr.df$Signature <- names(sig.matrix)[i]
    
    wv.df <- rbind(wv.df, temp.wv.df)
    test.df <- rbind(test.df, temp.test.df)
    drug.corr.df <- rbind(drug.corr.df, temp.drug.corr.df)
    frac.corr.df <- rbind(frac.corr.df, temp.frac.corr.df)
  }
  p.df <- test.df[test.df$variable == "p.value",]
  
  return(list(wv = wv.df, p = p.df, 
              drug.corr = drug.corr.df, frac.corr = frac.corr.df))
}

compareSigs <- function(global.df100, frac.df, sigs, 
                        value.var = "Log2FC", BeatAML, types=rep("rna",length(sigs)), gmt) {
  # combine signatures into matrix
  filtered.sigs.df <- data.table::rbindlist(sigs, use.names = TRUE, idcol = "Signature")
  sig.matrix <- reshape2::dcast(filtered.sigs.df, Gene ~ Signature, mean,
                                value.var = value.var)
  rownames(sig.matrix) <- sig.matrix$Gene
  sig.matrix$Gene <- NULL
  
  # test signature matrix
  sigResults <- evalMonoSig(global.df100, frac.df, sig.matrix, BeatAML,types,gmt)
  wv.df <- sigResults$wv
  p.df <- sigResults$p
  drug.corr.df <- sigResults$drug.corr
  frac.corr.df <- sigResults$frac.corr
  
  # plot results
  p.df$Significance <- "p > 0.05"
  p.df$value <- as.numeric(p.df$value)
  if (any(p.df$value <= 0.05)) {
    p.df[p.df$value <= 0.05,]$Significance <- "p <= 0.05" 
  }
  p.df$Significance <- factor(p.df$Significance,
                                      levels=c("p <= 0.05", "p > 0.05"))
  p.df$minusLogFDR <- 1E-4
  if (any(p.df$value != 0)) {
    p.df[p.df$value != 0,]$minusLogFDR <- -log(p.df$value, base=10) 
  }
  sigOrder <- p.df[order(p.df$value),]$Signature
  ggplot(p.df, aes(x=Signature, y=-log(value, base=10), 
                   fill = Significance)) + 
    geom_col() + theme_classic() + ylab("-Log(P-value)") + 
    ggplot2::scale_x_discrete(limits = sigOrder) +
    scale_fill_manual(values=c("red","grey"), breaks=c("p <= 0.05", "p > 0.05"))+
    ggtitle("T-test: Monocyte scores are higher in CD14+ samples")
  ggsave("pValue_DIA_WV.pdf", width = 5, height = 5)
  
  if ("Drug" %in% colnames(drug.corr.df)) {
    doi <- c("Azacytidine", "Venetoclax", "Azacytidine - Venetoclax")
    drug.corr.df$`Drug Treatment` <- NA
    drug.corr.df[drug.corr.df$Drug == "Azacytidine",]$`Drug Treatment` <- "Aza"
    drug.corr.df[drug.corr.df$Drug == "Azacytidine - Venetoclax",]$`Drug Treatment` <- "Aza + Ven"
    drug.corr.df[drug.corr.df$Drug == "Venetoclax",]$`Drug Treatment` <- "Ven"
    drug.corr.df$Significance <- "Adjusted p > 0.05"
    if (any(drug.corr.df$Pearson.q <= 0.05)){
      drug.corr.df[drug.corr.df$Pearson.q <= 0.05,]$Significance <- "Adjusted p <= 0.05" 
    }
    drug.corr.df$Significance <- factor(drug.corr.df$Significance,
                                        levels=c("Adjusted p <= 0.05", "Adjusted p > 0.05"))
    plot.df <- drug.corr.df[drug.corr.df$Drug %in% doi,]
    av.df <- drug.corr.df[drug.corr.df$Drug == "Azacytidine - Venetoclax",]
    #sigOrder <- unique(plot.df[order(plot.df$Pearson.est, decreasing=TRUE),]$Signature)
    sigOrder <- av.df[order(av.df$Pearson.est, decreasing=TRUE),]$Signature
    ggplot(plot.df, aes(x=Signature, y=Pearson.est, fill = Significance)) + 
      geom_col() + theme_minimal() + ylab("Pearson Correlation Estimate") + 
      facet_wrap(~ `Drug Treatment`) +
      ggplot2::scale_x_discrete(limits = sigOrder) +
      theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
      scale_fill_manual(values=c("red","grey"), breaks=c("Adjusted p <= 0.05", "Adjusted p > 0.05"))+
      ggtitle("Monocytic signatures predict drug response")
    ggsave("drugCorr_DIA_WV.pdf", width = 5, height = 5)
    
    plot.df <- drug.corr.df[drug.corr.df$Drug == "Venetoclax",]
    sigOrder <- unique(plot.df[order(plot.df$Pearson.est, decreasing=TRUE),]$Signature)
    ggplot(plot.df, aes(x=Signature, y=Pearson.est, fill = Significance)) + 
      geom_col() + theme_classic() + ylab("Pearson Correlation Estimate") + 
      ggplot2::scale_x_discrete(limits = sigOrder) +
      scale_fill_manual(values=c("red","grey"), breaks=c("Adjusted p <= 0.05", "Adjusted p > 0.05"))+
      ggtitle("Monocytic signatures predict Venetoclax response")
    ggsave("VenCorr_DIA_WV.pdf", width = 5, height = 5)
    
    plot.df <- drug.corr.df[drug.corr.df$Drug == "Azacytidine - Venetoclax",]
    sigOrder <- unique(plot.df[order(plot.df$Pearson.est, decreasing=TRUE),]$Signature)
    ggplot(plot.df, aes(x=Signature, y=Pearson.est, fill = Significance)) + 
      geom_col() + theme_classic() + ylab("Pearson Correlation Estimate") + 
      ggplot2::scale_x_discrete(limits = sigOrder) +
      scale_fill_manual(values=c("red","grey"), breaks=c("Adjusted p <= 0.05", "Adjusted p > 0.05"))+
      ggtitle("Monocytic signatures predict Aza + Ven response")
    ggsave("AzaVenCorr_DIA_WV.pdf", width = 5, height = 5)
    
    plot.df <- drug.corr.df[drug.corr.df$Drug == "Azacytidine",]
    sigOrder <- unique(plot.df[order(plot.df$Pearson.est, decreasing=TRUE),]$Signature)
    ggplot(drug.corr.df[drug.corr.df$Drug == "Azacytidine",], 
           aes(x=Signature, y=Pearson.est, 
               fill = Significance)) + 
      geom_col() + theme_classic() + ylab("Pearson Correlation Estimate") + 
      ggplot2::scale_x_discrete(limits = sigOrder) +
      scale_fill_manual(values=c("red","grey"), breaks=c("Adjusted p <= 0.05", "Adjusted p > 0.05"))+
      ggtitle("Monocytic signatures predict Azacytidine response")
    ggsave("AzaCorr_DIA_WV.pdf", width = 5, height = 5)
  }
  
  frac.corr.df$Significance <- "p > 0.05"
  if (any(frac.corr.df$Pearson.p <= 0.05)) {
    frac.corr.df[frac.corr.df$Pearson.p <= 0.05,]$Significance <- "p <= 0.05" 
  }
  frac.corr.df$Significance <- factor(frac.corr.df$Significance,
                                      levels=c("p <= 0.05", "p > 0.05"))
  sigOrder <- frac.corr.df[order(frac.corr.df$Pearson.est, decreasing=TRUE),]$Signature
  ggplot(frac.corr.df, aes(x=Signature, y=Pearson.est,
                           fill = Significance)) + 
    geom_col() + theme_classic() + ylab("Pearson Correlation Estimate") + 
    ggplot2::scale_x_discrete(limits = sigOrder) +
    scale_fill_manual(values=c("red","grey"), breaks=c("p <= 0.05", "p > 0.05"))+
    ggtitle("Monocytic signatures predict monocyte fraction")
  ggsave("fracCorr_DIA_WV_Pearson.pdf", width = 5, height = 5)
  
  frac.corr.df$Significance <- "p > 0.05"
  if (any(frac.corr.df$Spearman.p <= 0.05)) {
    frac.corr.df[frac.corr.df$Spearman.p <= 0.05,]$Significance <- "p <= 0.05" 
  }
  frac.corr.df$Significance <- factor(frac.corr.df$Significance,
                                      levels=c("p <= 0.05", "p > 0.05"))
  sigOrder <- frac.corr.df[order(frac.corr.df$Spearman.est, decreasing=TRUE),]$Signature
  ggplot(frac.corr.df, aes(x=Signature, y=Spearman.est,
                           fill = Significance)) + 
    geom_col() + theme_classic() + ylab("Spearman Correlation Estimate") + 
    ggplot2::scale_x_discrete(limits = sigOrder) +
    scale_fill_manual(values=c("red","grey"), breaks=c("p <= 0.05", "p > 0.05"))+
    ggtitle("Monocytic signatures predict monocyte fraction")
  ggsave("fracCorr_DIA_WV_Spearman.pdf", width = 5, height = 5)
  
  return(list(wv = wv.df, p = p.df, drug.corr = drug.corr.df, 
              frac.corr = frac.corr.df))
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
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
setwd("data")
sorted.patients <- c("18-00105", "21-00839", "22-00571", "22-00117", "16-01184",
                     "19-00074", "18-00103", "21-00432", "17-01060", "22-00251")
BeatAML <- load_not_norm_BeatAML_for_DMEA3(exclude.samples = sorted.patients)
gmt.drug <- readRDS("gmt_BeatAML_drug_MOA_2024-02-22.rds")

# load sorted proteomics
synapser::synLogin()
global.df <- read.csv(synapser::synGet("syn58895933")$path) # DIA
global.df <- global.df[ , which(colMeans(!is.na(global.df)) >= 0.75)] # 36 out of 48 samples are kept
outliers <- c("X00839_CD34plusFlow", "X00117_CD34plus", 
              "X00432_CD14plus", "X00251_CD14plus", "X00105_CD14plusFlow")
# syn.test <- synapser::synGet("syn58914135") # for some reason, can't access pre-filtered version ???
rownames(global.df) <- global.df$Gene
global.df$Gene <- NULL
global.df <- as.data.frame(t(global.df))
global.df$Sample <- rownames(global.df)
global.df <- global.df[,c("Sample", colnames(global.df)[1:(ncol(global.df)-1)])]
global.df <- global.df[!(global.df$Sample %in% outliers),]
global.df100 <- global.df[,colSums(is.na(global.df)) == 0]
global.df100 <- global.df100[!grepl("flow",global.df100$Sample, ignore.case=TRUE),] # 17 samples

# load sorted proteomics signature
sig.paths <- list("Sorted" = "analysis/DIA_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_noNA.csv",
                  "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
                  "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HSC_protein-coding.csv")
#cd14.sig <- na.omit(read.csv(synapser::synGet("syn64543462")$path)) # DIA

# import signatures and filter
sigs <- list()
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
for (i in names(sig.paths)) {
  sigs[[i]] <- read.csv(sig.paths[[i]])
  sigs[[i]] <- na.omit(sigs[[i]][sigs[[i]]$adj.P.Val <= 0.05,c("Gene","Log2FC")])
}

# load cell fraction data
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/data/clinical_metadata/"
setwd(base.path)
patient.key <- readxl::read_xlsx("lab_dbgap_key.xlsx")
patient.meta <- readxl::read_xlsx("Table_S1.xlsx",sheet=1)
patient.meta <- as.data.frame(patient.meta)
#rownames(patient.meta) <- patient.meta[,1]

# does metadata correlate Ven response?
key.cols <- colnames(patient.meta)[colnames(patient.meta) %in% colnames(patient.key)]
num.cols <- colnames(dplyr::select_if(patient.meta, is.numeric))
num.meta <- dplyr::distinct(patient.meta[,c(key.cols, num.cols)])
num.meta <- merge(patient.key, num.meta, by = key.cols)
num.meta <- num.meta[,c("labId", num.cols[2:length(num.cols)])] # leave out dbgap_subject_id from numeric columns
frac.meta <- num.meta[,c("labId", "%.Monocytes.in.PB")]
frac.meta$labId <- sub(".*-","X",frac.meta$labId)
sortedIDs <- unique(sub("_.*","",global.df100$Sample))
frac.meta2 <- frac.meta[frac.meta$labId %in% sortedIDs,] # 7
frac.meta2 <- na.omit(frac.meta2) # 4

gmt <- readRDS("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/data/gmt_BeatAML_drug_MOA_2024-02-22.rds")
setwd("/Users/gara093/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create("Monocyte_vs_progenitor_signatures_beadOnly")
setwd("Monocyte_vs_progenitor_signatures_beadOnly")

# evaluate signature
evalResults <- compareSigs(global.df100, frac.meta2, sigs, BeatAML = BeatAML, types=c("global", "rna", "rna", "rna"), gmt = gmt.drug) 
write.csv(evalResults$wv, "wv.csv", row.names = FALSE)
write.csv(evalResults$p, "pValues.csv", row.names = FALSE)
write.csv(evalResults$drug.corr, "drugCorrelations.csv", row.names = FALSE)
write.csv(evalResults$frac.corr, "cellFractionCorrelations.csv", row.names = FALSE)
