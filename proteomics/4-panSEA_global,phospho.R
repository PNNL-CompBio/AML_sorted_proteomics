# Differential expression & enrichment analyses: global & phospho
# PTRC2: exp 24 TMT
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-02-09

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr); library(dplyr); library(Biobase)

setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
source("panSEA_helper_20240508.R")

# overview
# 1. Import metadata & crosstabs
# 2. Import BeatAML data formatted for DMEA
# 3. Run panSEA across omics for each contrast

#### 1. Import metadata & crosstabs ####
### TMT
setwd("data")
meta.df <- readxl::read_excel("Exp24metadataTable_TMT.xlsx") 
meta.df$MeasurementName <- as.character(rownames(meta.df))
rownames(meta.df) <- paste0("X", rownames(meta.df))
meta.df$row.name <- rownames(meta.df)

# add other drug info & make sure sensitivity is correctly labeled
sens.info <- read.csv("Exp24_drug_sensitivity_20240209.csv")
meta.df <- merge(meta.df, sens.info)
meta.df$X <- NULL
# add 'X' to MeasurementName to match colnames of data
meta.df$id <- paste0('X', meta.df$MeasurementName)

# add other metadata for contrasts
meta.df$CD14 <- "Neg"
meta.df[grepl("CD14+", meta.df$SampleType),]$CD14_Pos <- "Neg"

meta.df$CD34 <- "Neg"
meta.df[grepl("CD34+", meta.df$SampleType),]$CD34_Pos <- "Pos"

meta.df$MSC <- "Non_MSC"
meta.df[meta.df$SampleType == "MSC Flow",]$MSC <- "MSC"

meta.df$'Sort Type' <- "Bead"
meta.df[grepl("Flow", meta.df$SampleType),]$'Sort Type' <- "Flow"

meta.df$Pooled <- meta.df$SampleType
meta.df[meta.df$CD14=="Pos",]$Pooled <- "CD14+"
meta.df[meta.df$CD34=="Pos",]$Pooled <- "CD34+"
meta.df[meta.df$MSC=="MSC",]$Pooled <- "MSC"

global.df <- read.table(
  "global_data/Exp24_crosstab_global_gene_corrected.txt", 
  sep = "\t")
phospho.df <- read.table(
  "phospho_data/Exp24_crosstab_phospho_SiteID_corrected.txt", 
  sep = "\t")

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)
phospho.df$SUB_SITE <- rownames(phospho.df)
tmt <- list("meta" = meta.df,
            "global" = global.df, # 5169 gene symbols
            "phospho" = phospho.df) # 794 phospho-sites

### DIA
meta.df <- readxl::read_excel("Exp24metadataTable_DIA.xlsx") 
meta.df$id <- stringr::str_split_i(meta.df$patient, "-", -1)
meta.df$id <- paste0("X", meta.df$id)
meta.df[meta.df$`sample type` == "CD14+", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD14+", ]$id, "_CD14plus")
meta.df[meta.df$`sample type` == "CD34+", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD34+", ]$id, "_CD34plus")
meta.df[meta.df$`sample type` == "CD14+ Flow", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD14+ Flow", ]$id, "_CD14plusFlow")
meta.df[meta.df$`sample type` == "CD34+ Flow", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD34+ Flow", ]$id, "_CD34plusFlow")
meta.df[meta.df$`sample type` == "MSC Flow", ]$id <- paste0(meta.df[meta.df$`sample type` == "MSC Flow", ]$id, "_MSCflow")
rownames(meta.df) <- meta.df$id

# add other drug info & make sure sensitivity is correctly labeled
meta.df <- merge(meta.df, sens.info)
meta.df$X <- NULL

# add other metadata for contrasts
meta.df$CD14 <- "Neg"
meta.df[grepl("CD14+", meta.df$`sample type`),]$CD14_Pos <- "Neg"

meta.df$CD34 <- "Neg"
meta.df[grepl("CD34+", meta.df$`sample type`),]$CD34_Pos <- "Pos"

meta.df$MSC <- "Non_MSC"
meta.df[meta.df$`sample type` == "MSC Flow",]$MSC <- "MSC"

meta.df$'Sort Type' <- "Bead"
meta.df[grepl("Flow", meta.df$`sample type`),]$'Sort Type' <- "Flow"

meta.df$Pooled <- meta.df$SampleType
meta.df[meta.df$CD14=="Pos",]$Pooled <- "CD14+"
meta.df[meta.df$CD34=="Pos",]$Pooled <- "CD34+"
meta.df[meta.df$MSC=="MSC",]$Pooled <- "MSC"

synapser::synLogin()
globalFile <- synapser::synGet("syn55234888") # Samantha's unprocessed version
global.df <- read.table(
  globalFile$path, 
  sep = "\t") # 8897 proteins, 48 samples

# require proteins to be quantified in at least half of samples
global.df <- global.df[which(rowSums(is.na(global.df)) < ncol(global.df)/2),] # 6115 proteins, 48 samples

# require samples to have at least 50% of proteins quantified
#global.df <- global.df[ , which(colSums(is.na(global.df)) < nrow(global.df)/2)] # 42 out of 48 samples are kept
global.df75 <- global.df[ , which(colMeans(!is.na(global.df)) >= 0.75)] # 36 out of 48 samples are kept
global.df <- NULL

# log2-transform DIA data
#global.df <- log(global.df, 2)

global.df75 <- log(global.df75[,colnames(global.df75)!="Gene"],2)

# subtract row (protein) medians
#global_row_medians <- apply(global.df, 1, median, na.rm = T)
#global.df <- sweep(global.df, 1, global_row_medians, FUN = '-')

global_row_medians75 <- apply(global.df75, 1, median, na.rm = T)
global.df75 <- sweep(global.df75, 1, global_row_medians75, FUN = '-')

# subtract column (sample) medians
#global_sample_coef <- apply(global.df, 2, median, na.rm = T)
#global.df <- sweep(global.df, 2, global_sample_coef, FUN = '-')

global_sample_coef75 <- apply(global.df75, 2, median, na.rm = T)
global.df75 <- sweep(global.df75, 2, global_sample_coef75, FUN = '-')

# add column for feature names and later make it the first column
#global.df$Gene <- rownames(global.df) # 6092 gene symbols
global.df75$Gene <- rownames(global.df75) # 6092 gene symbols
#write.csv(global.df, "Exp24_DIA_crosstab_global_gene_corrected.csv", row.names = FALSE)
#write.csv(global.df75, "Exp24_DIA_75PercentCoverage_crosstab_global_gene_corrected.csv", row.names = FALSE)
# dia <- list("meta" = meta.df[meta.df$id %in% colnames(global.df), ],
#             "global" = global.df)
dia75 <- list("meta" = meta.df[meta.df$id %in% colnames(global.df75), ],
            "global" = global.df75)
#synfile <- synapser::File("Exp24_DIA_crosstab_global_gene_corrected.csv", "syn54821995")
#synapser::synStore(synfile)

### combine DIA & TMT
#dia$meta$method <- "DIA"
dia75$meta$method <- "DIA"
tmt$meta$method <- "TMT"

# make sure dia & tmt meta data have same columns
#dia$meta[, colnames(dia$meta)[!(colnames(dia$meta) %in% colnames(tmt$meta))]] <- NULL
tmt$meta[, colnames(tmt$meta)[!(colnames(tmt$meta) %in% colnames(dia75$meta))]] <- NULL
dia75$meta[, colnames(dia75$meta)[!(colnames(dia75$meta) %in% colnames(tmt$meta))]] <- NULL

dia.tmt75 <- list("meta" = rbind(dia75$meta, tmt$meta),
                "global" = merge(dia75$global, tmt$global, all = TRUE, 
                                 suffixes = c("_DIA", "_TMT")))

#### 2. Histograms and PCA ####
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
setwd(base.path)

phenos <- c("Plex", "id", "patient", "SampleType", "Pooled", "Aza", "Ven", "Aza.Ven", "Flow", "method")
library(MSnSet.utils)
setwd(base.path)

if (file.exists("gmt_BeatAML_drug_MOA.rds")) {
  gmt.drug <- readRDS("gmt_BeatAML_drug_MOA.rds")
} else {
  gmt.drug <- DMEA::as_gmt(moa.BeatAML, sep = ", ")
  saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA.rds")
}

dir.create("histograms_and_PCA")
setwd("histograms_and_PCA")
dir.create("DIAnotSampleCentered")
setwd("DIAnotSampleCentered")
for (k in 1:length(method.data)) {
  if (grepl("DIA", names(method.data)[k])) {
    omics <- list("Global" = method.data[[k]]$global)
  } else {
    omics <- list("Global" = method.data[[k]]$global,
                  "Phospho" = method.data[[k]]$phospho)
  }
  temp.meta <- method.data[[k]]$meta
  temp.phenos <- phenos[phenos %in% colnames(temp.meta)]
  
  for (i in 1:length(omics)) {
    # histogram
    melted.df <- reshape2::melt(omics[[i]])
    xlab <- paste("Normalized", names(omics)[i], "Expression")
    title <- names(method.data)[k]
    pdf(file.path(paste0(names(method.data)[k], "_", names(omics)[i], "_histogram_", Sys.Date(), ".pdf")))
    hist(melted.df$value, xlab = xlab, main = title)
    dev.off()
    
    # pca
    if (names(method.data)[k] != "DIA_&_TMT") {
      #sample.names <- temp.meta$id[temp.meta$id %in% colnames(omics[[i]])]
      sample.names <- colnames(omics[[i]])[colnames(omics[[i]]) %in% temp.meta$id]
      pca.data <- MSnSet(exprs = omics[[i]][, sample.names] %>% as.matrix(),
                         pData = temp.meta[sample.names,])
      for (j in 1:length(temp.phenos)) {
        MSnSet.utils::plot_pca(pca.data, phenotype = temp.phenos[j], label = "patient") + ggtitle(paste(names(method.data)[k], names(omics)[i], "PCA"))
        ggsave(paste0(names(method.data)[k], "_", names(omics)[i], "_PCA_", temp.phenos[j], "_", Sys.Date(), ".pdf")) # 3018 DIA (49%), 3975 (77%) TMT; wo outliers1: 3261 DIA (53%), 3986 TMT (77%), 17 TMT phospho (2%)
      } 
    }
  }
}
# wo outliers2: 3368 DIA (55%), 4793 TMT (93%), 85 TMT phospho (11%)
# wo outliers2a: 3408 DIA (56%), 4860 TMT (94%), 122 TMT phospho (15%)
# not sample centered DIA, no outliers removed: 3018 DIA (49%), 3975 (77%) TMT, 14 TMT phospho (2%)
# 1917 complete rows for DIA, 3975 complete rows for TMT

# determine outliers based on PCAs and iterate
outliers <- c("X00839_CD34plusFlow", "X00117_CD34plus", 
              "X00432_CD14plus", "X00251_CD14plus")
outlier.ids <- c("X11", "X17", "X7", "X23")
possible.outlier.id <- c("X28")
possible.outlier <- c("X00839_CD14plus")

dia.wo.out <- list("meta" = dia75$meta[!(dia75$meta$id %in% outliers),],
                   "global" = dia75$global[,colnames(dia75$global)[!(colnames(dia75$global) %in% outliers)]]) #6115 proteins (same), 32 samples down from 36
tmt.wo.out <- list("meta" = tmt$meta[!(tmt$meta$id %in% outlier.ids),],
                   "global" = tmt$global[,colnames(tmt$global)[!(colnames(tmt$global) %in% outlier.ids)]],
                   "phospho" = tmt$phospho[,colnames(tmt$phospho)[!(colnames(tmt$phospho) %in% outlier.ids)]]) #5169 proteins (same), 794 phospho sites (same), 24 samples down from 28
dia.tmt.wo.out <- list("meta" = rbind(dia.wo.out$meta, tmt.wo.out$meta),
                       "global" = merge(dia.wo.out$global, tmt.wo.out$global, all = TRUE, 
                                        suffixes = c("_DIA", "_TMT")))
method.data <- list("DIA" = dia.wo.out,
                    "DIA_&_TMT" = dia.tmt.wo.out,
                    "TMT" = tmt.wo.out)

outliers2 <- c("X00103_CD14plusFlow", "X00103_CD14plus", "X00839_CD14plus", "X01184_MSCflow")
outlier.ids2 <- c("X4", "X28", "X30")
dia.wo.out2 <- list("meta" = dia.wo.out$meta[!(dia.wo.out$meta$id %in% outliers2),],
                   "global" = dia.wo.out$global[,colnames(dia.wo.out$global)[!(colnames(dia.wo.out$global) %in% outliers2)]]) #6115 proteins (same), 29 samples down from 36
tmt.wo.out2 <- list("meta" = tmt.wo.out$meta[!(tmt.wo.out$meta$id %in% outlier.ids2),],
                   "global" = tmt.wo.out$global[,colnames(tmt.wo.out$global)[!(colnames(tmt.wo.out$global) %in% outlier.ids2)]],
                   "phospho" = tmt.wo.out$phospho[,colnames(tmt.wo.out$phospho)[!(colnames(tmt.wo.out$phospho) %in% outlier.ids2)]]) #5169 proteins (same), 794 phospho sites (same), 21 samples down from 28
dia.tmt.wo.out2 <- list("meta" = rbind(dia.wo.out$meta, tmt.wo.out$meta),
                       "global" = merge(dia.wo.out$global, tmt.wo.out$global, all = TRUE, 
                                        suffixes = c("_DIA", "_TMT")))
method.data <- list("DIA" = dia.wo.out2,
                    "DIA_&_TMT" = dia.tmt.wo.out2,
                    "TMT" = tmt.wo.out2)

# missed 18-00105 CD14+ Flow
outliers2a <- c("X00103_CD14plusFlow", "X00103_CD14plus", "X00839_CD14plus", 
                "X01184_MSCflow", "X00105_CD14plusFlow")
outlier.ids2a <- c("X4", "X28", "X30", "X26")
dia.wo.out2a <- list("meta" = dia.wo.out$meta[!(dia.wo.out$meta$id %in% outliers2a),],
                    "global" = dia.wo.out$global[,colnames(dia.wo.out$global)[!(colnames(dia.wo.out$global) %in% outliers2a)]]) #6115 proteins (same), 29 samples down from 36
tmt.wo.out2a <- list("meta" = tmt.wo.out$meta[!(tmt.wo.out$meta$id %in% outlier.ids2a),],
                    "global" = tmt.wo.out$global[,colnames(tmt.wo.out$global)[!(colnames(tmt.wo.out$global) %in% outlier.ids2a)]],
                    "phospho" = tmt.wo.out$phospho[,colnames(tmt.wo.out$phospho)[!(colnames(tmt.wo.out$phospho) %in% outlier.ids2a)]]) #5169 proteins (same), 794 phospho sites (same), 21 samples down from 28
dia.tmt.wo.out2a <- list("meta" = rbind(dia.wo.out2a$meta, tmt.wo.out2a$meta),
                        "global" = merge(dia.wo.out2a$global, tmt.wo.out2a$global, all = TRUE, 
                                         suffixes = c("_DIA", "_TMT")))
method.data <- list("DIA" = dia.wo.out2a,
                    "TMT" = tmt.wo.out2a)

#### 3. Run panSEA for each contrast & combination ####
synapse_id <- "syn53606820"
all.degs <- data.frame()
contrasts <- c("CD14", "CD34", "MSC", "Aza", "Ven", "Aza.Ven", "Sort Type")
for (k in 1:length(method.data)) {
  setwd(base.path)
  method.path <- file.path(base.path, names(method.data)[k])
  dir.create(names(method.data)[k])
  setwd(names(method.data)[k])
  methodFolder <- 
    synapser::synStore(synapser::Folder(names(method.data)[k],
                                        parent = synapse_id))
  
  meta.df <- method.data[[k]]$meta
  
  # run contrast combos
  if (grepl("DIA", names(method.data)[k])) {
    omics <- list("global" = method.data[[k]]$global)
    feature.names <- "Gene"
    temp.expr <- list(BeatAML.data$global)
  } else {
    omics <- list("global" = method.data[[k]]$global,
                  "phospho" = method.data[[k]]$phospho)
    feature.names <- c("Gene", "SUB_SITE")
    temp.expr <- list(BeatAML.data$global, BeatAML.data$phospho)
  }
  names(temp.expr) <- names(omics)
  panSEA2_combos(contrasts, meta.df = meta.df, 
                      omics = omics,
                      expr = temp.expr,
                      gmt.drug = gmt.drug, drug.sens = BeatAML.data$drug,
                      base.path = base.path,
                      temp.path = method.path,
                      synapse_id = methodFolder)
  
  # get compiled DEGs
  methodDEGs <- as.list(synapser::synGetChildren(methodFolder, list("file"), sortBy = 'NAME'))
  if (length(methodDEGs) > 0) {
    if (methodDEGs[[1]]$name == "Differential_expression_results.csv") {
      methodFile <- synapser::synGet(methodDEGs[[1]]$id)
      methodDEG <- read.csv(methodFile$path)
      methodDEG$method <- names(method.data)[k]
      all.degs <- rbind(all.degs, methodDEG)
    }
  }
}
setwd(base.path)
all.DEG.files <- list("Differential_expression_results.csv" = 
                        all.degs,
                      "Differential_expression_results_max_5_percent_FDR.csv" = 
                        all.degs[all.degs$adj.P.Val <= 0.05, ])
save_to_synapse(all.DEG.files, synapse_id)

#### 4. Cell type predictions using weighted voting of CD14+ vs. CD34+ ####
#### CD14 vs CD34: same as MSC_Non_MSC: CD14_Pos_vs_Neg
#### Weighted voting: can CD14 vs CD34 signature rank data when not normalized across samples?
### prep raw DIA
synapser::synLogin()
globalFileDIA <- synapser::synGet("syn55234888") # Samantha's unprocessed version
global.df.DIA <- read.table(
  globalFileDIA$path, 
  sep = "\t") # 8897 proteins, 48 samples

# require proteins to be quantified in at least half of samples
global.df.DIA <- global.df.DIA[which(rowSums(is.na(global.df.DIA)) < ncol(global.df.DIA)/2),] # 6115 proteins, 48 samples

# require samples to have at least 75% of proteins quantified
global.df.DIA <- global.df.DIA[ , which(colMeans(!is.na(global.df.DIA)) >= 0.75)] # 36 out of 48 samples are kept

# remove outliers
global.df.DIA <- global.df.DIA[,colnames(global.df.DIA)[!(colnames(global.df.DIA) %in% c(outliers, outliers2a))]]

# log2-transform DIA data
global.df.DIA <- log(global.df.DIA[,colnames(global.df.DIA)!="Gene"],2)

# subtract column (sample) medians
global_sample_coef.DIA <- apply(global.df.DIA, 2, median, na.rm = T)
global.df.DIA <- sweep(global.df.DIA, 2, global_sample_coef.DIA, FUN = '-')

# transform to prep for weighted voting
global.df.DIA.trans <- t(global.df.DIA)

### prep raw TMT
globalFileTMT <- synapser::synGet("syn53493077") # unprocessed version
global.df.TMT <- read.table(
  globalFileTMT$path, 
  sep = "\t") # 8897 proteins, 48 samples

# require proteins to be quantified in at least half of samples
global.df.TMT <- global.df.TMT[which(rowSums(is.na(global.df.TMT)) < ncol(global.df.TMT)/2),] # 6115 proteins, 48 samples

# require samples to have at least 75% of proteins quantified
global.df.TMT <- global.df.TMT[ , which(colMeans(!is.na(global.df.TMT)) >= 0.75)] # 36 out of 48 samples are kept

# remove outliers
global.df.TMT <- global.df.TMT[,colnames(global.df.TMT)[!(colnames(global.df.TMT) %in% c(outlier.ids, outlier.ids2a))]]

# subtract column (sample) meTMTns
global_sample_coef.TMT <- apply(global.df.TMT, 2, median, na.rm = T)
global.df.TMT <- sweep(global.df.TMT, 2, global_sample_coef.TMT, FUN = '-')

# transform to prep for weighted voting
global.df.TMT.trans <- t(global.df.TMT)

### load signature of CD14+ vs. CD34+ cells
synapser::synLogin()
globalFileDIA <- synapser::synGet("syn59429685") # filtered for max FDR of 0.05
global.sig.DIA <- read.csv(globalFileDIA$path) # 1842 proteins

globalFileTMT <- synapser::synGet("syn59438946") # filtered for max FDR of 0.05
global.sig.TMT <- read.csv(globalFileTMT$path) # 703 proteins

### add sample names
global.df.DIA.trans <- as.data.frame(global.df.DIA.trans)
global.df.DIA.trans$id <- rownames(global.df.DIA.trans)
global.df.TMT.trans <- as.data.frame(global.df.TMT.trans)
global.df.TMT.trans$id <- rownames(global.df.TMT.trans)
global.df.DIA.trans <- global.df.DIA.trans[,c("id", global.sig.DIA$Gene[global.sig.DIA$Gene %in% colnames(global.df.DIA.trans)])]
global.df.TMT.trans <- global.df.TMT.trans[,c("id", global.sig.TMT$Gene[global.sig.TMT$Gene %in% colnames(global.df.TMT.trans)])]

### only consider proteins which have 100% coverage
global.df.DIA.trans100 <- global.df.DIA.trans[,colSums(is.na(global.df.DIA.trans)) == 0] # 1120 proteins
global.df.TMT.trans100 <- global.df.TMT.trans[,colSums(is.na(global.df.TMT.trans)) == 0] # 679 proteins
global.df.DIA.trans100$id <- rownames(global.df.DIA.trans100)
global.df.DIA.trans100 <- global.df.DIA.trans100[,c("id", colnames(global.df.DIA.trans100)[colnames(global.df.DIA.trans100)!="id"])]
global.df.TMT.trans100$id <- rownames(global.df.TMT.trans100)
global.df.TMT.trans100 <- global.df.TMT.trans100[,c("id", colnames(global.df.TMT.trans100)[colnames(global.df.TMT.trans100)!="id"])]

### run WV
DIA.WV <- DMEA::WV(global.df.DIA.trans100, global.sig.DIA, sample.names = "id")
TMT.WV <- DMEA::WV(global.df.TMT.trans100, global.sig.TMT, sample.names = "id")

# merge with meta.df
WV.results <- rbind(DIA.WV$scores, TMT.WV$scores)
WV.df <- merge(WV.results, dia.tmt.wo.out2a$meta)
WV.df$Pooled <- "MSC"
WV.df[WV.df$Pooled_CD14_Pos, ]$Pooled <- "CD14+"
WV.df[WV.df$Pooled_CD34_Pos, ]$Pooled <- "CD34+"
WV.df$SampleType <- WV.df$Pooled
WV.df[WV.df$Flow, ]$SampleType <- paste(WV.df[WV.df$Flow, ]$SampleType, "Flow")
WV.df <- na.omit(WV.df)

# DIA p-values
WV.p.CD14.DIA <- t.test(WV.df[WV.df$Pooled_CD14_Pos & WV.df$method == "DIA",]$WV, WV.df[!WV.df$Pooled_CD14_Pos & WV.df$method == "DIA",]$WV)$p.value
WV.p.CD34.DIA <- t.test(WV.df[WV.df$Pooled_CD34_Pos & WV.df$method == "DIA",]$WV, WV.df[!WV.df$Pooled_CD34_Pos & WV.df$method == "DIA",]$WV)$p.value
WV.p.MSC.DIA <- t.test(WV.df[WV.df$MSC_Flow & WV.df$method == "DIA",]$WV, WV.df[!WV.df$MSC_Flow & WV.df$method == "DIA",]$WV)$p.value

# TMT p-values
WV.p.CD14.TMT <- t.test(WV.df[WV.df$Pooled_CD14_Pos & WV.df$method == "TMT",]$WV, WV.df[!WV.df$Pooled_CD14_Pos & WV.df$method == "TMT",]$WV)$p.value
WV.p.CD34.TMT <- t.test(WV.df[WV.df$Pooled_CD34_Pos & WV.df$method == "TMT",]$WV, WV.df[!WV.df$Pooled_CD34_Pos & WV.df$method == "TMT",]$WV)$p.value
WV.p.MSC.TMT <- t.test(WV.df[WV.df$MSC_Flow & WV.df$method == "TMT",]$WV, WV.df[!WV.df$MSC_Flow & WV.df$method == "TMT",]$WV)$p.value

### violin plot: WV score vs. Sample Type
setwd(base.path)
library(ggplot2)
marker.violin <- ggplot2::ggplot(WV.df, 
                                 aes(fill = method, x=Pooled, y=WV)) + 
  geom_violin(position=position_dodge(width=0.4), alpha=0.5) + 
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0.5) + 
  bg.theme3 + xlab("Sample Type") + ylab("CD14+ vs. CD34+ Score")
ggsave(paste0("CD14_vs_CD34_WV_100PercentCoverage_by_pooled_sample_type_", Sys.Date(), ".pdf"), marker.violin)
marker.violin <- ggplot2::ggplot(WV.df, 
                                 aes(fill = method, x=SampleType, y=WV)) + 
  geom_violin(position=position_dodge(width=0.4), alpha=0.5) + 
  geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0.5) + 
  bg.theme3 + xlab("Sample Type") + ylab("CD14+ vs. CD34+ Score")
ggsave(paste0("CD14_vs_CD34_WV_100PercentCoverage_by_sample_type_", Sys.Date(), ".pdf"), marker.violin)

#### 5. Correlations between DIA & TMT ####
corr.scatter <- function(df, rank.var, value, xlab, ylab, title, Pearson.p, Pearson.est, 
                         se = TRUE, position.x = "min", position.y = "max", 
                         shape = NULL, color = NULL, symmetrical = TRUE) {
  # load themes for plots
  ng.theme <- ggplot2::theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black"),
    #legend.title = element_blank(),
    axis.title.y = element_text(size = 8, colour = "black")
  )
  
  bg.theme <- ggplot2::theme(
    legend.background = element_rect(), legend.position = "right",
    legend.text = element_text(size = 14), 
    legend.title = element_text(size=16),
    legend.key = element_blank(),
    axis.title.x = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16),
    plot.title = element_text(lineheight = .8, face = "bold", size = 36)
  )
  
  # set plot parameters
  min.x <- min(df[, c(rank.var)])
  max.x <- max(df[, c(rank.var)])
  mid.x <- 0.5 * (min.x + max.x)
  min.y <- min(df[, c(value)])
  max.y <- max(df[, c(value)])
  mid.y <- 0.5 * (min.y + max.y)
  if (position.x == "min") {
    pos.x <- min.x
  } else if (position.x == "mid") {
    pos.x <- mid.x
  } else if (position.x == "max") {
    pos.x <- max.x
  } else if (is.numeric(position.x)) {
    pos.x <- position.x
  }
  if (position.y == "min") {
    pos.y <- min.y
  } else if (position.y == "mid") {
    pos.y <- mid.y
  } else if (position.y == "max") {
    pos.y <- max.y
  } else if (is.numeric(position.y)) {
    pos.y <- position.y
  }
  
  stats_pearson <- substitute(
    r == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(Pearson.est, digits = 3),
      p = format(Pearson.p, digits = 3)
    )
  )
  
  scatter.plot <- ggplot2::ggplot(data = df,
                                  aes_string(x = rank.var, y = value)) +
    ggplot2::geom_point(aes_string(
      shape = shape, color = color)) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::ggtitle(title) +
    ggplot2::geom_smooth(method = "lm", size = 1.5,
                         linetype = "solid", color = "blue",
                         se = se, na.rm = TRUE) +
    ggplot2::geom_text(
      x = pos.x, y = pos.y, vjust = "inward", hjust = "inward",
      colour = "blue", parse = TRUE,
      label = as.character(as.expression(stats_pearson)), size = 8
    ) +
    ng.theme +
    bg.theme + theme_minimal()
  
  if (symmetrical) {
    abs.max <- max(abs(c(min.x, max.x, min.y, max.y)))
    
    if (position.y == "min") {
      pos.y <- -abs.max
    } else if (position.y == "mid") {
      pos.y <- 0
    } else {pos.y <- abs.max}
    
    if (position.x == "min") {
      pos.x <- -abs.max
    } else if (position.x == "mid") {
      pos.x <- 0
    } else {pos.x <- abs.max}
    
    scatter.plot <- ggplot2::ggplot(data = df,
                                    aes_string(x = rank.var, y = value)) +
      ggplot2::geom_point(aes_string(
        shape = shape, color = color)) +
      ggplot2::labs(x = xlab, y = ylab) +
      ggplot2::ggtitle(title) +
      ggplot2::geom_smooth(method = "lm", size = 1.5,
                           linetype = "solid", color = "blue",
                           se = se, na.rm = TRUE) +
      ggplot2::geom_text(
        x = pos.x, y = pos.y, vjust = "inward", hjust = "inward",
        colour = "blue", parse = TRUE,
        label = as.character(as.expression(stats_pearson)), size = 8
      ) +
      ng.theme +
      bg.theme + theme_minimal() + xlim(c(-abs.max, abs.max)) + 
      ylim(c(-abs.max, abs.max)) + geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) + theme(axis.title.x = element_text(size = 20),
                                         axis.text.x = element_text(size = 16),
                                         axis.title.y = element_text(size = 20),
                                         axis.text.y = element_text(size = 16),
                                         plot.title = element_text(lineheight = .8, face = "bold", size = 36))
  }
  
  return(plot = scatter.plot)
}

#### raw DIA vs. TMT paired for each sample in both
# try giving TMT same sample names as DIA to correlate them
rownames(meta.df) <- meta.df$id
meta.wo.out2a <- meta.df[meta.df$id %in% colnames(global.df.TMT),]
meta.wo.out2a <- meta.wo.out2a[colnames(global.df.TMT),]
dia.names <- meta.wo.out2a$DIA_id
colnames(global.df.TMT) <- dia.names

# make sure we are comparing the same genes in the same orders
gene.names <- rownames(global.df.TMT)
gene.names <- gene.names[gene.names %in% rownames(global.df.DIA)]
matching.DIA <- global.df.DIA[gene.names,]
matching.TMT <- global.df.TMT[gene.names,]

# make sure we have the same samples in the same order in DIA & TMT
shared.dia.names <- dia.names[dia.names %in% colnames(matching.DIA)]
matching.DIA <- matching.DIA[,shared.dia.names]
matching.TMT <- matching.TMT[,shared.dia.names]
DIA.TMT.correlations <- corrr::correlate(matching.DIA, matching.TMT, use="pairwise.complete.obs", diagonal = 1)
rownames(DIA.TMT.correlations) <- DIA.TMT.correlations$term
heatmap(as.matrix(dplyr::select_if(DIA.TMT.correlations, is.numeric)))
pheatmap::pheatmap(as.matrix(dplyr::select_if(DIA.TMT.correlations, is.numeric)), cluster_rows = FALSE, cluster_cols = FALSE)

#### WV scores: DIA vs. TMT
library(ggplot2)
# prepare data frame
DIA.WV.results <- merge(DIA.WV$scores, meta.df, by.x="id", by.y = "DIA_id", suffixes = c("_DIA", "_TMT"))
WV.df <- merge(TMT.WV$scores, DIA.WV.results, by.x = "id", by.y = "id_TMT", suffixes = c("_TMT", "_DIA"))
WV.df$'Sample' <- WV.df$SampleType
WV.df$'Patient' <- WV.df$patient

# run correlation
corr.DIA.TMT <- cor.test(WV.df$WV_DIA, WV.df$WV_TMT)
p <- corr.DIA.TMT$p.value
est <- as.numeric(corr.DIA.TMT$estimate)

# create scatter plot
WV.DIA.TMT.plot <- corr.scatter(WV.df, "WV_DIA", "WV_TMT", "DIA", "TMT", 
                                "CD14+ vs. CD34+ Score", p, est, 
                                shape = "Sample", color = "Patient")
ggsave("DIA_vs_TMT_CD14_vs_CD34_score_v3.pdf", WV.DIA.TMT.plot, height=7, width=7)

#### differentially expressed proteins based on adjusted p-values <= 0.05: DIA vs. TMT
# prepare data frame
global.sig <- merge(global.sig.DIA, global.sig.TMT, by="Gene", suffixes=c("_DIA", "_TMT")) # 525

# run correlation
corr.DIA.TMT <- cor.test(global.sig$Log2FC_DIA, global.sig$Log2FC_TMT)
p <- corr.DIA.TMT$p.value
est <- as.numeric(corr.DIA.TMT$estimate)

# create scatter plot
sig.DIA.TMT.plot <- corr.scatter(global.sig, "Log2FC_DIA", "Log2FC_TMT", 
                                 "DIA", "TMT", "CD14+ vs. CD34+ Log2FC", p, est)
ggsave("DIA_vs_TMT_CD14_vs_CD34_Log2FC_v3.pdf", sig.DIA.TMT.plot, height=7, width=7)

#### Venn diagrams: DIA vs. TMT
### differentially expressed proteins based on adjusted p-values <= 0.05
venn.data <- list("TMT" = global.sig.TMT$Gene,
                  "DIA" = global.sig.DIA$Gene)
sig.venn <- ggvenn::ggvenn(venn.data, set_name_size = 8, text_size = 8)
ggsave("Venn_diagram_DIA_vs_TMT_CD14_vs_CD34_differential_expression.pdf", sig.venn, height=8, width=8)

### all proteins quantified after filters
all.venn.data <- list("TMT" = rownames(global.df.TMT), 
                      "DIA" = rownames(global.df.DIA))
all.venn <- ggvenn::ggvenn(all.venn.data, set_name_size = 8, text_size = 8)
ggsave("Venn_diagram_DIA_vs_TMT_proteins_quantified_after_filters.pdf", all.venn, height=8, width=8)

##### 6. Plot cell markers and DEGS #####
# plot markers
# MSCs (Mesenchymal Stem Cells): should be high in CD73, CD90, CD105, CD106, CD146, STRO-1 and low in CD14, CD34, CD45, HLA-DR
# AML monocytes should be high in CD4, CD11c, CD14, and CD64; source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5907644/#:~:text=Surface%20antigens%20that%20indicate%20leukemic,mostly%20expressed%20on%20mature%20monocytes.
# Leukemia stem cells (LSCs):  CD34, CD38, CD123, TIM3, CD25, CD32 and CD96; source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5388677/
# LSCs are emphasized in recent PTRC paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10950354/

markers <- unique(c("CD73", "CD90", "CD105", "CD106", "CD146", "STRO-1", "CD14", 
                    "CD34", "CD45", "HLA-DR", "CD4", "CD11c", "CD14", "CD64", 
                    "CD34", "CD38", "CD123", "TIM3", "CD25", "CD32", "CD96"))
markers <- unique(c("NT5E", "THY1", "ENG", "VCAM1", "MCAM", "CD14", "CD34", 
                    "PTPRC", "CD4", "ITGAX", "ITGAX", "FCGR1A", "CD38", "IL3RA", 
                    "HAVCR2", "IL2RA", "ISG20", "FCGR2A", "FCGR2B", "CD96"))
markers.from.Anupriya <- c("CD3", "HLA-DR", "CD1A", "CD4", "CD5", "ITGAL", 
                           "ITGAM", "CD14", "FUT4", "ITGB2", "CD19", "IL2RA",
                           "ISG20", "CD38", "NCAM1", "SELE", "SELP", "IL3RA",
                           "CDH5", "FASLG", "CD9", "ITGB1", "CD44", "CD46",
                           "CD47", "ITGA1", "ITGA2", "ITGA5", "CD58", "CD59",
                           "ITGB3", "CD63", "NT5E", "CD81", "THY1", "SLC3A2",
                           "SLC7A5", "BSG", "CD151", "CD200", "HLA-A", "HLA-B",
                           "HLA-C", "ITGA3", "ITGAV", "FAS", "ENG", "CD13", 
                           "ANPEP", "CD33")
all.markers <- unique(c(markers, markers.from.Anupriya))
dia.tmt.markers <- dia.tmt.wo.out2a$global[dia.tmt.wo.out2a$global$Gene %in% markers,
                                    c("Gene", colnames(dia.tmt.wo.out2a$global)[colnames(dia.tmt.wo.out2a$global) != "Gene"])]
write.csv(dia.tmt.markers, "Cell_markers_global_DIA_TMT_75percentCoverage_2024-05-04.csv", row.names = FALSE)

### violin plots of cell markers
library(ggplot2)
# load theme for plots
bg.theme3 <- ggplot2::theme(
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.border = element_rect(fill = NA),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"), 
  axis.title.x = element_text(size = 16, colour = "black"),
  axis.title.y = element_text(size = 16, colour = "black"),
  axis.text.x = element_text(colour = "black", size =16), 
  axis.text.y = element_text(colour = "black", size = 16),
  axis.ticks.x = element_line(colour = "black"), 
  axis.ticks.y = element_line(colour = "black"),
  legend.title = element_blank(), legend.background = element_rect(), 
  legend.position = "top",
  legend.text = element_text(size = 14), legend.key = element_blank(),
  plot.title = element_text(lineheight = .8, face = "bold", size = 36)
)

# prepare long data frame
long.global <- reshape2::melt(dia.tmt.wo.out2a$global, variable.name = "id")
long.global <- merge(long.global, dia.tmt.wo.out2a$meta, by = "id")

dir.create("markers")
setwd("markers")
for (i in 1:length(markers)) {
  marker.df <- long.global[long.global$Gene == markers[i],]
  if (nrow(marker.df) > 0) {
    marker.violin <- ggplot2::ggplot(marker.df, 
                                     aes(fill = method, x=Pooled, y=value)) + 
      geom_violin(position=position_dodge(width=0.4), alpha=0.5) + 
      geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0.5) + 
      bg.theme3 + xlab("Sample Type") + ylab("Normalized Protein Expression")
    ggsave(paste0(markers[i],"_75PercentCoverage_by_pooled_sample_type_", Sys.Date(), ".pdf"), marker.violin)
    marker.violin <- ggplot2::ggplot(marker.df, 
                                     aes(fill = method, x=SampleType, y=value)) + 
      geom_violin(position=position_dodge(width=0.4), alpha=0.5) + 
      geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0.5) + 
      bg.theme3 + xlab("Sample Type") + ylab("Normalized Protein Expression")
    ggsave(paste0(markers[i],"_75PercentCoverage_by_sample_type_", Sys.Date(), ".pdf"), marker.violin)
  }
}

### pathway data for heatmaps
# repeat for KRAS pathway hits based on Jeff's paper
jeff.markers <- c("KRAS", "PTPN11", "BCLXL", "MCL1", "CD40", "CD14", "CLEC7A", 
             "TRAF2", "IRAK1", "NFKB", "BCL2A1", "BCL2", "BAK", "BAX")

# Myc targets
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "H")
myc.targets <- unique(msigdb.info[grepl("MYC_TARGETS", msigdb.info$gs_name, 
                                        ignore.case = TRUE), ]$gene_symbol) # 240
myc.targetsv1 <- unique(msigdb.info[grepl("MYC_TARGETS_V1", msigdb.info$gs_name, 
                                          ignore.case = TRUE), ]$gene_symbol) # 200
myc.targetsv2 <- unique(msigdb.info[grepl("MYC_TARGETS_V2", msigdb.info$gs_name, 
                                          ignore.case = TRUE), ]$gene_symbol) # 58

# JAK/STAT
jak.stat <- unique(msigdb.info[grepl("JAK_STAT", msigdb.info$gs_name, 
                                     ignore.case = TRUE), ]$gene_symbol) 

# KRAS
kras <- unique(msigdb.info[grepl("KRAS", msigdb.info$gs_name, 
                                 ignore.case = TRUE), ]$gene_symbol) 

# TGF beta
msigdb.info <- msigdbr::msigdbr("Homo sapiens", "C2", "CP:KEGG")
tgf <- unique(msigdb.info[grepl("TGF_BETA", msigdb.info$gs_name, 
                                ignore.case = TRUE), ]$gene_symbol) 

genesets <- list('Myc_targets' = myc.targets,
                 'Myc_targets_v1' = myc.targetsv1,
                 'Myc_targets_v2' = myc.targetsv2,
                 "JAK-STAT" = jak.stat,
                 "KRAS" = kras,
                 "TGF_Beta" = tgf,
                 "Jeff" = jeff.markers)

# prep annotations for heatmaps
meta.df <- dia.tmt.wo.out2a$meta
cc.df <- meta.df[,c("Sort Type", "Sample Type", "Patient")]

# create heatmaps
DIA.heatmaps <- make_heatmaps(dia.wo.out2a$global, cc.df, top.gmt = genesets, fontsize=6)
TMT.heatmaps <- make_heatmaps(tmt.wo.out2a$global, cc.df, top.gmt = genesets, fontsize=6)
heatmaps <- list("Global_DIA" = DIA.heatmaps,
                 "Global_TMT" = TMT.heatmaps)

# save heatmaps
dir.create("heatmaps_of_interest")
setwd("heatmaps_of_interest")
dir.create("fontsize_6")
setwd("fontsize_6")
save_to_synapse(heatmaps)

# instead just look at individual data points
dir.create("markers")
setwd("markers")
for (i in 1:length(markers)) {
  marker.df <- na.omit(long.global[long.global$Gene == markers[i],])
  if (nrow(marker.df) > 0) {
    marker.violin <- ggplot2::ggplot(marker.df, 
                                     aes(color = patient, x=Pooled, y=value, shape = method)) + 
      geom_point(position=position_dodge(width=0.4)) +  
      bg.theme3 + xlab("Sample Type") + ylab("Normalized Protein Expression")
    ggsave(paste0(markers[i],"_75PercentCoverage_by_pooled_sample_type_dotplot_", Sys.Date(), ".pdf"), marker.violin, width = 11)
    marker.violin <- ggplot2::ggplot(marker.df, 
                                     aes(color = patient, x=SampleType, y=value, shape = method)) + 
      geom_point(position=position_dodge(width=0.4)) +  
      bg.theme3 + xlab("Sample Type") + ylab("Normalized Protein Expression")
    ggsave(paste0(markers[i],"_75PercentCoverage_by_sample_type_dotplot_", Sys.Date(), ".pdf"), marker.violin, width = 11)
  }
}

