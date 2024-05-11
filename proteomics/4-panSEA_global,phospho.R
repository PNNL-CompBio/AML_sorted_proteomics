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
#### 1. Import metadata & crosstabs ####
#### 2. Import BeatAML data formatted for DMEA ####
#### 3. Run panSEA across omics for each contrast ####

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
meta.df$Aza_Sensitive <- FALSE
meta.df[meta.df$Aza == "Sensitive",]$Aza_Sensitive <- TRUE

meta.df$Aza.Ven_Sensitive <- FALSE
meta.df[meta.df$Aza.Ven == "Sensitive",]$Aza.Ven_Sensitive <- TRUE

meta.df$Ven_Sensitive <- FALSE
meta.df[meta.df$Ven == "Sensitive",]$Ven_Sensitive <- TRUE

meta.df$Pooled_CD14_Pos <- FALSE
meta.df[grepl("CD14+", meta.df$SampleType),]$Pooled_CD14_Pos <- TRUE

meta.df$Pooled_CD34_Pos <- FALSE
meta.df[grepl("CD34+", meta.df$SampleType),]$Pooled_CD34_Pos <- TRUE

# meta.df$CD14_Pos <- FALSE
# meta.df[meta.df$SampleType == "CD14+",]$CD14_Pos <- TRUE
# 
# meta.df$CD34_Pos <- FALSE
# meta.df[meta.df$SampleType == "CD34+",]$CD34_Pos <- TRUE
# 
# meta.df$CD14_Pos_Flow <- FALSE
# meta.df[meta.df$SampleType == "CD14+ Flow",]$CD14_Pos_Flow <- TRUE
# 
# meta.df$CD34_Pos_Flow <- FALSE
# meta.df[meta.df$SampleType == "CD34+ Flow",]$CD34_Pos_Flow <- TRUE

meta.df$MSC_Flow <- FALSE
meta.df[meta.df$SampleType == "MSC Flow",]$MSC_Flow <- TRUE

meta.df$Flow <- FALSE
meta.df[grepl("Flow", meta.df$SampleType),]$Flow <- TRUE

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
# test.tmt.global <- tmt$global[ , which(colMeans(!is.na(tmt$global)) >= 0.5)] # no samples removed (all samples had at least half of proteins measured)
# test.tmt.global75 <- tmt$global[ , which(colMeans(!is.na(tmt$global)) >= 0.75)] # no samples removed (all samples had at least 3/4 of proteins measured)
# test.tmt.global80 <- tmt$global[ , which(colMeans(!is.na(tmt$global)) >= 0.80)] # 1 sample removed (requiring at least 4/5 of proteins measured)
# test.tmt.global90 <- tmt$global[ , which(colMeans(!is.na(tmt$global)) >= 0.90)] # 1 sample removed (requiring at least 9/10 of proteins measured)
# tmt80 <- list("meta" = tmt$meta[tmt$meta$id %in% colnames(test.tmt.global80), ],
#               "global" = test.tmt.global80)
# tmt90 <- list("meta" = tmt$meta[tmt$meta$id %in% colnames(test.tmt.global90), ],
#               "global" = test.tmt.global90)

### DIA
meta.df <- readxl::read_excel("Exp24metadataTable_DIA.xlsx") 
meta.df$id <- stringr::str_split_i(meta.df$patient, "-", -1)
meta.df$id <- paste0("X", meta.df$id)
meta.df[meta.df$`sample type` == "CD14+", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD14+", ]$id, "_CD14plus")
meta.df[meta.df$`sample type` == "CD34+", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD34+", ]$id, "_CD34plus")
meta.df[meta.df$`sample type` == "CD14+ Flow", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD14+ Flow", ]$id, "_CD14plusFlow")
meta.df[meta.df$`sample type` == "CD34+ Flow", ]$id <- paste0(meta.df[meta.df$`sample type` == "CD34+ Flow", ]$id, "_CD34plusFlow")
meta.df[meta.df$`sample type` == "MSC Flow", ]$id <- paste0(meta.df[meta.df$`sample type` == "MSC Flow", ]$id, "_MSCflow")

# add other drug info & make sure sensitivity is correctly labeled
meta.df <- merge(meta.df, sens.info)
meta.df$X <- NULL

# add other metadata for contrasts
meta.df$Aza_Sensitive <- FALSE
meta.df[meta.df$Aza == "Sensitive",]$Aza_Sensitive <- TRUE

meta.df$Aza.Ven_Sensitive <- FALSE
meta.df[meta.df$Aza.Ven == "Sensitive",]$Aza.Ven_Sensitive <- TRUE

meta.df$Ven_Sensitive <- FALSE
meta.df[meta.df$Ven == "Sensitive",]$Ven_Sensitive <- TRUE

meta.df$Pooled_CD14_Pos <- FALSE
meta.df[grepl("CD14+", meta.df$`sample type`),]$Pooled_CD14_Pos <- TRUE

meta.df$Pooled_CD34_Pos <- FALSE
meta.df[grepl("CD34+", meta.df$`sample type`),]$Pooled_CD34_Pos <- TRUE

# meta.df$CD14_Pos <- FALSE
# meta.df[meta.df$`sample type` == "CD14+",]$CD14_Pos <- TRUE
# 
# meta.df$CD34_Pos <- FALSE
# meta.df[meta.df$`sample type` == "CD34+",]$CD34_Pos <- TRUE
# 
# meta.df$CD14_Pos_Flow <- FALSE
# meta.df[meta.df$`sample type` == "CD14+ Flow",]$CD14_Pos_Flow <- TRUE
# 
# meta.df$CD34_Pos_Flow <- FALSE
# meta.df[meta.df$`sample type` == "CD34+ Flow",]$CD34_Pos_Flow <- TRUE

meta.df$MSC_Flow <- FALSE
meta.df[meta.df$`sample type` == "MSC Flow",]$MSC_Flow <- TRUE

meta.df$Flow <- FALSE
meta.df[grepl("Flow", meta.df$`sample type`),]$Flow <- TRUE

synapser::synLogin()
#globalFile <- synapser::synGet("syn55233852") # Samantha's processed version
#globalFile <- synapser::synGet("syn55271973") # Camilo's processed version
globalFile <- synapser::synGet("syn55234888") # Samantha's unprocessed version
global.df <- read.table(
  globalFile$path, 
  sep = "\t") # 8897 proteins, 48 samples

# require proteins to be quantified in at least half of samples
global.df <- global.df[which(rowSums(is.na(global.df)) < ncol(global.df)/2),] # 6115 proteins, 48 samples

# require samples to have at least 50% of proteins quantified
#global.df <- global.df[ , which(colSums(is.na(global.df)) < nrow(global.df)/2)] # 42 out of 48 samples are kept
global.df75 <- global.df[ , which(colMeans(!is.na(global.df)) >= 0.75)] # 36 out of 48 samples are kept
# global.df80 <- global.df[ , which(colMeans(!is.na(global.df)) >= 0.80)] # 34 out of 48 samples are kept
# global.df90 <- global.df[ , which(colMeans(!is.na(global.df)) >= 0.90)] # 28 out of 48 samples are kept
# global.df100 <- global.df[ , which(colMeans(!is.na(global.df)) >= 1)] # 0 out of 48 samples are kept
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

# check PCA
library(MSnSet.utils)
library(ggplot2)
rownames(meta.df) <- meta.df$id
meta.df$`Pooled Sample Type` <- meta.df$`sample type`
meta.df[meta.df$`sample type` == "CD14+ Flow",]$`Pooled Sample Type` <- "CD14+"
meta.df[meta.df$`sample type` == "CD34+ Flow",]$`Pooled Sample Type` <- "CD34+"
meta.df[meta.df$`sample type` == "MSC Flow",]$`Pooled Sample Type` <- "MSC"
m_global <- MSnSet(exprs = global.df %>% as.matrix(), 
                   pData = meta.df[colnames(global.df), ])
phenos <- c("id", "patient", "sample type", "Pooled Sample Type", "Aza", "Ven", "Aza.Ven")
for (i in 1:length(phenos)) {
  plot_pca(m_global, phenotype = phenos[i]) + ggtitle("Global PCA") # uses 1917 complete rows (proteins) out of 6115
  ggsave(paste0("exp24_DIA_global_PCA_by_", phenos[i], ".pdf"))
}
m_global <- MSnSet(exprs = global.df75 %>% as.matrix(), 
                   pData = meta.df[colnames(global.df75), ])
phenos <- c("id", "patient", "sample type", "Pooled Sample Type", "Aza", "Ven", "Aza.Ven")
for (i in 1:length(phenos)) {
  plot_pca(m_global, phenotype = phenos[i], label = "patient") + ggtitle("Global PCA") # uses 1917 complete rows (proteins) out of 6115
  ggsave(paste0("exp24_DIA_75percentCoverage_notSampleCentered_global_PCA_by_", phenos[i], ".pdf")) #3018 complete rows
}

# add column for feature names and later make it the first column
#global.df$Gene <- rownames(global.df) # 6092 gene symbols
global.df75$Gene <- rownames(global.df75) # 6092 gene symbols
# global.df80$Gene <- rownames(global.df80) # 6092 gene symbols
# global.df90$Gene <- rownames(global.df90) # 6092 gene symbols
#write.csv(global.df, "Exp24_DIA_crosstab_global_gene_corrected.csv", row.names = FALSE)
#write.csv(global.df75, "Exp24_DIA_75PercentCoverage_crosstab_global_gene_corrected.csv", row.names = FALSE)
# dia <- list("meta" = meta.df[meta.df$id %in% colnames(global.df), ],
#             "global" = global.df)
dia75 <- list("meta" = meta.df[meta.df$id %in% colnames(global.df75), ],
            "global" = global.df75)
# dia80 <- list("meta" = dia$meta[dia$meta$id %in% colnames(global.df80), ],
#               "global" = global.df80)
# dia90 <- list("meta" = dia$meta[dia$meta$id %in% colnames(global.df90), ],
#               "global" = global.df90)
#synfile <- synapser::File("Exp24_DIA_crosstab_global_gene_corrected.csv", "syn54821995")
#synapser::synStore(synfile)
#dia.wo.outliers <- list("meta" = )
### combine DIA & TMT
#dia$meta$method <- "DIA"
dia75$meta$method <- "DIA"
# dia80$meta$method <- "DIA"
# dia90$meta$method <- "DIA"
# tmt80$meta$method <- "TMT"
# tmt90$meta$method <- "TMT"
tmt$meta$method <- "TMT"

# make sure dia & tmt meta data have same columns
#dia$meta[, colnames(dia$meta)[!(colnames(dia$meta) %in% colnames(tmt$meta))]] <- NULL
tmt$meta[, colnames(tmt$meta)[!(colnames(tmt$meta) %in% colnames(dia75$meta))]] <- NULL
dia75$meta[, colnames(dia75$meta)[!(colnames(dia75$meta) %in% colnames(tmt$meta))]] <- NULL
# tmt80$meta[, colnames(tmt80$meta)[!(colnames(tmt80$meta) %in% colnames(dia80$meta))]] <- NULL
# dia80$meta[, colnames(dia80$meta)[!(colnames(dia80$meta) %in% colnames(tmt80$meta))]] <- NULL
# tmt90$meta[, colnames(tmt90$meta)[!(colnames(tmt90$meta) %in% colnames(dia90$meta))]] <- NULL
# dia90$meta[, colnames(dia90$meta)[!(colnames(dia90$meta) %in% colnames(tmt90$meta))]] <- NULL

# dia.tmt <- list("meta" = rbind(dia$meta, tmt$meta),
#                 "global" = merge(dia$global, tmt$global, all = TRUE, 
#                                  suffixes = c("_DIA", "_TMT")),
#                 "phospho" = tmt$phospho)
dia.tmt75 <- list("meta" = rbind(dia75$meta, tmt$meta),
                "global" = merge(dia75$global, tmt$global, all = TRUE, 
                                 suffixes = c("_DIA", "_TMT")))
# dia.tmt80 <- list("meta" = rbind(dia80$meta, tmt80$meta),
#                   "global" = merge(dia80$global, tmt80$global, all = TRUE, 
#                                    suffixes = c("_DIA", "_TMT")))
# dia.tmt90 <- list("meta" = rbind(dia90$meta, tmt90$meta),
#                   "global" = merge(dia90$global, tmt90$global, all = TRUE, 
#                                    suffixes = c("_DIA", "_TMT")))
# don't need to change IDs because they were distinct for DIA & TMT
# # add _DIA and _TMT to end of ids
# dia.tmt$method$id <- paste0(dia.tmt$method$id, "_", dia.tmt$method$method)
# old.colnames <- colnames(dia.tmt$phospho)[1:(ncol(dia.tmt$phospho)-1)]
# new.colnames <- paste0(old.colnames, "_TMT")
# colnames(dia.tmt$phospho)[1:(ncol(dia.tmt$phospho)-1)] <- new.colnames

### PCAs
library(MSnSet.utils)
library(ggplot2)


#### 2. Import BeatAML data formatted for DMEA ####
# import drug MOA annotations
moa.BeatAML <- utils::read.csv(
  "~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
  stringsAsFactors = FALSE, fileEncoding = "latin1")

# load data from Synapse
BeatAML.data <- load_BeatAML_for_DMEA("BeatAML_DMEA_inputs")

#### 3. Run panSEA for each contrast & combination ####
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
setwd(base.path)
# prepare set annotations
if (file.exists("gmt_BeatAML_drug_MOA.rds")) {
  gmt.drug <- readRDS("gmt_BeatAML_drug_MOA.rds")
} else {
  gmt.drug <- DMEA::as_gmt(moa.BeatAML, sep = ", ")
  saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA.rds")
}

setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create("analysis")
setwd("analysis")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"

synapse_id <- "syn53606820"
## run contrasts without filters
contrasts <- colnames(dia.tmt75$meta)[10:(ncol(dia.tmt75$meta)-1)]
priority.contrasts <- c("Flow", "Pooled_CD14_Pos", "Pooled_CD34_Pos", "MSC_Flow")
contrasts <- c(priority.contrasts, contrasts[!(contrasts %in% priority.contrasts)])
contrasts <- c("Sort Type", contrasts[contrasts != "Flow"])
method.data <- list("DIA" = dia75,
                "TMT" = tmt,
                "DIA_&_TMT" = dia.tmt75)
method.data <- list("DIA" = dia75,
                    "DIA_&_TMT" = dia.tmt75,
                    "TMT" = tmt)
# method.data <- list("DIA_75_Percent_Coverage" = dia75,
#                     "TMT_75_Percent_Coverage" = tmt,
#                     "DIA_&_TMT_75_Percent_Coverage" = dia.tmt75)
# method.data <- list("DIA_75_Percent_Coverage" = dia75,
#                     "TMT_75_Percent_Coverage" = tmt,
#                     "DIA_&_TMT_75_Percent_Coverage" = dia.tmt75)
# method.data <- list("DIA_80_Percent_Coverage" = dia80,
#                     "TMT_80_Percent_Coverage" = tmt80,
#                     "DIA_&_TMT_80_Percent_Coverage" = dia.tmt80)
all.degs <- data.frame()

# look at histograms and PCA first
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
  row.names(temp.meta) <- temp.meta$id
  temp.meta$Pooled <- NA
  temp.meta[temp.meta$Pooled_CD14_Pos,]$Pooled <- "CD14+"
  temp.meta[temp.meta$Pooled_CD34_Pos,]$Pooled <- "CD34+"
  temp.meta[temp.meta$MSC_Flow,]$Pooled <- "MSC"
  temp.meta$SampleType <- NA
  temp.meta[temp.meta$Pooled_CD14_Pos,]$SampleType <- "CD14+"
  if (names(method.data)[k] != "TMT") {
    temp.meta[temp.meta$Pooled_CD14_Pos &
                temp.meta$Flow,]$SampleType <- "CD14+ Flow" 
  }
  temp.meta[temp.meta$Pooled_CD34_Pos,]$SampleType <- "CD34+"
  temp.meta[temp.meta$Pooled_CD34_Pos &
              temp.meta$Flow,]$SampleType <- "CD34+ Flow"
  temp.meta[temp.meta$MSC_Flow,]$SampleType <- "MSC Flow"
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
# for DIA_&_TMT PCAs:
# Error in `featureNames<-`(`*tmp*`, value = featureNames(featureData)) : 
#   'value' length (0) must equal feature number in AssayData (6609)

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
                    "TMT" = tmt.wo.out2a,
                    "DIA_&_TMT" = dia.tmt.wo.out2a)
method.data <- list("DIA" = dia.wo.out2a,
                    "TMT" = tmt.wo.out2a)
# determine cc.df
#cc.df <- meta.df[,c("patient", "SampleType")] # also want to include contrast column

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
  row.names(meta.df) <- meta.df$id
  
  meta.df$'Sample Type' <- NA
  meta.df[meta.df$Pooled_CD14_Pos,]$'Sample Type' <- "CD14+"
  meta.df[meta.df$Pooled_CD34_Pos,]$'Sample Type' <- "CD34+"
  meta.df[meta.df$MSC_Flow,]$'Sample Type' <- "MSC"
  
  meta.df$Patient <- meta.df$patient
  meta.df$patient <- NULL
  
  meta.df$'Sort Type' <- "Bead"
  meta.df[meta.df$Flow,]$'Sort Type' <- "Flow"
  
  meta.df$CD14 <- meta.df$Pooled_CD14_Pos
  meta.df[meta.df$Pooled_CD14_Pos,]$CD14 <- "Pos"
  meta.df[!meta.df$Pooled_CD14_Pos,]$CD14 <- "Neg"
  meta.df$Pooled_CD14_Pos <- NULL
  
  meta.df$CD34 <- meta.df$Pooled_CD34_Pos
  meta.df[meta.df$Pooled_CD34_Pos,]$CD34 <- "Pos"
  meta.df[!meta.df$Pooled_CD34_Pos,]$CD34 <- "Neg"
  meta.df$Pooled_CD34_Pos <- NULL
  
  meta.df$MSC <- meta.df$MSC_Flow
  meta.df[meta.df$MSC_Flow,]$MSC <- "MSC"
  meta.df[!meta.df$MSC_Flow,]$MSC <- "Non_MSC"
  meta.df$MSC_Flow <- NULL
  
  # run TF contrast combos
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

# 4th contrast no filter, probably diffexp heatmap:
# Error in hclust(d, method = method) : 
#   NA/NaN/Inf in foreign function call (arg 10)

# # why aren't the Aza_Sensitive_FALSE in dia$global columns?
# aza.res.dia <- c("X01184_CD34plus", "X01184_CD14plus", "X01184_CD14plusFlow", 
#                  "X01184_MSCflow", "X01060_CD34plus", "X01060_CD14plus", 
#                  "X01060_CD34plusFlow", "X01060_CD14plusFlow", "X01060_MSCflow", 
#                  "X00105_CD34plus", "X00105_CD14plus", "X00105_CD34plusFlow", 
#                  "X00105_CD14plusFlow", "X00105_MSCflow", "X00432_CD34plus", 
#                  "X00432_CD14plus", "X00432_CD34plusFlow", 
#                  "X00432_CD14plusFlow", "X00432_MSCflow", "X00839_CD34plus", 
#                  "X00839_CD14plus", "X00839_CD34plusFlow", 
#                  "X00839_CD14plusFlow", "X00839_MSCflow", "X00117_CD34plus",
#                  "X00117_CD14plus", "X00117_CD34plusFlow", 
#                  "X00117_CD14plusFlow", "X00117_MSCflow", "X00251_CD34plus", 
#                  "X00251_CD14plus", "X00251_CD34plusFlow", 
#                  "X00251_CD14plusFlow", "X00251_MSCflow", "X00571_CD34plus",
#                  "X00571_CD14plus", "X00571_CD34plusFlow",
#                  "X00571_CD14plusFlow", "X00571_MSCflow")
# missing.aza.res.dia <- aza.res.dia[!(aza.res.dia %in% colnames(dia$global))] # X00251_MSCflow
# missing.dia.meta <- dia$meta[!(dia$meta$id %in% colnames(dia$global)),] # just X00251_MSCflow

# for TMT run_TF_contrast_combos_global_phospho_human:
# Running ssGSEA using phospho_ksdb data
# Running enrichment analysis...
# Error in `$<-.data.frame`(`*tmp*`, "N_drugs", value = NA) : 
#   replacement has 1 row, data has 0

setwd(base.path)
all.DEG.files <- list("Differential_expression_results.csv" = 
                        all.degs,
                      "Differential_expression_results_max_5_percent_FDR.csv" = 
                        all.degs[all.degs$adj.P.Val <= 0.05, ])
save_to_synapse(all.DEG.files, synapse_id)

##### 4. Plot cell markers and DEGS #####
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
dia.tmt.markers <- dia.tmt75$global[dia.tmt75$global$Gene %in% markers,
                                  c("Gene", colnames(dia.tmt75$global)[colnames(dia.tmt75$global) != "Gene"])]
dia.tmt.markers <- dia.tmt80$global[dia.tmt80$global$Gene %in% markers,
                                    c("Gene", colnames(dia.tmt80$global)[colnames(dia.tmt80$global) != "Gene"])]
write.csv(dia.tmt.markers, "Cell_markers_global_DIA_TMT_75percentCoverage_2024-05-04.csv", row.names = FALSE)
dia.tmt.markers <- dia.tmt90$global[dia.tmt90$global$Gene %in% markers,
                                    c("Gene", colnames(dia.tmt90$global)[colnames(dia.tmt90$global) != "Gene"])]
write.csv(dia.tmt.markers, "Cell_markers_global_DIA_TMT_90percentCoverage_2024-04-17.csv", row.names = FALSE)
dia.tmt.markers <- dia.tmt$global[dia.tmt$global$Gene %in% markers,
                                    c("Gene", colnames(dia.tmt$global)[colnames(dia.tmt$global) != "Gene"])]
write.csv(dia.tmt.markers, paste0("Cell_markers_global_DIA_TMT_", Sys.Date(), ".csv"), row.names = FALSE)

### violin plots of CD14, CD34
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
long.global <- reshape2::melt(dia.tmt75$global, variable.name = "id")
long.global <- merge(long.global, dia.tmt75$meta, by = "id")
long.global <- reshape2::melt(dia.tmt80$global, variable.name = "id")
long.global <- merge(long.global, dia.tmt80$meta, by = "id")
long.global <- reshape2::melt(dia.tmt90$global, variable.name = "id")
long.global <- merge(long.global, dia.tmt90$meta, by = "id")

long.global <- reshape2::melt(dia.tmt$global, variable.name = "id")
long.global <- merge(long.global, dia.tmt$meta, by = "id")

long.global <- reshape2::melt(dia.tmt.wo.out2a$global, variable.name = "id")
long.global <- merge(long.global, dia.tmt.wo.out2a$meta, by = "id")

long.global$Pooled <- NA
long.global[long.global$Pooled_CD14_Pos,]$Pooled <- "CD14+"
long.global[long.global$Pooled_CD34_Pos,]$Pooled <- "CD34+"
long.global[long.global$MSC_Flow,]$Pooled <- "MSC"
long.global$SampleType <- NA
long.global[long.global$Pooled_CD14_Pos,]$SampleType <- "CD14+"
long.global[long.global$Pooled_CD14_Pos & 
              long.global$Flow,]$SampleType <- "CD14+ Flow"
long.global[long.global$Pooled_CD34_Pos,]$SampleType <- "CD34+"
long.global[long.global$Pooled_CD34_Pos & 
              long.global$Flow,]$SampleType <- "CD34+ Flow"
long.global[long.global$MSC_Flow,]$SampleType <- "MSC Flow"

# for CD14: X00105_CD34plusFlow is unusually high, X00074_MSCflow is pretty high and so is X25
# for CD34: X01184_CD14plusFlow and X00251_CD14plusFlow are unusually high

#long.global.markers <- reshape2::melt(dia.tmt.markers)
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

# instead just look at individual data points
# color by sample ID to help identify outliers
dir.create("markers")
setwd("markers")
for (i in 1:length(markers)) {
  marker.df <- na.omit(long.global[long.global$Gene == markers[i],])
  if (nrow(marker.df) > 0) {
    marker.violin <- ggplot2::ggplot(marker.df, 
                                     aes(fill = method, x=Pooled, y=value)) + 
      geom_violin(position=position_dodge(width=0.4), alpha=0.5) + 
      geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0.5) + 
      geom_dotplot(binwidth=0.1, binaxis = "y", position=position_dodge(width=0.4), alpha=0.5, dotsize = 0.5, fill = marker.df$patient) +  
      bg.theme3 + xlab("Sample Type") + ylab("Normalized Protein Expression")
    ggsave(paste0(markers[i],"_75PercentCoverage_by_pooled_sample_type_", Sys.Date(), ".pdf"), marker.violin)
    marker.violin <- ggplot2::ggplot(marker.df, 
                                     aes(fill = method, x=SampleType, y=value)) + 
      geom_violin(position=position_dodge(width=0.4), alpha=0.5) + 
      geom_boxplot(width=0.1, position = position_dodge(width=0.4), alpha=0.5) + 
      geom_dotplot(binwidth=0.1, binaxis = "y", position=position_dodge(width=0.4), alpha=0.5, dotsize = 0.5, fill = marker.df$patient) + 
      bg.theme3 + xlab("Sample Type") + ylab("Normalized Protein Expression")
    ggsave(paste0(markers[i],"_75PercentCoverage_by_sample_type_", Sys.Date(), ".pdf"), marker.violin)
  }
}

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

setwd(file.path(base.path, "DIA", "no_filter"))
coi <- c("CD14_Pos_True_vs_False", "CD34_Pos_True_vs_False")
# source: https://www.novusbio.com/research-areas/stem-cells/mesenchymal-stem-cell-markers#:~:text=Sets%20of%20cell%20surface%20markers,%2C%20CD79a%20and%20HLA%2DDR.
# monocytes should be high in CD14; source: https://www.abcam.com/primary-antibodies/immune-cell-markers-poster
# erythrocytes should be high in CD235a; source: https://www.abcam.com/primary-antibodies/immune-cell-markers-poster
# progenitors should be high in CD34; source: https://www.abcam.com/primary-antibodies/immune-cell-markers-poster
# HSCs (hematopoietic stem cells): CD34+, CD38-, CD45RA-, CD49+, CD90/Thy1+; source: https://www.abcam.com/primary-antibodies/immune-cell-markers-poster
# more reading: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4616896/#:~:text=Hematopoietic%20stem%20cells%20(HSC)%20and,adipocytes%2C%20endothelial%20cells%20and%20myocyte.

sig.degs <- read.csv("Differential_expression_results_max_5_percent_FDR.csv")
filters <- unique(sig.degs$filter)
omics <- list("global" = dia.tmt$global,
              "phospho" = dia.tmt$phospho)
n.degs <- 50
for (k in 1:length(filters)) {
  setwd(base.path)
  setwd("DIA_&_TMT")
  if (filters[k] == "_NA") {
    filter.name <- "no_filter"
  } else {
    filter.name <- filters[k]
  }
  setwd(filter.name)
  for (j in 1:length(contrasts)) {
    sig.degs <- sig.degs.no.filter[sig.degs.no.filter$Contrast == contrasts[j], ]
    if (nrow(sig.degs) > 0) {
      contrast.name <- paste0(contrasts[j], "_TRUE_vs_FALSE")
      if (file.exists(contrast.name)) {
        setwd(contrast.name)
        # get top 25 & bottom 25 degs
        top.sig.degs <- sig.degs %>% slice_max(Log2FC, n.degs/2)
        bot.sig.degs <- sig.degs %>% slice_min(Log2FC, n.degs/2)
        
        for (i in 1:length(omics)) {
          if (i == 1) {
            feature.name <- "Gene"
          } else {
            feature.name <- "SUB_SITE"
          }
          # get data for top & bottom degs
          top.bot.df <- omics[[i]][omics[[i]][,feature.name] %in% top.sig.degs[,feature.name],
                                   c(feature.name, colnames(omics[[i]])[colnames(omics[[i]]) != feature.name])] 
          write.csv(top.bot.df, paste0(names(omics)[i], "_", contrasts[j], 
          "_top_", n.degs/2, "_bottom_", n.degs/2, "_diffexp_max5percentFDR.csv"),
                    row.names = NULL)
        }
      } 
    } 
  }
}


