# Differential expression & enrichment analyses: global & phospho
# PTRC2: exp 24 TMT
# Author: Belinda B. Garana
# Created: 2023-12-06
# Last edit: 2024-02-09

library(readxl); library(panSEA); library(synapser)
library(stringr); library(tidyr); library(dplyr)

setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
source("panSEA_helper.R")

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

meta.df$CD14_Pos <- FALSE
meta.df[meta.df$SampleType == "CD14+",]$CD14_Pos <- TRUE

meta.df$CD34_Pos <- FALSE
meta.df[meta.df$SampleType == "CD34+",]$CD34_Pos <- TRUE

meta.df$CD14_Pos_Flow <- FALSE
meta.df[meta.df$SampleType == "CD14+ Flow",]$CD14_Pos_Flow <- TRUE

meta.df$CD34_Pos_Flow <- FALSE
meta.df[meta.df$SampleType == "CD34+ Flow",]$CD34_Pos_Flow <- TRUE

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
            "global" = global.df,
            "phospho" = phospho.df)

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

meta.df$CD14_Pos <- FALSE
meta.df[meta.df$`sample type` == "CD14+",]$CD14_Pos <- TRUE

meta.df$CD34_Pos <- FALSE
meta.df[meta.df$`sample type` == "CD34+",]$CD34_Pos <- TRUE

meta.df$CD14_Pos_Flow <- FALSE
meta.df[meta.df$`sample type` == "CD14+ Flow",]$CD14_Pos_Flow <- TRUE

meta.df$CD34_Pos_Flow <- FALSE
meta.df[meta.df$`sample type` == "CD34+ Flow",]$CD34_Pos_Flow <- TRUE

meta.df$MSC_Flow <- FALSE
meta.df[meta.df$`sample type` == "MSC Flow",]$MSC_Flow <- TRUE

meta.df$Flow <- FALSE
meta.df[grepl("Flow", meta.df$`sample type`),]$Flow <- TRUE

globalFile <- synapser::synGet("syn54926061")
global.df <- read.table(
  globalFile$path, 
  sep = "\t")

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)

dia <- list("meta" = meta.df,
            "global" = global.df)

### combine DIA & TMT
dia$meta$method <- "DIA"
tmt$meta$method <- "TMT"

# make sure dia & tmt meta data have same columns
tmt$meta[, colnames(tmt$meta)[!(colnames(tmt$meta) %in% colnames(dia$meta))]] <- NULL
dia$meta[, colnames(dia$meta)[!(colnames(dia$meta) %in% colnames(tmt$meta))]] <- NULL

dia.tmt <- list("meta" = rbind(dia$meta, tmt$meta),
                "global" = merge(dia$global, tmt$global, all = TRUE, 
                                 suffixes = c("_DIA", "_TMT")),
                "phospho" = tmt$phospho)

# don't need to change IDs because they were distinct for DIA & TMT
# # add _DIA and _TMT to end of ids
# dia.tmt$method$id <- paste0(dia.tmt$method$id, "_", dia.tmt$method$method)
# old.colnames <- colnames(dia.tmt$phospho)[1:(ncol(dia.tmt$phospho)-1)]
# new.colnames <- paste0(old.colnames, "_TMT")
# colnames(dia.tmt$phospho)[1:(ncol(dia.tmt$phospho)-1)] <- new.colnames

#### 2. Import BeatAML data formatted for DMEA ####
# import drug MOA annotations
moa.BeatAML <- utils::read.csv(
  "~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
  stringsAsFactors = FALSE, fileEncoding = "latin1")

# login to Synapse
synapser::synLogin()

# load data from Synapse
BeatAML.data <- load_BeatAML_for_DMEA("BeatAML_DMEA_inputs")

#### 3. Run panSEA for each contrast & combination ####
## set up comparisons
sens.types <- c("Sensitive", "Resistant")
drug.types <- c("Aza", "Ven", "Aza.Ven")

# prepare set annotations
if (file.exists("gmt_BeatAML_drug_MOA.rds")) {
  gmt.drug <- readRDS("gmt_BeatAML_drug_MOA.rds")
} else {
  gmt.drug <- DMEA::as_gmt(moa.BeatAML, sep = ", ")
  saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA.rds")
}

beatAML <- list("global" = BeatAML.data$global,
                "phospho" = BeatAML.data$phospho)

setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create("analysis")
setwd("analysis")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"

synapse_id <- "syn53653091"
## run contrasts without filters
contrasts <- colnames(meta.df)[4:ncol(meta.df)] # not sure if 4 is right?
method.data <- list("DIA" = dia,
                "TMT" = tmt,
                "DIA_&_TMT" = dia.tmt)
all.degs <- data.frame()
for (k in 1:length(method.data)) {
  setwd(base.path)
  method.path <- file.path(base.path, names(method.data)[k])
  dir.create(names(method.data)[k])
  setwd(names(method.data)[k])
  methodFolder <- 
    synapser::synStore(synapser::Folder(names(method.data)[k],
                                        parent = synapse_id))
  
  meta.df <- method.data[[k]]$meta
  
  # run TF contrast combos
  if (names(method.data)[k] == "TMT") {
    omics <- list("global" = method.data[[k]]$global)
    run_TF_contrast_combos_global_human(contrasts, "id", meta.df, omics,
                                                base.path = method.path, 
                                                synapse_id = methodFolder)
  } else {
    omics <- list("global" = method.data[[k]]$global,
                  "phospho" = method.data[[k]]$phospho)
    run_TF_contrast_combos_global_phospho_human(contrasts, "id", meta.df, omics,
                                                base.path = method.path, 
                                                synapse_id = methodFolder) 
  }
  
  # get compiled DEGs
  methodDEGs <- synapser::synGetChildren(methodFolder, list("file"), sortBy = 'NAME')
  methodFile <- synapser::synGet(methodDEGs[[1]]$id)
  methodDEG <- read.csv(methodFile$path)
  methodDEG$method <- names(method.data)[k]
  all.degs <- rbind(all.degs, methodDEG)
}
setwd(base.path)
all.DEG.files <- list("Differential_expression_results.csv" = 
                        all.degs,
                      "Differential_expression_results_max_5_percent_FDR.csv" = 
                        all.degs[all.degs$adj.P.Val <= 0.05, ])
save_to_synapse(all.DEG.files, synapse_id)