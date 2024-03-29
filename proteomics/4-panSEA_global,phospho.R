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

synapse_id <- "syn53606820"
## run contrasts without filters
contrasts <- colnames(dia.tmt$meta)[10:(ncol(dia.tmt$meta)-1)]
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
  if (grepl("DIA", names(method.data)[k])) {
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

##### 4. Plot cell markers and DEGS #####
# plot markers
# MSCs (Mesenchymal Stem Cells): should be high in CD73, CD90, CD105, CD106, CD146, STRO-1 and low in CD14, CD34, CD45, HLA-DR
# AML monocytes should be high in CD4, CD11c, CD14, and CD64; source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5907644/#:~:text=Surface%20antigens%20that%20indicate%20leukemic,mostly%20expressed%20on%20mature%20monocytes.
# Leukemia stem cells (LSCs):  CD34, CD38, CD123, TIM3, CD25, CD32 and CD96; source: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5388677/
# LSCs are emphasized in recent PTRC paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10950354/

markers <- unique(c("CD73", "CD90", "CD105", "CD106", "CD146", "STRO-1", "CD14", 
                    "CD34", "CD45", "HLA-DR", "CD4", "CD11c", "CD14", "CD64", 
                    "CD34", "CD38", "CD123", "TIM3", "CD25", "CD32" "CD96"))
dia.tmt.markers <- dia.tmt$global[dia.tmt$global$Gene %in% markers,
                                  c("Gene", colnames(dia.tmt$global)[colnames(dia.tmt$global) != "Gene"])]
write.csv(dia.tmt.markers, "Cell_markers_global_DIA_TMT.csv", row.names = NULL)

# source: https://www.novusbio.com/research-areas/stem-cells/mesenchymal-stem-cell-markers#:~:text=Sets%20of%20cell%20surface%20markers,%2C%20CD79a%20and%20HLA%2DDR.
# monocytes should be high in CD14; source: https://www.abcam.com/primary-antibodies/immune-cell-markers-poster
# erythrocytes should be high in CD235a; source: https://www.abcam.com/primary-antibodies/immune-cell-markers-poster
# progenitors should be high in CD34; source: https://www.abcam.com/primary-antibodies/immune-cell-markers-poster
# HSCs (hematopoietic stem cells): CD34+, CD38-, CD45RA-, CD49+, CD90/Thy1+; source: https://www.abcam.com/primary-antibodies/immune-cell-markers-poster
# more reading: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4616896/#:~:text=Hematopoietic%20stem%20cells%20(HSC)%20and,adipocytes%2C%20endothelial%20cells%20and%20myocyte.

sig.degs.no.filter <- all.degs[is.na(all.degs$filter) & 
                                 all.degs$adj.P.Val <= 0.05, ]
for (j in 1:length(contrasts)) {
  sig.degs <- sig.degs.no.filter[sig.degs.no.filter$Contrast == contrasts[j], ]
  
  # get top 25 & bottom 25 degs
  top.sig.degs <- sig.degs %>% slice_max(Log2FC, 25)
  bot.sig.degs <- sig.degs %>% slice_min(Log2FC, 25)
  
  for (i in 1:length(omics)) {
    if (i == 1) {
      feature.name <- "Gene"
    } else {
      feature.name <- "SUB_SITE"
    }
    # get data for top & bottom degs
    top.bot.df <- omics[[i]][omics[[i]][,feature.name] %in% top.sig.degs[,feature.name],
                             c(feature.name, colnames(omics[[i]])[colnames(omics[[i]]) != feature.name])] 
    write.csv(top.bot.df, paste0(names(omics)[i], "_", contrasts[j], "_top_25_bottom_25_diffexp_max5percentFDR.csv"),
              row.names = NULL)
  }
}