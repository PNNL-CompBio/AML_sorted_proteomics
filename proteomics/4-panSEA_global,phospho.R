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
setwd("data")
meta.df <- readxl::read_excel("Exp24metadataTable_TMT.xlsx") 
meta.df$MeasurementName <- as.character(rownames(meta.df))
rownames(meta.df) <- paste0("X", rownames(meta.df))
meta.df$row.name <- rownames(meta.df)

# add other drug info & make sure sensitivity is correctly labeled
sens.info <- read.csv("Exp24_drug_sensitivity_20240209.csv")
meta.df <- merge(meta.df, sens.info)
meta.df$X <- NULL

global.df <- read.table(
  "global_data/Exp24_crosstab_global_gene_corrected.txt", 
  sep = "\t")
phospho.df <- read.table(
  "phospho_data/Exp24_crosstab_phospho_SiteID_corrected.txt", 
  sep = "\t")

# add column for feature names and later make it the first column
global.df$Gene <- rownames(global.df)
phospho.df$SUB_SITE <- rownames(phospho.df)

#### 2. Import BeatAML data formatted for DMEA ####
# import drug MOA annotations
moa.BeatAML <- utils::read.csv(
  "~/OneDrive - PNNL/Documents/PTRC2/BeatAML_single_drug_moa.csv",
  stringsAsFactors = FALSE, fileEncoding = "latin1")

# login to Synapse
#synapser::synLogin()

# load data from Synapse
BeatAML.data <- load_BeatAML_for_DMEA("BeatAML_DMEA_inputs")

#### 3. Run panSEA across sort.types for each omics, sens.type, drug.type ####
## set up comparisons
sort.types <- unique(na.omit(meta.df$SampleType))
sens.types <- unique(na.omit(meta.df$Drug))
drug.types <- c("Aza", "Ven", "Aza.Ven")

# sort.type contrasts without spaces or + in names
meta.df$CellType <- meta.df$SampleType
meta.df[meta.df$SampleType == "CD14+",]$CellType <- "CD14_Pos"
meta.df[meta.df$SampleType == "CD34+",]$CellType <- "CD34_Pos"
meta.df[meta.df$SampleType == "CD14+ Flow",]$CellType <- "CD14_Pos_Flow"
meta.df[meta.df$SampleType == "CD34+ Flow",]$CellType <- "CD34_Pos_Flow"
meta.df[meta.df$SampleType == "MSC Flow",]$CellType <- "MSC_Flow"
sort.types4 <- na.omit(unique(meta.df$CellType))
sort.contrasts4 <- list()
counter <- 0
for (i in 1:length(sort.types4)) {
  for (j in 1:length(sort.types4)) {
    if (i != j) {
      counter <- counter + 1
      alpha.pair <- sort(c(sort.types4[i], sort.types4[j]))
      sort.contrasts4[[counter]] <- alpha.pair
    }
  }
}
sort.contrasts4 <- unique(sort.contrasts4)

# drug.sens contrasts
sens.contrasts <- list()
counter <- 0
for (i in drug.types) {
  sens.contrasts[[i]] <- c("Sensitive", "Resistant")
}

# get set annotations for pathway analyses
## prepare set annotations
# generate gmt.features beforehand to save time
gmt <- readRDS("gmt_MSigDB_Homo-sapiens_C2_CP_KEGG.rds")

msigdb.info <- msigdbr::msigdbr("Homo sapiens", "H")

# extract necessary info into data frame
msigdb.info <- as.data.frame(msigdb.info[, c(
  "gene_symbol",
  "gs_name",
  "gs_description"
)])

# gmt.H <- DMEA::as_gmt(
#   msigdb.info, "gene_symbol", "gs_name", min.per.set = 6,
#   descriptions = "gs_description"
# ) 
# saveRDS(gmt.H, "gmt_MSigDB_Homo-sapiens_H.rds")
gmt.H <- readRDS("gmt_MSigDB_Homo-sapiens_H.rds")

gmt.ksdb <- readRDS("gmt_ksdb_human_PNNL.rds")

# only create substrate gmt the first time
# SUB_SITE <- phospho.df$SUB_SITE
# phospho.ref <- data.frame(SUB_SITE)
# phospho.ref <- phospho.ref %>% tidyr::extract(SUB_SITE, "KINASE",
#                                               remove = FALSE)
# SUB_SITE <- NULL
# gmt.sub <- DMEA::as_gmt(phospho.ref, "SUB_SITE", "KINASE", min.per.set = 6)
# saveRDS(gmt.sub, "gmt_PNNL_kinase-substrate_PTRC2_exp24.rds")
gmt.sub <- readRDS("gmt_PNNL_kinase-substrate_PTRC2_exp24.rds")

#gmt.drug <- DMEA::as_gmt(moa.BeatAML, sep = ", ")
#saveRDS(gmt.drug, "gmt_BeatAML_drug_MOA_2024-02-22.rds")
gmt.drug <- readRDS("gmt_BeatAML_drug_MOA_2024-02-22.rds")
omics <- list("global" = global.df,
              "phospho" = phospho.df)
beatAML <- list("global" = BeatAML.data$global,
                "phospho" = BeatAML.data$phospho)

# prepare universal inputs
gmt.features <- list(gmt, gmt.ksdb)
gmt.features2 <- list(gmt.H, gmt.sub)
gmt.features3a <- list(gmt, gmt) # for phospho GSEA^2
gmt.features3b <- list(gmt.H, gmt.H) # for phospho GSEA^2

setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
dir.create("analysis")
setwd("analysis")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
# add 'X' to MeasurementName to match colnames of data
meta.df$id <- paste0('X', meta.df$MeasurementName)

synapse_id <- "syn53653091"

# run sens vs. res contrasts without cell type filter
for (j in drug.types) {
  setwd(base.path)
  run_contrasts_global_phospho(list(c("Sensitive", "Resistant")), j, 'id',
                               meta.df, omics, beatAML, gmt.features, 
                               gmt.features2, gmt.features3a, gmt.features3b, 
                               gmt.drug, BeatAML.data$drug, base.path,
                               synapse_id)
}
for (i in sort.types) {
  # filter meta.df
  meta.filtered <- meta.df[meta.df$SampleType == i, ]
  
  # create folder for each cell type
  setwd(base.path)
  dir.create(file.path(i))
  setwd(file.path(i))
  temp.base.path <- file.path(base.path, i)
  dataFolder <- 
    synapser::synStore(synapser::Folder(file.path(i),
                                        parent = synapse_id))
  
  # run sens vs. res contrasts for each cell type
  for (j in drug.types) {
    if (nrow(meta.filtered[meta.filtered[ , j] == "Sensitive", ]) > 0 &
        nrow(meta.filtered[meta.filtered[ , j] == "Resistant", ]) > 0) {
      run_contrasts_global_phospho(list(c("Sensitive", "Resistant")), j, 'id', 
                                   meta.filtered, omics, beatAML, gmt.features, 
                                   gmt.features2, gmt.features3a, gmt.features3b, 
                                   gmt.drug, BeatAML.data$drug, temp.base.path,
                                   synapse_id = dataFolder) 
    }
  } 
}

# run cell type contrasts without drug sensitivity filter
setwd(base.path)
contrasts.w.KSEA <- sort.contrasts4[c(1:2, 4:5, 7, 9)]
contrasts.wo.KSEA <- sort.contrasts4[c(3, 6, 8, 10)]
run_contrasts_global_phospho(contrasts.w.KSEA, "CellType",'id', meta.df,
                             omics, beatAML, gmt.features, gmt.features2, 
                             gmt.features3a, gmt.features3b, 
                             gmt.drug, BeatAML.data$drug, base.path,
                             synapse_id)

run_contrasts_global_phospho(contrasts.wo.KSEA, "CellType",'id', meta.df,
                             omics, beatAML, gmt.features, gmt.features2, 
                             gmt.features3a, gmt.features3b, 
                             gmt.drug, BeatAML.data$drug, base.path,
                             synapse_id, KSEA = FALSE)

for (j in drug.types) {
  #j = "Ven"
  #j = "Aza.Ven"
  for (i in sens.types) {
    #i = "Resistant"
    # filter meta.df
    meta.filtered <- meta.df[meta.df[ , j] == i, ]
    
    # create folder for each drug sensitivity type
    drug.sens.type <- paste0(j, "_", i)
    setwd(base.path)
    dir.create(drug.sens.type)
    setwd(drug.sens.type)
    temp.base.path <- file.path(base.path, drug.sens.type)
    dataFolder <- 
      synapser::synStore(synapser::Folder(drug.sens.type,
                                          parent = synapse_id))
    
    # run cell type contrasts for each drug sensitivity
    if (i == "Resistant") {
      contrasts.w.KSEA <- sort.contrasts4[c(1:2, 4:5, 7, 9)]
      if (j == "Aza") {
        contrasts.wo.KSEA <- sort.contrasts4[c(3, 6, 8, 10)]
      } else {
        contrasts.wo.KSEA <- sort.contrasts4[c(3, 6, 8)] 
      }
      run_contrasts_global_phospho(contrasts.w.KSEA, "CellType", 'id',
                                   meta.filtered, omics, beatAML, gmt.features,
                                   gmt.features2, gmt.features3a, gmt.features3b,
                                   gmt.drug, BeatAML.data$drug, temp.base.path,
                                   synapse_id = dataFolder)
      run_contrasts_global_phospho(contrasts.wo.KSEA, "CellType", 'id',
                                   meta.filtered, omics, beatAML, gmt.features,
                                   gmt.features2, gmt.features3a, gmt.features3b,
                                   gmt.drug, BeatAML.data$drug, temp.base.path,
                                   synapse_id = dataFolder, KSEA = FALSE)
    } else {
      if (i == "Sensitive" & j == "Aza") {
        # not enough samples for Aza_Sensitive items 2, 3, 4 of sort.contrasts4:
        # Error in .ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim,  : 
        #                    No residual degrees of freedom in linear model fits
        run_contrasts_global_phospho(sort.contrasts4[c(1,7)], "CellType", 'id', 
                                     meta.filtered, omics, beatAML, gmt.features, 
                                     gmt.features2, gmt.features3a, gmt.features3b, 
                                     gmt.drug, BeatAML.data$drug, temp.base.path,
                                     synapse_id = dataFolder) 
      } else {
        run_contrasts_global_phospho(sort.contrasts4[1:length(sort.contrasts4)], "CellType", 'id', 
                                     meta.filtered, omics, beatAML, gmt.features, 
                                     gmt.features2, gmt.features3a, gmt.features3b, 
                                     gmt.drug, BeatAML.data$drug, temp.base.path,
                                     synapse_id = dataFolder)  
      }
    }
    
    # no KSEA for Aza_Resistant sort.contrasts4 items 3, 6, 8, 10
    # no KSEA for Aza.Ven_Resistant sort.contrasts4 items 3, 6, 8
    # no KSEA for Ven_Resistant sort.contrasts4 items 3, 6, 8
    # not enough samples for Aza.Ven_Resistant sort.contrasts4 item 10 or Ven_Resistant sort.contrasts4 item 10:
    # Error in .ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim,  : 
    #                    No residual degrees of freedom in linear model fits
    # run_contrasts_global_phospho(sort.contrasts4[8], "CellType", 'id',
    #                              meta.filtered, omics, beatAML, gmt.features,
    #                              gmt.features2, gmt.features3a, gmt.features3b,
    #                              gmt.drug, BeatAML.data$drug, temp.base.path,
    #                              synapse_id = dataFolder, KSEA = FALSE)
  } 
}

# pool flow and non-flow samples
meta.df$Pooled <- meta.df$SampleType
meta.df[grepl("CD14+", meta.df$SampleType),]$Pooled <- "CD14_Pos"
meta.df[grepl("CD34+", meta.df$SampleType),]$Pooled <- "CD34_Pos"
meta.df[meta.df$SampleType == "MSC Flow",]$Pooled <- "MSC_Flow"
sort.types.pooled <- na.omit(unique(meta.df$Pooled))
sort.pooled <- list()
counter <- 0
for (i in 1:length(sort.types.pooled)) {
  for (j in 1:length(sort.types.pooled)) {
    if (i != j) {
      counter <- counter + 1
      alpha.pair <- sort(c(sort.types.pooled[i], sort.types.pooled[j]))
      sort.pooled[[counter]] <- alpha.pair
    }
  }
}
sort.pooled <- unique(sort.pooled)

# filter for pooled cell types
for (i in sort.types.pooled) {
  # filter meta.df
  meta.filtered <- meta.df[meta.df$Pooled == i, ]
  
  # create folder for each cell type
  setwd(base.path)
  folder.name <- paste0("Pooled_", i)
  dir.create(folder.name)
  setwd(folder.name)
  temp.base.path <- file.path(base.path, folder.name)
  dataFolder <- 
    synapser::synStore(synapser::Folder(folder.name,
                                        parent = synapse_id))
  
  # run sens vs. res contrasts for each cell type
  for (j in drug.types) {
    if (nrow(meta.filtered[meta.filtered[ , j] == "Sensitive", ]) > 0 &
        nrow(meta.filtered[meta.filtered[ , j] == "Resistant", ]) > 0) {
      run_contrasts_global_phospho(list(c("Sensitive", "Resistant")), j, 'id', 
                                   meta.filtered, omics, beatAML, gmt.features, 
                                   gmt.features2, gmt.features3a, gmt.features3b, 
                                   gmt.drug, BeatAML.data$drug, temp.base.path,
                                   synapse_id = dataFolder) 
    }
  } 
}

# run pooled cell type contrasts without drug sensitivity filter
setwd(base.path)
contrasts.w.KSEA <- sort.pooled[1:length(sort.pooled)]
run_contrasts_global_phospho(contrasts.w.KSEA, "Pooled",'id', meta.df,
                             omics, beatAML, gmt.features, gmt.features2, 
                             gmt.features3a, gmt.features3b, 
                             gmt.drug, BeatAML.data$drug, base.path,
                             synapse_id)


for (j in drug.types) {
  #j = "Ven"
  #j = "Aza.Ven"
  for (i in sens.types) {
    #i = "Resistant"
    # filter meta.df
    meta.filtered <- meta.df[meta.df[ , j] == i, ]
    
    # create folder for each drug sensitivity type
    drug.sens.type <- paste0(j, "_", i)
    setwd(base.path)
    dir.create(drug.sens.type)
    setwd(drug.sens.type)
    temp.base.path <- file.path(base.path, drug.sens.type)
    dataFolder <- 
      synapser::synStore(synapser::Folder(drug.sens.type,
                                          parent = synapse_id))
    
    # run cell type contrasts for each drug sensitivity
    run_contrasts_global_phospho(sort.pooled, "Pooled", 'id',
                                 meta.filtered, omics, beatAML, gmt.features,
                                 gmt.features2, gmt.features3a, gmt.features3b,
                                 gmt.drug, BeatAML.data$drug, temp.base.path,
                                 synapse_id = dataFolder)

  }
} 

#### 4. Compile differential expression (DEGs) ####
all.degs <- data.frame()
# get degs from sens vs. res contrasts without cell type filter
for (j in drug.types) {
  contrast.type <- paste0(j, "_Sensitive_vs_Resistant")
  setwd(file.path(base.path, contrast.type)) 
  
  # load degs
  global.degs <- read.csv("global/Differential_expression/Differential_expression_results.csv")
  phospho.degs <- read.csv("phospho/Differential_expression/Differential_expression_results.csv")
  
  # add feature_type
  global.degs$Feature_type <- colnames(global.degs)[1]
  phospho.degs$Feature_type <- colnames(phospho.degs)[1]
  colnames(global.degs)[1] <- "Feature"
  colnames(phospho.degs)[1] <- "Feature"
  
  # rbind and add contrast information
  temp.degs <- na.omit(rbind(global.degs, phospho.degs))
  temp.degs$Cell_type_filter <- NA
  temp.degs$Drug_sensitivity_filter <- NA
  temp.degs$Contrast <- contrast.type
  all.degs <- rbind(all.degs, temp.degs)
}

# get degs from sens vs. res contrasts for each cell type
for (i in sort.types) {
  for (j in drug.types) {
    contrast.type <- paste0(j, "_Sensitive_vs_Resistant")
    setwd(base.path)
    setwd(i) # change to sort.type path separately because of potential space
    #dir.create(contrast.type) # in case it wasn't created before
    #setwd(contrast.type)
    
    if (file.exists(file.path(contrast.type, "global/Differential_expression/Differential_expression_results.csv"))) {
      setwd(contrast.type)
      
      # load degs
      global.degs <- read.csv("global/Differential_expression/Differential_expression_results.csv")
      phospho.degs <- read.csv("phospho/Differential_expression/Differential_expression_results.csv")
      
      # add feature_type
      global.degs$Feature_type <- colnames(global.degs)[1]
      phospho.degs$Feature_type <- colnames(phospho.degs)[1]
      colnames(global.degs)[1] <- "Feature"
      colnames(phospho.degs)[1] <- "Feature"
      
      # rbind and add contrast information
      temp.degs <- na.omit(rbind(global.degs, phospho.degs))
      temp.degs$Cell_type_filter <- i
      temp.degs$Drug_sensitivity_filter <- NA
      temp.degs$Contrast <- contrast.type
      all.degs <- rbind(all.degs, temp.degs)
    }
  }
}

# get degs from cell type contrasts without drug sensitivity filter
for (i in 1:length(sort.contrasts4)) {
    contrast.type <- paste0(sort.contrasts4[[i]][1], "_vs_", sort.contrasts4[[i]][2])
    contrast.type.path <- paste0("CellType_", contrast.type)
    setwd(file.path(base.path, contrast.type.path))
    
    if (file.exists("global/Differential_expression/Differential_expression_results.csv")) {
      # load degs
      global.degs <- read.csv("global/Differential_expression/Differential_expression_results.csv")
      phospho.degs <- read.csv("phospho/Differential_expression/Differential_expression_results.csv")
      
      # add feature_type
      global.degs$Feature_type <- colnames(global.degs)[1]
      phospho.degs$Feature_type <- colnames(phospho.degs)[1]
      colnames(global.degs)[1] <- "Feature"
      colnames(phospho.degs)[1] <- "Feature"
      
      # rbind and add contrast information
      temp.degs <- na.omit(rbind(global.degs, phospho.degs))
      temp.degs$Cell_type_filter <- NA
      temp.degs$Drug_sensitivity_filter <- NA
      temp.degs$Contrast <- contrast.type
      all.degs <- rbind(all.degs, temp.degs)
    }
}

# get degs from cell type contrasts for each drug sensitivity
for (j in drug.types) {
  for (k in sens.types) {
    drug.sens.type <- paste0(j, "_", k)
    for (i in 1:length(sort.contrasts4)) {
      contrast.type <- paste0(sort.contrasts4[[i]][1], "_vs_", sort.contrasts4[[i]][2])
      contrast.type.path <- paste0("CellType_", contrast.type)
      setwd(base.path)
      #dir.create(drug.sens.type)
      setwd(drug.sens.type)
      #dir.create(contrast.type.path)
      #setwd(contrast.type.path)
      
      if (file.exists(file.path(contrast.type.path, "global/Differential_expression/Differential_expression_results.csv"))) {
        setwd(contrast.type.path)
        
        # load degs
        global.degs <- read.csv("global/Differential_expression/Differential_expression_results.csv")
        phospho.degs <- read.csv("phospho/Differential_expression/Differential_expression_results.csv")
        
        # add feature_type
        global.degs$Feature_type <- colnames(global.degs)[1]
        phospho.degs$Feature_type <- colnames(phospho.degs)[1]
        colnames(global.degs)[1] <- "Feature"
        colnames(phospho.degs)[1] <- "Feature"
        
        # rbind and add contrast information
        temp.degs <- na.omit(rbind(global.degs, phospho.degs))
        temp.degs$Cell_type_filter <- NA
        temp.degs$Drug_sensitivity_filter <- drug.sens.type
        temp.degs$Contrast <- contrast.type
        all.degs <- rbind(all.degs, temp.degs)
      }
    }
  }
}

# save locally
setwd(base.path)
write.csv(all.degs, "Exp24_TMT_differential_expression.csv", row.names = FALSE) # 379,420 features
write.csv(all.degs[all.degs$adj.P.Val <= 0.05, ], "Exp24_TMT_differential_expression_max_5_percent_FDR.csv", row.names = FALSE) # 584 features
filtered.degs <- all.degs[all.degs$adj.P.Val <= 0.05, ]

CSV.files <- list.files(pattern = ".*.csv", full.names = TRUE)
if (length(CSV.files) > 0) {
  # save to synapse
  CSVs <- lapply(as.list(CSV.files), synapser::File,
                 parent = synapse_id)
  lapply(CSVs, synapser::synStore)
}

# add on pooled analyses
setwd(base.path)
all.degs <- read.csv("Exp24_TMT_differential_expression.csv")

# get degs filtered for pooled cell types
for (i in sort.types.pooled) {
  for (j in drug.types) {
    contrast.type <- paste0(j, "_Sensitive_vs_Resistant")
    setwd(base.path)
    folder.name <- paste0("Pooled_", i)
    setwd(folder.name) # change to sort.type path separately because of potential space
    #dir.create(contrast.type) # in case it wasn't created before
    #setwd(contrast.type)
    
    if (file.exists(file.path(contrast.type, "global/Differential_expression/Differential_expression_results.csv"))) {
      setwd(contrast.type)
      
      # load degs
      global.degs <- read.csv("global/Differential_expression/Differential_expression_results.csv")
      phospho.degs <- read.csv("phospho/Differential_expression/Differential_expression_results.csv")
      
      # add feature_type
      global.degs$Feature_type <- colnames(global.degs)[1]
      phospho.degs$Feature_type <- colnames(phospho.degs)[1]
      colnames(global.degs)[1] <- "Feature"
      colnames(phospho.degs)[1] <- "Feature"
      
      # rbind and add contrast information
      temp.degs <- na.omit(rbind(global.degs, phospho.degs))
      temp.degs$Cell_type_filter <- folder.name
      temp.degs$Drug_sensitivity_filter <- NA
      temp.degs$Contrast <- contrast.type
      all.degs <- rbind(all.degs, temp.degs)
    }
  }
}

# get degs from pooled cell type contrasts without drug sensitivity filter
for (i in 1:length(sort.pooled)) {
  contrast.type <- paste0("Pooled_", sort.pooled[[i]][1], "_vs_", sort.pooled[[i]][2])
  setwd(file.path(base.path, contrast.type))
  
  if (file.exists("global/Differential_expression/Differential_expression_results.csv")) {
    # load degs
    global.degs <- read.csv("global/Differential_expression/Differential_expression_results.csv")
    phospho.degs <- read.csv("phospho/Differential_expression/Differential_expression_results.csv")
    
    # add feature_type
    global.degs$Feature_type <- colnames(global.degs)[1]
    phospho.degs$Feature_type <- colnames(phospho.degs)[1]
    colnames(global.degs)[1] <- "Feature"
    colnames(phospho.degs)[1] <- "Feature"
    
    # rbind and add contrast information
    temp.degs <- na.omit(rbind(global.degs, phospho.degs))
    temp.degs$Cell_type_filter <- NA
    temp.degs$Drug_sensitivity_filter <- NA
    temp.degs$Contrast <- contrast.type
    all.degs <- rbind(all.degs, temp.degs)
  }
}

# get degs from pooled cell type contrasts for each drug sensitivity
for (j in drug.types) {
  for (k in sens.types) {
    drug.sens.type <- paste0(j, "_", k)
    for (i in 1:length(sort.pooled)) {
      contrast.type <- paste0("Pooled_", sort.pooled[[i]][1], "_vs_", sort.pooled[[i]][2])
      setwd(base.path)
      setwd(drug.sens.type)
      
      if (file.exists(file.path(contrast.type, "global/Differential_expression/Differential_expression_results.csv"))) {
        setwd(contrast.type)
        
        # load degs
        global.degs <- read.csv("global/Differential_expression/Differential_expression_results.csv")
        phospho.degs <- read.csv("phospho/Differential_expression/Differential_expression_results.csv")
        
        # add feature_type
        global.degs$Feature_type <- colnames(global.degs)[1]
        phospho.degs$Feature_type <- colnames(phospho.degs)[1]
        colnames(global.degs)[1] <- "Feature"
        colnames(phospho.degs)[1] <- "Feature"
        
        # rbind and add contrast information
        temp.degs <- na.omit(rbind(global.degs, phospho.degs))
        temp.degs$Cell_type_filter <- NA
        temp.degs$Drug_sensitivity_filter <- drug.sens.type
        temp.degs$Contrast <- contrast.type
        all.degs <- rbind(all.degs, temp.degs)
      }
    }
  }
}

# save locally
setwd(base.path)
write.csv(all.degs, "Exp24_TMT_differential_expression.csv", row.names = FALSE) # now 539449 rows; was 379,420 before pooled analyses
write.csv(all.degs[all.degs$adj.P.Val <= 0.05, ], "Exp24_TMT_differential_expression_max_5_percent_FDR.csv", row.names = FALSE) # now 2128 rows; was 584 before pooled analyses
filtered.degs <- all.degs[all.degs$adj.P.Val <= 0.05, ]

CSV.files <- c("Exp24_TMT_differential_expression.csv", "Exp24_TMT_differential_expression_max_5_percent_FDR.csv")
if (length(CSV.files) > 0) {
  # save to synapse
  CSVs <- lapply(as.list(CSV.files), synapser::File,
                 parent = synapse_id)
  lapply(CSVs, synapser::synStore)
}

# create heatmaps for top 30-50 differentially expressed features
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
my.contrast <- "Pooled_CD14_Pos_vs_CD34_Pos"
n.features <- c(10, 20, 30, 50)
omics <- c("global", "phospho")
max.abs.Log2FC <- 0 # 2.94455 so use 3
max.minusLogFDR <- 0 # 4.343 so use 4.5
for (i in omics) {
  setwd(file.path(base.path, my.contrast, i, "Differential_expression"))
  # load differential expression results & filter for q <= 0.05
  degs <- read.csv("Differential_expression_results.csv")
  degs <- degs[degs$adj.P.Val <= 0.05, ]
  for (j in n.features) {
    # identify top features
    top.degs <- degs %>% dplyr::slice_max(abs(Log2FC), n = j)
    top.degs$minusLogFDR <- -log(top.degs$adj.P.Val, 10)
    
    # save top differential expression results
    write.csv(top.degs, paste0("top_", j, "_diffexp_features_maxFDR0.05.csv"), row.names = FALSE)
    if (i == "global") {
      write.csv(top.degs[,c("Gene", "Log2FC")], paste0("top_", j, "_diffexp_features_Log2FC_maxFDR0.05.csv"), row.names = FALSE)
      write.csv(top.degs[,c("Gene", "minusLogFDR")], paste0("top_", j, "_diffexp_features_minusLogFDR_maxFDR0.05.csv"), row.names = FALSE)
    } else {
      write.csv(top.degs[,c("SUB_SITE", "Log2FC")], paste0("top_", j, "_diffexp_features_Log2FC_maxFDR0.05.csv"), row.names = FALSE)
      write.csv(top.degs[,c("SUB_SITE", "minusLogFDR")], paste0("top_", j, "_diffexp_features_minusLogFDR_maxFDR0.05.csv"), row.names = FALSE)
    }
    
    # keep track of max abs(Log2FC) and minusLogFDR
    temp.max.abs.Log2FC <- max(abs(top.degs$Log2FC))
    temp.max.minusLogFDR <- max(top.degs$minusLogFDR, na.rm = TRUE)
    max.abs.Log2FC <- ifelse(temp.max.abs.Log2FC > max.abs.Log2FC, 
                             temp.max.abs.Log2FC, max.abs.Log2FC)
    max.minusLogFDR <- ifelse(temp.max.minusLogFDR > max.minusLogFDR, 
                             temp.max.minusLogFDR, max.minusLogFDR)
  } 
}

# repeat but looking at individual sample expression of top diffexp proteins/phospho-sites
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
my.contrast <- "Pooled_CD14_Pos_vs_CD34_Pos"
n.features <- c(10, 20, 30, 50)
omics <- c("global", "phospho")
max.expr <- 0 
global.df <- global.df[,c("Gene", colnames(global.df)[1:(ncol(global.df)-1)])]
phospho.df <- phospho.df[,c("SUB_SITE", colnames(phospho.df)[1:(ncol(phospho.df)-1)])]
expr <- list("global" = global.df, "phospho" = phospho.df)
for (i in omics) {
  setwd(file.path(base.path, my.contrast, i, "Differential_expression"))
  temp.expr <- expr[[i]]
  if (i == "global") {
    feature.name <- "Gene"
  } else {
    feature.name <- "SUB_SITE"
  }

  for (j in n.features) {
    # identify top features
    top.degs <- read.csv(paste0("top_", j, "_diffexp_features_maxFDR0.05.csv"))
    
    # get expression data for top features
    top.expr <- temp.expr[rownames(temp.expr) %in% top.degs[ , feature.name], ]
    
    # save expression data for top features
    write.csv(top.expr, 
              paste0("top_", j, "_", feature.name, "_maxFDR0.05.csv"), 
              row.names = FALSE)
    
    # keep track of max abs(expression)
    temp.max.expr <- max(abs(top.expr[ , colnames(top.expr) != feature.name]))
    max.expr <- ifelse(temp.max.expr > max.expr, temp.max.expr, max.expr)
  } 
}

## repeat but balance up & down by pulling half from top & half from bottom based on Log2FC
# create heatmaps for top 30-50 differentially expressed features
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
my.contrast <- "Pooled_CD14_Pos_vs_CD34_Pos"
n.features <- c(10, 20, 30, 50)
omics <- c("global", "phospho")
max.abs.Log2FC <- 0 # 2.94455 so use 3
max.minusLogFDR <- 0 # 4.089 so use 4.5
for (i in omics) {
  setwd(file.path(base.path, my.contrast, i, "Differential_expression"))
  # load differential expression results & filter for q <= 0.05
  degs <- read.csv("Differential_expression_results.csv")
  degs <- degs[degs$adj.P.Val <= 0.05, ]
  
  if (i == "global") {
    feature.name <- "Gene"
  } else {
    feature.name <- "SUB_SITE"
  }
  
  for (j in n.features) {
    # identify top features
    top.degs <- degs %>% dplyr::slice_max(Log2FC, n = j/2)
    bottom.degs <- degs %>% dplyr::slice_min(Log2FC, n = j/2)
    top.degs$minusLogFDR <- -log(top.degs$adj.P.Val, 10)
    bottom.degs$minusLogFDR <- -log(top.degs$adj.P.Val, 10)
    bal.degs <- rbind(top.degs, bottom.degs)
    
    # save top differential expression results
    write.csv(bal.degs, paste0("top_", j, "_balanced_diffexp_", feature.name, "_maxFDR0.05.csv"), row.names = FALSE)
    write.csv(bal.degs[,c(feature.name, "Log2FC")], paste0("top_", j, "_balanced_diffexp_", feature.name, "_Log2FC_maxFDR0.05.csv"), row.names = FALSE)
    write.csv(bal.degs[,c(feature.name, "minusLogFDR")], paste0("top_", j, "_balanced_diffexp_", feature.name, "_minusLogFDR_maxFDR0.05.csv"), row.names = FALSE)
    
    # keep track of max abs(Log2FC) and minusLogFDR
    temp.max.abs.Log2FC <- max(abs(bal.degs$Log2FC))
    temp.max.minusLogFDR <- max(bal.degs$minusLogFDR, na.rm = TRUE)
    max.abs.Log2FC <- ifelse(temp.max.abs.Log2FC > max.abs.Log2FC, 
                             temp.max.abs.Log2FC, max.abs.Log2FC)
    max.minusLogFDR <- ifelse(temp.max.minusLogFDR > max.minusLogFDR, 
                              temp.max.minusLogFDR, max.minusLogFDR)
  } 
}

# repeat but looking at individual sample expression of top diffexp proteins/phospho-sites
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
my.contrast <- "Pooled_CD14_Pos_vs_CD34_Pos"
n.features <- c(10, 20, 30, 50)
omics <- c("global", "phospho")
max.expr <- 0 # 6.755 so using 7
expr <- list("global" = global.df, "phospho" = phospho.df)
for (i in omics) {
  setwd(file.path(base.path, my.contrast, i, "Differential_expression"))
  temp.expr <- expr[[i]]
  if (i == "global") {
    feature.name <- "Gene"
  } else {
    feature.name <- "SUB_SITE"
  }
  
  for (j in n.features) {
    # identify top features
    top.degs <- read.csv(paste0("top_", j, "_balanced_diffexp_", feature.name, "_maxFDR0.05.csv"))
    
    # get expression data for top features
    top.expr <- temp.expr[rownames(temp.expr) %in% top.degs[ , feature.name], ]
    
    # save expression data for top features
    write.csv(top.expr, 
              paste0("top_", j, "_balanced_", feature.name, "_maxFDR0.05.csv"), 
              row.names = FALSE)
    
    # keep track of max abs(expression)
    temp.max.expr <- max(abs(top.expr[ , colnames(top.expr) != feature.name]), na.rm = TRUE)
    max.expr <- ifelse(temp.max.expr > max.expr, temp.max.expr, max.expr)
  } 
}

## repeat but only for CD14/CD34 pooled samples - balanced up & down by pulling half from top & half from bottom based on Log2FC
# repeat but looking at individual sample expression of top diffexp proteins/phospho-sites
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis"
my.contrast <- "Pooled_CD14_Pos_vs_CD34_Pos"
n.features <- c(10, 20, 30, 50)
omics <- c("global", "phospho")
max.expr <- 0 # 5.24 so using 5.25

global.df <- global.df[,c("Gene", colnames(global.df)[1:(ncol(global.df)-1)])]
phospho.df <- phospho.df[,c("SUB_SITE", colnames(phospho.df)[1:(ncol(phospho.df)-1)])]

# identify pooled CD14/CD34 samples
CD14.pool <- meta.df[meta.df$Pooled == "CD14_Pos", ]$row.name
CD34.pool <- meta.df[meta.df$Pooled == "CD34_Pos", ]$row.name
sample.pool <- sort(c(CD14.pool, CD34.pool))

# get patient ID annotations
patient.IDs <- meta.df[sample.pool, ]$patient

global.pool.samples <- colnames(global.df)[colnames(global.df) %in% sample.pool]
phospho.pool.samples <- colnames(phospho.df)[colnames(phospho.df) %in% sample.pool]
pool.global <- global.df[ , c("Gene", global.pool.samples)]
pool.phospho <- phospho.df[ , c("SUB_SITE", phospho.pool.samples)]

expr <- list("global" = pool.global, "phospho" = pool.phospho)
for (i in omics) {
  setwd(file.path(base.path, my.contrast, i, "Differential_expression"))
  temp.expr <- expr[[i]]
  if (i == "global") {
    feature.name <- "Gene"
  } else {
    feature.name <- "SUB_SITE"
  }
  
  for (j in n.features) {
    # identify top features
    top.degs <- read.csv(paste0("top_", j, "_balanced_diffexp_", feature.name, "_maxFDR0.05.csv"))
    
    # get expression data for top features
    top.expr <- temp.expr[rownames(temp.expr) %in% top.degs[ , feature.name], ]
    
    # save expression data for top features
    write.csv(top.expr, 
              paste0("top_", j, "_pooled_balanced_", feature.name, "_maxFDR0.05.csv"), 
              row.names = FALSE)
    
    # keep track of max abs(expression)
    temp.max.expr <- max(abs(top.expr[ , colnames(top.expr) != feature.name]), na.rm = TRUE)
    max.expr <- ifelse(temp.max.expr > max.expr, temp.max.expr, max.expr)
  } 
}
