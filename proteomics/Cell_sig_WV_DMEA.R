# trying DMEA using Beat AML not normalized across patients
setwd("~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
source("panSEA_helper_20240508.R")

# load Beat AML data
BeatAML <- load_not_norm_BeatAML_for_DMEA()

# load signature of CD14+ vs. CD34+ cells
synapser::synLogin()
globalFileDIA <- synapser::synGet("syn59429685") # filtered for max FDR of 0.05
global.sig.DIA <- read.csv(globalFileDIA$path) # 1842 proteins

globalFileTMT <- synapser::synGet("syn59438946") # filtered for max FDR of 0.05
global.sig.TMT <- read.csv(globalFileTMT$path) # 703 proteins

phosphoFileTMT <- synapser::synGet("syn59439059") # filtered for max FDR of 0.05
phospho.sig.TMT <- read.csv(phosphoFileTMT$path) # 33 phospho-sites

# filter for genes with 100% coverage in Beat AML database
global100 <- BeatAML$global[,colSums(is.na(BeatAML$global)) == 0] # 6885 proteins out of originally 9413
phospho100 <- BeatAML$phospho[,colSums(is.na(BeatAML$phospho)) == 0] # 1308 sites out of originally 71808

sig <- list("Global_DIA" = global.sig.DIA,
            "Global_TMT" = global.sig.TMT,
            "Phospho_TMT" = phospho.sig.TMT)
# run mDMEA
# mDMEA.results <- panSEA::mDMEA(BeatAML$drug, gmt = BeatAML$gmt,
#                                expression = list(BeatAML$global, BeatAML$global, BeatAML$phospho),
#                                weights = sig, types = names(sig), 
#                                feature.names = c("Gene", "Gene", "SUB_SITE"),
#                                ylab = "Drug AUC")
mDMEA.results <- panSEA::mDMEA(BeatAML$drug, gmt = BeatAML$gmt,
                               expression = list(global100, global100, phospho100),
                               weights = sig, types = names(sig), 
                               feature.names = c("Gene", "Gene", "SUB_SITE"),
                               ylab = "Drug AUC")
global.DIA.DMEA <- extract_DMEA_files(mDMEA.results, 1)
global.TMT.DMEA <- extract_DMEA_files(mDMEA.results, 2)
phospho.TMT.DMEA <- extract_DMEA_files(mDMEA.results, 3)
DMEA.files <- list("Global_DIA" = global.DIA.DMEA,
                   "Global_TMT" = global.TMT.DMEA,
                   "Phospho_TMT" = phospho.TMT.DMEA)

setwd("analysis")
dir.create("Raw_DMEA_CD14_vs_CD34")
setwd("Raw_DMEA_CD14_vs_CD34")
save_to_synapse(DMEA.files)
