# correlations between sorted expression and bulk AUC
library(synapser);library(ggplot2);library(DMEA)
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis")
synapser::synLogin()
base.path <- getwd()
dir.create("correlations")
setwd("correlations")

#### prep data ####
# get sorted data
dia.wo.out <- readRDS("DIA_2batches_noOutliers.rds")

# get bulk drug sensitivity data
drug.BeatAML <- read.csv(synapser::synGet("syn65677650")$path) # other drug sensitivity data: syn51674470
ven <- drug.BeatAML[drug.BeatAML$Specimen..Lab.ID %in% dia.wo.out$meta$patient & 
                      drug.BeatAML$Inhibitor.Panel.Definition..Drug=="Venetoclax",
                            c("Specimen..Lab.ID","Probit.Interpretation..Area.Under.Curve")] # 20 / 23
azaVen <- drug.BeatAML[drug.BeatAML$Specimen..Lab.ID %in% dia.wo.out$meta$patient & 
                      drug.BeatAML$Inhibitor.Panel.Definition..Drug=="Azacitidine - Venetoclax",
                    c("Specimen..Lab.ID","Probit.Interpretation..Area.Under.Curve")] # 20 / 23
# any samples lost?
lost.samples <- unique(dia.wo.out$meta$patient[!(dia.wo.out$meta$patient %in% sorted.drug$Specimen..Lab.ID)])
# 3 / 23 patients lost: "16-01184" "21-00839" "22-00251"
lost.samples.ven <- unique(dia.wo.out$meta$patient[!(dia.wo.out$meta$patient %in% ven$Specimen..Lab.ID)])
# 3 / 23 patients lost: "16-01184" "21-00839" "22-00251"
lost.samples.av <- unique(dia.wo.out$meta$patient[!(dia.wo.out$meta$patient %in% azaVen$Specimen..Lab.ID)])
# 3 / 23 patients lost: "16-01184" "21-00839" "22-00251"

# add data from another file to get info for these 3 patients
other.pts <- read.csv(synapser::synGet("syn53627410")$path)
colnames(ven) <- c("sample_id","auc")
colnames(azaVen) <- c("sample_id","auc")
other.ven <- other.pts[other.pts$patient %in% lost.samples,c("patient","Ven_AUC")]
other.av <- other.pts[other.pts$patient %in% lost.samples,c("patient","Aza.Ven_AUC")]
colnames(other.ven) <- colnames(ven)
colnames(other.av) <- colnames(ven)
ven <- rbind(ven, other.ven)
azaVen <- rbind(azaVen, other.av)

# format proteomics for correlations
prot <- as.data.frame(t(dia.wo.out$global))
prot$sample_id <- rownames(prot)
prot <- prot[,c("sample_id",colnames(prot)[1:(ncol(prot)-1)])]
prot$sample_id <- sub("PTRC_","",prot$sample_id)
prot$sample_id <- sub("_.*","",prot$sample_id)
prot$sample_id <- sub("[.]","-",prot$sample_id)
id.map <- unique(prot$sample_id[grepl("X",prot$sample_id)])
names(id.map) <- c("16-01184","17-01060","18-00103","18-00105","19-00074",
                   "21-00432","21-00839","22-00117","22-00251","22-00571")
for (i in 1:length(id.map)) {
  prot[prot$sample_id==id.map[i],]$sample_id <- names(id.map)[i]
}

# filter proteomics for sorting method, cell type
prot.bead <- prot[!grepl("_f_",rownames(prot)) & 
                    !grepl("Flow",rownames(prot), ignore.case=TRUE),]
prot.flow <- prot[grepl("_f_",rownames(prot)) | 
                    grepl("Flow",rownames(prot), ignore.case=TRUE),]
prot.cd14 <- prot[grepl("cd14",rownames(prot),ignore.case=TRUE),]
prot.cd14.bead <- prot.bead[grepl("cd14",rownames(prot.bead),ignore.case=TRUE),]
prot.cd14.flow <- prot.flow[grepl("cd14",rownames(prot.flow),ignore.case=TRUE),]
prot.cd34 <- prot[grepl("cd34",rownames(prot),ignore.case=TRUE),]
prot.cd34.bead <- prot.bead[grepl("cd34",rownames(prot.bead),ignore.case=TRUE),]
prot.cd34.flow <- prot.flow[grepl("cd34",rownames(prot.flow),ignore.case=TRUE),]
prot.msc.flow <- prot.flow[grepl("msc",row.names(prot.flow),ignore.case=TRUE),]

#### run correlations ####
inputs <- list("Overall" = prot, "Bead" = prot.bead, "Flow" = prot.flow,
               "CD14" = prot.cd14, "CD14_Bead" = prot.cd14.bead,
               "CD14_Flow" = prot.cd14.flow, "CD34" = prot.cd34,
               "CD34_Bead" = prot.cd34.bead, "CD34_Flow" = prot.cd34.flow,
               "MSC_Flow" = prot.msc.flow)
for (i in names(inputs)) {
  # prep data
  ven.prot <- merge(ven, inputs[[i]], by="sample_id")
  ven.prot$auc <- as.numeric(ven.prot$auc)
  ven.prot <- ven.prot[!is.na(ven.prot$auc),]
  av.prot <- merge(azaVen, inputs[[i]], by="sample_id")
  av.prot$auc <- as.numeric(av.prot$auc)
  av.prot <- av.prot[!is.na(av.prot$auc),]
  
  # run correlation
  ven.prot.corr <- DMEA::rank_corr(ven.prot, variable="Protein",value="normAbudance") # no q<0.05
  av.prot.corr <- DMEA::rank_corr(av.prot, variable="Protein",value="normAbudance") # no q<0.05
  
  # save results
  write.csv(ven.prot.corr$result,
            paste0("venAUC_correlationsWithSortedProteomicsDIA_",i,".csv"),row.names=FALSE)
  write.csv(av.prot.corr$result,
            paste0("azaVenAUC_correlationsWithSortedProteomicsDIA_",i,".csv"),row.names=FALSE)
  if (is.list(ven.prot.corr$scatter.plots)) {
    ggsave(paste0("venAUC_correlationsWithSortedProteomicsDIA_",i,".pdf"),
           ven.prot.corr$scatter.plots, width=5,height=5)
  }
  if (is.list(av.prot.corr$scatter.plots)) {
    ggsave(paste0("azaVenAUC_correlationsWithSortedProteomicsDIA_",i,".pdf"),
           av.prot.corr$scatter.plots, width=5,height=5)
  }
}

#### histograms ####
setwd(file.path(base.path,"correlations"))
dir.create("histograms")
setwd("histograms")
for (i in names(inputs)) {
  # prep data
  setwd(file.path(base.path,"correlations"))
  ven.prot.corr <- read.csv(paste0("venAUC_correlationsWithSortedProteomicsDIA_",i,".csv"))
  av.prot.corr <- read.csv(paste0("azaVenAUC_correlationsWithSortedProteomicsDIA_",i,".csv"))
  gsea.inputs <- list("Ven" = ven.prot.corr[,c("Protein","Spearman.est","Spearman.q")],
                      "AzaVen" = av.prot.corr[,c("Protein","Spearman.est","Spearman.q")])
  
  for (j in names(gsea.inputs)) {
    setwd(file.path(base.path,"correlations","histograms"))
    dir.create(j)
    setwd(j)
    temp.df <- na.omit(gsea.inputs[[j]])
    temp.df$Significance <- "Adjusted p > 0.05"
    if (any(temp.df$Spearman.q <= 0.05)) {
      temp.df[which(temp.df$Spearman.q <= 0.05),]$Significance <- "Adjusted p <= 0.05"
    }
    n.distinct <- length(unique(temp.df$Spearman.est))
    n.ties <- nrow(temp.df) - n.distinct
    perc.ties <- round(n.ties * 100 / nrow(temp.df),0)
    temp.plot <- ggplot2::ggplot(temp.df, aes(x=Spearman.est, 
                                              group=Significance, 
                                              fill=Significance)) +
      geom_histogram(position="identity", alpha=0.5) + theme_classic() +
      scale_fill_manual(values = c("#00BFC4", "#F8766D"), breaks = c("Adjusted p <= 0.05", "Adjusted p > 0.05")) +
      xlab("Spearman Correlation Estimates") + 
      ggtitle(paste0(j, " AUC correlations with DIA\n",i," Protein Expression\n(", perc.ties, "% tied)")) +
      theme(axis.text = element_text(size=12), axis.title = element_text(size=16), legend.title = element_text(size=16),
            legend.text=element_text(size=12), title = element_text(size = 24, hjust = 0.5),
            legend.position = "bottom")
    ggsave(paste0(j,"AUC_correlationsWithSortedProteomicsDIA_",i,".pdf"),
           temp.plot, width=7, height=7)
  }
}

#### run GSEA ####
setwd(file.path(base.path,"correlations"))
dir.create("GSEA")
setwd("GSEA")
for (i in names(inputs)) {
  # prep data
  setwd(file.path(base.path,"correlations"))
  ven.prot.corr <- read.csv(paste0("venAUC_correlationsWithSortedProteomicsDIA_",i,".csv"))
  av.prot.corr <- read.csv(paste0("azaVenAUC_correlationsWithSortedProteomicsDIA_",i,".csv"))
  gsea.inputs <- list("Ven" = ven.prot.corr[,c("Protein","Spearman.est")],
                      "AzaVen" = av.prot.corr[,c("Protein","Spearman.est")])
  
  # run GSEA with ties
  temp.gsea <- panSEA::mGSEA(gsea.inputs, 
                             feature.names = rep("Protein", length(gsea.inputs)),
                             rank.var = rep("Spearman.est", length(gsea.inputs)),
                             gmt=as.list(rep("msigdb_Homo sapiens_HS_H", length(gsea.inputs))),
                             ties=TRUE, types=names(gsea.inputs))
  
  # save results
  for (j in names(temp.gsea$all.results)) {
    setwd(file.path(base.path,"correlations","GSEA"))
    dir.create(j)
    setwd(j)
    write.csv(temp.gsea$all.results[[j]]$result,
              paste0(j,"AUC_gsea_WithSortedProteomicsDIA_",i,".csv"),row.names=FALSE)
    write.csv(temp.gsea$all.results[[j]]$result.w.ties,
              paste0(j,"AUC_gseaWithTies_WithSortedProteomicsDIA_",i,".csv"),row.names=FALSE)
    ggsave(paste0(j,"AUC_gseaVolcano_WithSortedProteomicsDIA_",i,".pdf"),
           temp.gsea$all.results[[j]]$volcano.plot, width=5, height=5)
    ggsave(paste0(j,"AUC_gseaBar_WithSortedProteomicsDIA_",i,".pdf"),
           temp.gsea$all.results[[j]]$bar.plot, width=5, height=5)
    ggsave(paste0(j,"AUC_gseaDot_WithSortedProteomicsDIA_",i,".pdf"),
           temp.gsea$all.results[[j]]$dot.plot, width=5, height=5)
    ggsave(paste0(j,"AUC_gseaDotSD_WithSortedProteomicsDIA_",i,".pdf"),
           temp.gsea$all.results[[j]]$dot.sd, width=5, height=5)
  }
}


