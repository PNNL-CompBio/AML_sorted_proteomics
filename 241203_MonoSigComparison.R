# Comparing monocytic signatures

# overview
# 1. Import signatures
# 2. Venn diagram
# 3. Correlation matrix
# 3a. Signature weights
# 3b. Beat AML expression of markers (RNA, protein)
# 3c. Beat AML scores using signatures

base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/"
setwd(base.path)

# path.list: named list of paths to monocytic signatures
# paths should point to CSV files formatted with columns: Gene, Log2FC, P.Value, adj.P.Val
# for Triana, there was no P.Value, so I set P.Value = adj.P.Val
compareSigs <- function(path.list, venn.names = names(path.list), 
                        fname = "Monocyte_vs_progenitor_signatures",
                        fillVals = RColorBrewer::brewer.pal(length(venn.names), "Set2")
                        #,
                        #relevant.genes = c("BCL2", "MAPK14", "CD14", "CD34")
                        ) {
  og.path <- getwd()
  library(plyr)
  dir.create(fname)
  # import signatures and generate venn data
  sigs <- list()
  venn.data <- list()
  sig.venn.data <- list()
  for (i in names(path.list)) {
    sigs[[i]] <- read.csv(path.list[[i]])
    venn.data[[i]] <- sigs[[i]]$Gene
    sig.venn.data[[i]] <- sigs[[i]][sigs[[i]]$adj.P.Val <= 0.05,]$Gene
  }

  # Venn diagram
  setwd(fname)
  my.venn <- ggvenn::ggvenn(venn.data, show_percentage=FALSE, text_size=6, 
                            fill_color=fillVals)
  ggplot2::ggsave("venn_diagram.pdf",my.venn, 
                  width = 7, height = 7)
  
  my.venn <- ggvenn::ggvenn(sig.venn.data, show_percentage=FALSE, text_size=6,
                            fill_color=fillVals)
  ggplot2::ggsave("venn_diagram_significant.pdf",my.venn, 
                  width = 7, height = 7)
  
  # # check if BCL2, MAPK14, CD14, CD34 are in these
  # Gene <- relevant.genes
  # goi <- data.frame(Gene)
  # goi[,names(sigs)] <- NA
  # for (i in 1:length(Gene)) {
  #   for (j in 1:length(sigs)) {
  #     if (Gene[i] %in% sigs[[j]]$Gene) {
  #       goi[i,j+1] <- sigs[[j]][sigs[[j]]$Gene == Gene[i],]$Value
  #     } 
  #   }
  # }
  # write.csv(goi, "relevant_genes.csv", row.names = FALSE)
  # Error in x[[jj]][iseq] <- vjj : replacement has length zero
  # Triana: has CD14 (up for monocytes as expected)
  
  # Correlation matrix
  ## combine signatures into one data frame
  # start with first signature
  sig.df <- sigs[[1]]
  sig.df <- sig.df[,1:2]
  colnames(sig.df)[2] <- names(sigs)[1]
  
  # merge in other signatures one at a time
  for (i in 2:length(sigs)) {
    temp.sig <- sigs[[i]]
    colnames(temp.sig)[2] <- names(sigs)[i]
    sig.df <- merge(sig.df, temp.sig[,1:2], by="Gene")
  }
  write.csv(sig.df, "signature_Log2FC.csv", row.names = FALSE)
  rownames(sig.df) <- sig.df$Gene
  
  if (length(path.list) == 2) {
    comparison <- paste0(names(sigs), collapse = "_and_")
    corr.df <- data.frame(comparison)
    corr.df[,c("Pearson.est", "Pearson.p", "Spearman.est", "Spearman.p")] <- NA
    pearson <- stats::cor.test(sig.df[,2], sig.df[,3], method="pearson")
    corr.df$Pearson.est <- pearson$estimate
    corr.df$Pearson.p <- pearson$p.value
    spearman <- stats::cor.test(sig.df[,2], sig.df[,3], method="spearman")
    corr.df$Spearman.est <- spearman$estimate
    corr.df$Spearman.p <- spearman$p.value
    write.csv(corr.df, "correlation_2signatures.csv", row.names = FALSE)
  }
  
  # create correlation matrix
  corr.mat <- stats::cor(as.matrix(sig.df[,2:ncol(sig.df)]))
  write.csv(corr.mat, "correlations.csv")
  
  # plot correlation matrix
  corr.mat.plot <- ggcorrplot::ggcorrplot(corr.mat)
  ggplot2::ggsave("correlation_matrix.pdf", 
                  corr.mat.plot, width = 7, height = 7)

  
  #try another approach
  # sorted: 2118, van Galen: 50, Triana: 424, Lasry: 10,299
  DEG.df <- data.table::rbindlist(sigs, use.names = TRUE, idcol = "Signature", fill = TRUE)
  DEG.df$minusLogP <- -log(DEG.df$P.Value, base = 10)
  DEG.df$minusLogFDR <- -log(DEG.df$adj.P.Val, base = 10)
  DEG.df$sig <- FALSE
  DEG.df[DEG.df$adj.P.Val < 0.05, ]$sig <- TRUE
  
  # bar plot
  #DEG.df$Significance <- DEG.df$sig
  DEG.df$Direction <- "Upregulated"
  if (nrow(DEG.df[!is.na(DEG.df$Log2FC) & 
                  DEG.df$Log2FC < 0,]) > 0) {
    DEG.df[!is.na(DEG.df$Log2FC) &
             DEG.df$Log2FC < 0,]$Direction <- "Downregulated"
  }
  DEG.df$Significance <- "Not Significant"
  if (any(DEG.df$sig)) {
    DEG.df[DEG.df$sig,]$Significance <- DEG.df[DEG.df$sig,]$Direction
  }
  library(ggplot2)
  bar.plot3log <- ggplot2::ggplot(DEG.df, aes(fill = Significance, x=forcats::fct_infreq(Signature))) + 
    geom_bar(position = "dodge", stat="count") + ggplot2::theme_classic() + 
    scale_y_continuous(trans='log10') +
    ggplot2::xlab("Signature Source") +
    ggplot2::ylab("Number of Quantified Features")
  ggplot2::ggsave("barPlot_logScale.pdf", bar.plot3log, width = 7, height = 7, device = "pdf")
  
  bar.plot4log <- ggplot2::ggplot(DEG.df[DEG.df$Significance != "Not Significant",], 
                                  aes(fill = Direction, x=forcats::fct_infreq(Signature))) + 
    geom_bar(position = "dodge", stat="count") + ggplot2::theme_classic() + 
    scale_y_continuous(trans='log10') +
    ggplot2::xlab("Signature Source") +
    ggplot2::ylab("Number of Differentially Expressed Features")
  ggplot2::ggsave("barPlot_significant_logScale.pdf", bar.plot4log, width = 7, height = 7, device = "pdf")
  
  mean.DEG.df <- plyr::ddply(DEG.df, .(Gene), summarize,
                             mean_Log2FC = mean(Log2FC),
                             sd_Log2FC = sd(Log2FC),
                             Fisher_p = metap::sumlog(na.omit(P.Value))$p,
                             types = paste0(Signature, collapse = ", "),
                             N_types = length(unique(Signature)),
                             N_sig = length(sig[sig]),
                             sig_types = paste0(Signature[sig], collapse = ", "))
  if (length(unique(mean.DEG.df$Fisher_p)) > 1) {
    mean.DEG.df$adj_Fisher_p <- 
      qvalue::qvalue(mean.DEG.df$Fisher_p, pi0=1)$qvalues
  } else {
    mean.DEG.df$adj_Fisher_p <- NA
  }
  write.csv(mean.DEG.df, "compiled_signatures.csv", row.names = FALSE)
  
  #DEG.df$sig <- NULL
  write.csv(DEG.df, "signatures.csv", row.names = FALSE)
  setwd(og.path)
}

#### 1. Mono vs Other ####
sig.paths <- list("Sorted" = "analysis/DIA/no_filter/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "van Galen" = "data/externalSignatures/formatted/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Other_protein-coding.csv",
                  "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_otherCellTypes_differentialExpression.csv",
                  "Lasry" = "data/externalSignatures/formatted/Lasry_metadata_clustering_w_header_upd_csv--Cell_type_identity--CD14pos--cluster--wilcoxon.csv")
compareSigs(sig.paths, fname = "Monocyte_vs_other_signatures")

sig.paths <- list("Sorted" = "analysis/DIA/MSC_Non_MSC/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "van Galen" = "data/externalSignatures/formatted/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Other_protein-coding.csv",
                  "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_otherCellTypes_differentialExpression.csv",
                  "Lasry" = "data/externalSignatures/formatted/Lasry_metadata_clustering_w_header_upd_csv--Cell_type_identity--CD14pos--cluster--wilcoxon.csv")
compareSigs(sig.paths, fname="Monocyte_vs_other_exceptSorted_signatures")

#### 2. Mono vs. progenitor ####
sig.paths <- list("Sorted" = "analysis/DIA/MSC_Non_MSC/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "van Galen" = "data/externalSignatures/formatted/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_noNA.csv",
                  "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
                  "Lasry" = "data/externalSignatures/formatted/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HSC.csv")
compareSigs(sig.paths)

sig.paths <- list("Sorted" = "analysis/DIA_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_noNA.csv",
                  "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
                  "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HSC_protein-coding.csv")
compareSigs(sig.paths, fname = "Monocyte_vs_progenitor_signatures_beadOnly")

sig.paths <- list("Sorted" = "analysis/DIA_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "Lasry" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_Lasry_AML_CD14PosMonocyte_vs_HSC_protein-coding.csv",
                  "Triana" = "data/externalSignatures/formatted/Triana_RNA_AML_100PercentCells_Classical-Monocytes_vs_HSCs-and-MPPs_differentialExpression.csv",
                  "van Galen" = "data/externalSignatures/formatted/notFilteredForMalignant/Differential_expression_van_Galen_AML_D0_Mono-like_vs_Prog-like_noNA.csv")
compareSigs(sig.paths, fname = "Monocyte_vs_progenitor_signatures_beadOnly_2025-01-20")

#### 3. flow type ####
# since MSC was only flow sorted, remove MSC from sort type comparisons
sig.paths <- list("Bead vs. Flow" = "analysis/DIA_noMSC/no_filter/Sort Type_Bead_vs_Flow/Differential_expression/Differential_expression_results.csv",
                  "CD14 vs. CD34" = "analysis/DIA_noMSC/no_filter/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "Bead: CD14 vs. CD34" = "analysis/DIA_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "Flow: CD14 vs. CD34" = "analysis/DIA_noMSC/Sort Type_Flow/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv")
compareSigs(sig.paths, fname = "Monocyte_vs_other_signatures_and_bead_vs_flow")

sig.paths <- list("Bead" = "analysis/DIA_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv",
                  "Flow" = "analysis/DIA_noMSC/Sort Type_Flow/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv")
compareSigs(sig.paths, fname = "Monocyte_vs_other_signatures_beadOrFlow")
