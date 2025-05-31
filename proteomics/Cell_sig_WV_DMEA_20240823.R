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

# run DMEA with all features adj p <= 0.05
setwd("analysis")
dir.create("Raw_DMEA_CD14_vs_CD34")
setwd("Raw_DMEA_CD14_vs_CD34")
base.path <- "~/OneDrive - PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/analysis/Raw_DMEA_CD14_vs_CD34"
# run mDMEA
mDMEA.results <- panSEA::mDMEA(BeatAML$drug, gmt = BeatAML$gmt,
                               expression = list(global100, global100, phospho100),
                               weights = sig, types = names(sig), 
                               feature.names = c("Gene", "Gene", "SUB_SITE"),
                               ylab = "Drug AUC", xlab = "CD14+ vs. CD34+ Score")
global.DIA.DMEA <- extract_DMEA_files(mDMEA.results, 1)
global.TMT.DMEA <- extract_DMEA_files(mDMEA.results, 2)
phospho.TMT.DMEA <- extract_DMEA_files(mDMEA.results, 3)
DMEA.files <- list("Global_DIA" = global.DIA.DMEA,
                   "Global_TMT" = global.TMT.DMEA,
                   "Phospho_TMT" = phospho.TMT.DMEA)
save_to_synapse(DMEA.files)

# run DMEA with top i features
n.features <- c(1, 2, 5, 10, 25, 50, 100, 250, 500, 1000, 2000)
corr.results <- list("Global_DIA" = list(),
                     "Global_TMT" = list(),
                     "Phospho_TMT" = list())
omics <- list("Global_DIA" = sig[["Global_DIA"]],
              "Global_TMT" = sig[["Global_TMT"]],
              "Phospho_TMT" = sig[["Phospho_TMT"]])
for (j in names(omics)) {
  sig <- omics[[j]]
  if (grepl("Phospho", j)) {
    expr <- phospho100
  } else {
    expr <- global100
  }
  omics.corr <- list()
  omics.WV <- list()
  for (i in n.features) {
    top.sig <- sig %>% slice_max(abs(Log2FC), n=i)
    if (nrow(sig) < i) {temp.n <- nrow(sig)} else {temp.n <- i}
    setwd(base.path)
    dir.create(paste0("top_", temp.n, "_features_by_absLog2FC"))
    setwd(paste0("top_", temp.n, "_features_by_absLog2FC"))
    # run mDMEA
    DMEA.results <- DMEA::DMEA(BeatAML$drug, gmt = BeatAML$gmt,
                                   expression = expr,
                                   weights = top.sig, 
                                   ylab = "Drug AUC", xlab = "CD14+ vs. CD34+ Score",
                                position.y = "max")
    temp.WV <- DMEA.results$WV.scores
    temp.corr <- DMEA.results$corr.result
    omics.corr[[as.character(temp.n)]] <- temp.corr[temp.corr$Drug == "Venetoclax" | 
                                                      temp.corr$Drug == "Azacytidine - Venetoclax",]
    mDMEA.format <- list()
    mDMEA.format[["all.results"]] <- list(DMEA.results)
    DMEA.result.files <- extract_DMEA_files(mDMEA.format, 1)
    DMEA.files <- list(DMEA.result.files)
    names(DMEA.files) <- j
    save_to_synapse(DMEA.files)  
    if (i >= nrow(sig)) {break}
  } 
  omics.corr.df <- data.table::rbindlist(omics.corr, use.names = TRUE, idcol = "N_features")
  corr.results[[j]] <- omics.corr.df
}
corr.results.df <- data.table::rbindlist(corr.results, use.names = TRUE, idcol = "Omics")
setwd(base.path)
write.csv(corr.results.df, "Correlation_results_AzaVen_Ven.csv", row.names = FALSE)

# plot Aza, Ven results for each omics type
drugs <- c("Azacytidine - Venetoclax", "Venetoclax")

# set theme for dot plot
library(ggplot2)
bg.theme <- ggplot2::theme(
  legend.background = element_rect(), legend.position = "top",
  legend.text = element_text(size = 14),
  legend.key = element_blank(),
  legend.title = element_text(size = 16),
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 16, colour = "black"),
  axis.title.y = element_text(size = 20),
  axis.text.y = element_text(size = 16, colour = "black"),
  plot.title = element_text(
    lineheight = .8, face = "bold", size = 36
  ),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(fill = NA),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks.x = element_line(colour = "black"),
  axis.ticks.y = element_line(colour = "black")
)

# for (j in names(omics)) {
#   omics.results <- corr.results.df[corr.results.df$Omics == j,]
#   # dot plot
#   dot.plot <- ggplot2::ggplot(
#     omics.results,
#     ggplot2::aes(
#       x = N_features, y = Drug, color = Pearson.est,
#       size = -log10(Pearson.q)
#     )
#   ) +
#     ggplot2::geom_point() +
#     # ggplot2::scale_y_discrete(limits = mean.DEG.df[
#     #   mean.DEG.df$feature %in% top.DEG.df$feature, ]$feature) +
#     viridis::scale_color_viridis() +
#     bg.theme +
#     ggplot2::labs(
#       x = "Number of Features",
#       #y = "Drug",
#       color = "Pearson Correlation", size = "-log(adjusted p-value)"
#     )
#   ggplot2::ggsave(paste0(j, "_Pearson_correlation_AzaVen_Ven_vs_N_features_dot_plot.pdf"))
#   
#   for (i in drugs) {
#     drug.results <- omics.results[omics.results$Drug == i,]
#     dot.plot <- ggplot2::ggplot(
#       drug.results,
#       ggplot2::aes(
#         x = as.numeric(N_features), y = Pearson.est,
#         size = -log10(Pearson.q)
#       )
#     ) +
#       ggplot2::geom_point() +
#       bg.theme + ylim(0, 1) +
#       ggplot2::labs(
#         x = "Number of Features",
#         y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
#       )
#     ggplot2::ggsave(paste0(j, "_", i, "_vs_N_features_scatter_plot_01.pdf"))
#     
#     dot.plot <- ggplot2::ggplot(
#       drug.results,
#       ggplot2::aes(
#         x = as.numeric(N_features), y = Pearson.est,
#         size = -log10(Pearson.q)
#       )
#     ) +
#       ggplot2::geom_point() +
#       bg.theme + ylim(min(corr.results.df$Pearson.est), max(corr.results.df$Pearson.est)) +
#       ggplot2::labs(
#         x = "Number of Features",
#         y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
#       )
#     ggplot2::ggsave(paste0(j, "_", i, "_vs_N_features_scatter_plot_minMax.pdf"))
#   }
# }
# 
# # plot DIA & TMT global results on same scatter plot
# 
# corr.results.df$Method <- sub("Global_", "", corr.results.df$Omics)
# for (i in drugs) {
#   drug.results <- na.omit(corr.results.df[corr.results.df$Drug == i,])
#   dot.plot <- ggplot2::ggplot(
#     drug.results,
#     ggplot2::aes(
#       x = as.numeric(N_features), y = Pearson.est,
#       size = -log10(Pearson.q), shape = Method, color = Method 
#     )
#   ) +
#     ggplot2::geom_point() +
#     bg.theme + ylim(0,1) +
#     ggplot2::labs(
#       x = "Number of Features",
#       y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
#     )
#   ggplot2::ggsave(paste0("Global_", i, "_vs_N_features_scatter_plot_01_withColor.pdf"))
#   dot.plot <- ggplot2::ggplot(
#     drug.results,
#     ggplot2::aes(
#       x = as.numeric(N_features), y = Pearson.est,
#       size = -log10(Pearson.q), shape = Method, color = Method
#     )
#   ) +
#     ggplot2::geom_point() +
#     bg.theme + ylim(min(corr.results.df$Pearson.est), max(corr.results.df$Pearson.est)) +
#     ggplot2::labs(
#       x = "Number of Features",
#       y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
#     )
#   ggplot2::ggsave(paste0("Global_", i, "_vs_N_features_scatter_plot_minMax_withColor.pdf"))
# }

# both drugs on same plot
dot.plot <- ggplot2::ggplot(
  corr.results.df,
  ggplot2::aes(
    x = as.numeric(N_features), y = Pearson.est,
    size = -log10(Pearson.q), shape = Method, color = Drug
  )
) +
  ggplot2::geom_point() +
  bg.theme + ylim(0,1) +
  ggplot2::labs(
    x = "Number of Features",
    y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
  ) + theme(legend.position = "right")
ggplot2::ggsave("Global_vs_N_features_scatter_plot_01_withColor_v2_wider.pdf", width = 12)
dot.plot <- ggplot2::ggplot(
  corr.results.df,
  ggplot2::aes(
    x = as.numeric(N_features), y = Pearson.est,
    size = -log10(Pearson.q), shape = Method, color = Drug
  )
) +
  ggplot2::geom_point() +
  bg.theme + ylim(min(corr.results.df$Pearson.est), max(corr.results.df$Pearson.est)) +
  ggplot2::labs(
    x = "Number of Features",
    y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
  ) + theme(legend.position = "right")
ggplot2::ggsave(paste0("Global_vs_N_features_scatter_plot_minMax_withColor_v2_wider.pdf"), width = 12)

# log-scaled x-axis
dot.plot <- ggplot2::ggplot(
  corr.results.df,
  ggplot2::aes(
    x = as.numeric(N_features), y = Pearson.est,
    size = -log10(Pearson.q), shape = Method, color = Drug
  )
) + scale_x_continuous(trans="log10") +
  ggplot2::geom_point() +
  bg.theme + ylim(0,1) +
  ggplot2::labs(
    x = "Number of Features",
    y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
  ) + theme(legend.position = "right")
ggplot2::ggsave("Global_vs_LogN_features_scatter_plot_01_withColor_v2_wider.pdf", width = 12)
dot.plot <- ggplot2::ggplot(
  corr.results.df,
  ggplot2::aes(
    x = as.numeric(N_features), y = Pearson.est,
    size = -log10(Pearson.q), shape = Method, color = Drug
  )
) + scale_x_continuous(trans="log10") +
  ggplot2::geom_point() +
  bg.theme + ylim(min(corr.results.df$Pearson.est), max(corr.results.df$Pearson.est)) +
  ggplot2::labs(
    x = "Number of Features",
    y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
  ) + theme(legend.position = "right")
ggplot2::ggsave(paste0("Global_vs_LogN_features_scatter_plot_minMax_withColor_v2_wider.pdf"), width = 12)

# now select proteins based on min adj p-value
for (j in names(omics)) {
  sig <- omics[[j]]
  if (grepl("Phospho", j)) {
    expr <- phospho100
  } else {
    expr <- global100
  }
  omics.corr <- list()
  for (i in n.features) {
    top.sig <- sig %>% slice_min(adj.P.Val, n=i)
    if (nrow(sig) < i) {temp.n <- nrow(sig)} else {temp.n <- i}
    setwd(base.path)
    dir.create(paste0("top_", temp.n, "_features_by_adjPval"))
    setwd(paste0("top_", temp.n, "_features_by_adjPval"))
    # run mDMEA
    DMEA.results <- DMEA::DMEA(BeatAML$drug, gmt = BeatAML$gmt,
                               expression = expr,
                               weights = top.sig, 
                               ylab = "Drug AUC", xlab = "CD14+ vs. CD34+ Score",
                               position.y = "max")
    temp.corr <- DMEA.results$corr.result
    omics.corr[[as.character(temp.n)]] <- temp.corr[temp.corr$Drug == "Venetoclax" | 
                                                      temp.corr$Drug == "Azacytidine - Venetoclax",]
    mDMEA.format <- list()
    mDMEA.format[["all.results"]] <- list(DMEA.results)
    DMEA.result.files <- extract_DMEA_files(mDMEA.format, 1)
    DMEA.files <- list(DMEA.result.files)
    names(DMEA.files) <- j
    save_to_synapse(DMEA.files)  
    if (i >= nrow(sig)) {break}
  } 
  omics.corr.df <- data.table::rbindlist(omics.corr, use.names = TRUE, idcol = "N_features")
  corr.results[[j]] <- omics.corr.df
}
corr.results.df <- data.table::rbindlist(corr.results, use.names = TRUE, idcol = "Omics")
setwd(base.path)
write.csv(corr.results.df, "Correlation_results_AzaVen_Ven_by_adjPval.csv", row.names = FALSE)

# both drugs on same plot
corr.results.df$Method <- NA
for (i in 1:nrow(corr.results.df)) {
  corr.results.df$Method[i] <- stringr::str_split(corr.results.df$Omics[i], "_")[[1]][2] 
}
corr.results.df <- na.omit(corr.results.df)
global.results.df <- corr.results.df[corr.results.df$Omics == "Global_DIA" | 
                                    corr.results.df$Omics == "Global_TMT",]
dot.plot <- ggplot2::ggplot(
  global.results.df,
  ggplot2::aes(
    x = as.numeric(N_features), y = Pearson.est,
    size = -log10(Pearson.q), shape = Method, color = Drug
  )
) +
  ggplot2::geom_point() +
  bg.theme + ylim(0,1) +
  ggplot2::labs(
    x = "Number of Features",
    y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
  ) + theme(legend.position = "right")
ggplot2::ggsave("Global_vs_N_features_by_adjPval_scatter_plot_01.pdf", width = 12)
dot.plot <- ggplot2::ggplot(
  global.results.df,
  ggplot2::aes(
    x = as.numeric(N_features), y = Pearson.est,
    size = -log10(Pearson.q), shape = Method, color = Drug
  )
) +
  ggplot2::geom_point() +
  bg.theme + ylim(min(global.results.df$Pearson.est), max(global.results.df$Pearson.est)) +
  ggplot2::labs(
    x = "Number of Features",
    y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
  ) + theme(legend.position = "right")
ggplot2::ggsave(paste0("Global_vs_N_features_by_adjPval_scatter_plot_minMax.pdf"), width = 12)

# log-scaled x-axis
dot.plot <- ggplot2::ggplot(
  global.results.df,
  ggplot2::aes(
    x = as.numeric(N_features), y = Pearson.est,
    size = -log10(Pearson.q), shape = Method, color = Drug
  )
) + scale_x_continuous(trans="log10") +
  ggplot2::geom_point() +
  bg.theme + ylim(0,1) +
  ggplot2::labs(
    x = "Number of Features",
    y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
  ) + theme(legend.position = "right")
ggplot2::ggsave("Global_vs_LogN_features_by_adjPval_scatter_plot_01.pdf", width = 12)
dot.plot <- ggplot2::ggplot(
  global.results.df,
  ggplot2::aes(
    x = as.numeric(N_features), y = Pearson.est,
    size = -log10(Pearson.q), shape = Method, color = Drug
  )
) + scale_x_continuous(trans="log10") +
  ggplot2::geom_point() +
  bg.theme + ylim(min(global.results.df$Pearson.est), max(global.results.df$Pearson.est)) +
  ggplot2::labs(
    x = "Number of Features",
    y = "Pearson Correlation Estimate", size = "-log(adjusted p-value)"
  ) + theme(legend.position = "right")
ggplot2::ggsave(paste0("Global_vs_LogN_features_by_adjPval_scatter_plot_minMax.pdf"), width = 12)

Log2FC.results <- read.csv("Correlation_results_AzaVen_Ven.csv")
global.results.df$Selection <- "Adjusted P-value"
global.results.df$N_features <- as.numeric(global.results.df$N_features)
Log2FC.results$Selection <- "Absolute Log2FC"
Log2FC.results$N_features <- as.numeric(Log2FC.results$N_features)
all.results <- merge(global.results.df, Log2FC.results, by=c("Omics", "N_features", "Drug"), suffixes = c("_adjPval", "_absLog2FC"))
# run correlation
corr.DIA.TMT <- cor.test(all.results$Pearson.est_absLog2FC, all.results$Pearson.est_adjPval)
p <- corr.DIA.TMT$p.value
est <- as.numeric(corr.DIA.TMT$estimate)
all.results$Pearson.est_absLog2FC <- as.numeric(all.results$Pearson.est_absLog2FC)
all.results$Pearson.est_adjPval <- as.numeric(all.results$Pearson.est_adjPval)
corr.scatter.plot <- corr.scatter(all.results, "Pearson.est_absLog2FC", "Pearson.est_adjPval", "Proteins Selected by Adjusted P-value",
                                  "Proteins Selected by Absolute Log2FC", "Pearson Correlation as Measure of Prediction Accuracy", 
                                  p, est, shape = "Method", color = "Drug", symmetric = FALSE)
df <- all.results
rank.var <- "Pearson.est_absLog2FC"
value <- "Pearson.est_adjPval"
xlab <- "Proteins Selected by Absolute Log2FC"
ylab <- "Proteins Selected by Adjusted P-value"
title <- "Pearson Correlation as Measure of Prediction Accuracy"
Pearson.p <- p
Pearson.est <- est
shape <- "Method"
color <- "Drug"


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
pos.x <- min(df[, c(rank.var)])
pos.y <- max(df[, c(value)])

stats_pearson <- substitute(
  r == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(Pearson.est, digits = 3),
    p = format(Pearson.p, digits = 3)
  )
)
df$Pearson.est_adjPval <- as.numeric(df$Pearson.est_adjPval)
df$Pearson.est_absLog2FC <- as.numeric(df$Pearson.est_absLog2FC)
scatter.plot <- ggplot2::ggplot(data = df,
                                aes_string(x = rank.var, y = value)) +
  ggplot2::geom_point(aes_string(
    shape = shape, color = color)) +
  ggplot2::labs(x = xlab, y = ylab) +
  ggplot2::ggtitle(title) + 
  #ylim(0,1) + xlim(0,1) +
  ggplot2::geom_smooth(method = "lm", size = 1.5,
                       linetype = "solid", color = "blue",
                       se = TRUE, na.rm = TRUE) +
  # ggplot2::geom_text(
  #   x = pos.x, y = pos.y, vjust = "inward", hjust = "inward",
  #   colour = "blue", parse = TRUE,
  #   label = as.character(as.expression(stats_pearson)), size = 8
  # ) +
  ng.theme +
  bg.theme + theme_minimal()
write.csv(all.results, "Correlation_results_AzaVen_Ven_by_adjPval_or_absLog2FC.csv", row.names = FALSE)
