# how does mono sig compare to ven res signature?
# DIA bead
setwd("~/Library/CloudStorage/OneDrive-PNNL/Documents/GitHub/Exp24_patient_cells/proteomics/")
mono <- read.csv("analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/Differential_expression/Differential_expression_results.csv")
ven <- read.csv("analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/Ven_Sensitive_vs_Resistant/global/Differential_expression/Differential_expression_results.csv")
ven$Log2FC <- -ven$Log2FC # make res vs sens instead of sens vs res
mono$Contrast <- "CD14+ vs. CD34+"
ven$Contrast <- "Ven Res vs. Sens"
venn.list <- list("CD14+ vs. CD34+" = mono[mono$adj.P.Val<=0.05,c("Contrast","Gene","Log2FC")]$Gene,
                  "Ven Res vs. Sens" = ven[ven$adj.P.Val<=0.05,c("Contrast","Gene","Log2FC")]$Gene)
sigs <- merge(mono[mono$adj.P.Val<=0.05,c("Contrast","Gene","Log2FC")], 
              ven[ven$adj.P.Val<=0.05,c("Contrast","Gene","Log2FC")], 
              by="Gene", suffixes=c(": CD14+ vs. CD34+", ": Ven Res vs. Sens"))
# only 2 genes overlapping?
ggvenn::ggvenn(venn.list) # only 2 were sig in ven res vs. sens but both are upregulated in both

av <- read.csv("analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/Aza.Ven_Sensitive_vs_Resistant/global/Differential_expression/Differential_expression_results.csv")
# not there? looks like it could be calculated

# what about GSEA?
mono <- read.csv("analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/CD14_Pos_vs_Neg/global/GSEA/GSEA_Hallmark/GSEA_results.csv")
ven <- read.csv("analysis/combined24-27/DIA_2batches_noOutliers_noMSC/Sort Type_Bead/Ven_Sensitive_vs_Resistant/global/GSEA/GSEA_Hallmark/GSEA_results.csv")
ven$NES <- -ven$NES # make res vs sens instead of sens vs res
mono$Contrast <- "CD14+ vs. CD34+"
ven$Contrast <- "Ven Res vs. Sens"
venn.list <- list("CD14+ vs. CD34+" = mono[mono$FDR_q_value<=0.25 & mono$p_value<=0.05,c("Contrast","Feature_set","NES")]$Feature_set,
                  "Ven Res vs. Sens" = ven[ven$FDR_q_value<=0.25 & ven$p_value<=0.05,c("Contrast","Feature_set","NES")]$Feature_set)
sigs <- merge(mono[mono$FDR_q_value<=0.25 & mono$p_value<=0.05,c("Contrast","Feature_set","NES")], 
              ven[ven$FDR_q_value<=0.25 & ven$p_value<=0.05,c("Contrast","Feature_set","NES")], 
              by="Feature_set", suffixes=c(": CD14+ vs. CD34+", ": Ven Res vs. Sens")) # 20 overlapping
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 5, text_size=12)
ggsave("GSEA_venn_mono_venRes.pdf",width=5,height=5)
ggvenn::ggvenn(venn.list, show_percentage = FALSE, set_name_size = 5, text_size=5)
ggsave("GSEA_venn_mono_venRes_size5.pdf",width=5,height=5)
gsea.cor <- cor.test(sigs$`NES: CD14+ vs. CD34+`, sigs$`NES: Ven Res vs. Sens`)
# Pearson's product-moment correlation
# 
# data:  sigs$`NES: CD14+ vs. CD34+` and sigs$`NES: Ven Res vs. Sens`
# t = 20.091, df = 18, p-value = 8.904e-14
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9451077 0.9916057
# sample estimates:
#       cor 
# 0.9784229 

Pearson.est <- gsea.cor$estimate
Pearson.p <- gsea.cor$p.value
stats_pearson <- substitute(
  r == est * "," ~ ~"p" ~ "=" ~ p,
  list(
    est = format(as.numeric(Pearson.est), digits = 3),
    p = format(Pearson.p, digits = 3)
  )
)
maxVal <- max(abs(c(sigs$`NES: Ven Res vs. Sens`, sigs$`NES: CD14+ vs. CD34+`)))
sigs$Feature_set <- sub("HALLMARK_","", sigs$Feature_set)
sigs$inMain <- FALSE
sigs[sigs$Feature_set %in% c("COMPLEMENT","COAGULATION","TNFA_SIGNALING_VIA_NFKB",
                             "INTERFERON_GAMMA_RESPONSE","INFLAMMATORY_RESPONSE",
                             "PROTEIN_SECRETION","XENOBIOTIC_METABOLISM",
                             "OXIDATIVE_PHOSPHORYLATION","E2F_TARGETS",
                             "MYC_TARGETS_V2"),]$inMain <- TRUE
ggplot(sigs, aes(x=`NES: CD14+ vs. CD34+`, y=`NES: Ven Res vs. Sens`)) + geom_point() + theme_minimal() + 
  scale_x_continuous(limits=c(-maxVal, maxVal)) + scale_y_continuous(limits=c(-maxVal, maxVal)) + 
  labs(#x="Lasry RNA-seq: LRCC25", y = "Sorted Proteomics: NCF2", 
    title="Gene Set Enrichment (n = 20)") +
  geom_smooth(method="lm", se=FALSE, linetype="dashed") + 
  ggrepel::geom_label_repel(data=subset(sigs, inMain), aes(label=Feature_set), size=2.5) + 
  theme(plot.title=element_text(hjust=0.5)) +
  ggplot2::geom_text(
    x = -Inf, y = 0.5, vjust = "inward", hjust = "inward",
    colour = "blue", parse = TRUE, 
    label = as.character(as.expression(stats_pearson)), size = 4.5
  ) 
ggsave("Mono_vs_VenRes_GSEA_NES.pdf", width=4, height=4)
