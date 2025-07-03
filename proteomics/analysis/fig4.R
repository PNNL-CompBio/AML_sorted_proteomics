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


av$Log2FC <- -av$Log2FC # make res vs sens instead of sens vs res
mono$Contrast <- "CD14+ vs. CD34+"
av$Contrast <- "Aza+Ven Res vs. Sens"
venn.list <- list("CD14+ vs. CD34+" = mono[mono$adj.P.Val<=0.05,c("Contrast","Gene","Log2FC")]$Gene,
                  "Aza+Ven Res vs. Sens" = av[av$adj.P.Val<=0.05,c("Contrast","Gene","Log2FC")]$Gene)
sigs <- merge(mono[mono$adj.P.Val<=0.05,c("Contrast","Gene","Log2FC")], 
              av[av$adj.P.Val<=0.05,c("Contrast","Gene","Log2FC")], 
              by="Gene", suffixes=c(": CD14+ vs. CD34+", ": Aza+Ven Res vs. Sens"))
# only 2 genes overlapping?
ggvenn::ggvenn(venn.list) # only 2 were sig in ven res vs. sens but both are upregulated in both
