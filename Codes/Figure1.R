## Figure 1
library(gplots)
library(ComplexHeatmap)

degs_4reg <- read.table("../Data/DEGs_4regions_edgeR.tsv",head=T,row.names = 1,stringsAsFactors=FALSE)
degs_4reg_fc <- select(degs_4reg,logFC_Amy,logFC_Hyp,logFC_Hipp,logFC_ParTemp)

pdf("../Figures/Figure1.pdf", w=4, h=6)
Heatmap(degs_4reg_fc, col = colorpanel(75,low="#839b5c",mid="#E7E6D5",high="#74325c"), 
        column_title = expression(paste("Log"[2],"(Fold Change)")), 
        name="Log2(Fold change)", row_names_gp = gpar(fontsize = 8))
dev.off()
