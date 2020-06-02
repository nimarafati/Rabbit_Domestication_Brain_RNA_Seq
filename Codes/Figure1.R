## Figure 1
library(VennDiagram)
library(gplots)
library(ComplexHeatmap)
library(ggplot2)
library(grid)
library(ggplotify)

degs1 <- read.table("../Data/DEGs_Amy_edgeR.tsv",head=T,stringsAsFactors=FALSE)
degs2 <- read.table("../Data/DEGs_Hyp_edgeR.tsv",head=T,stringsAsFactors=FALSE)
degs3 <- read.table("../Data/DEGs_Hipp_edgeR.tsv",head=T,stringsAsFactors=FALSE)
degs4 <- read.table("../Data/DEGs_ParTemp_edgeR.tsv",head=T,stringsAsFactors=FALSE)

v <- venn.diagram(list(Hyp = degs2$TranscriptID, Amy = degs1$TranscriptID, ParTemp = degs4$TranscriptID, Hipp = degs3$TranscriptID),
                  fill = c(alpha("#31688EFF", 0.5), alpha("#35B779FF",0.5), alpha("#FDE725FF", 0.5), alpha("#440154FF", 0.5)),
                  col = c("#31688EFF", "#35B779FF", "#FDE725FF", "#440154FF"),
                  alpha = c(0.5, 0.5, 0.5, 0.5), 
                  lty = 1, 
                  cex = 1,
                  cat.just=list(c(0,1) , c(0,0) , c(0,0) , c(1,1)),
                  cat.fontface = 1,
                  cat.cex = 1,
                  cat.fontfamily = "sans",
                  cat.default.pos = "outer",
                  fontfamily = "sans",
                  filename=NULL,
                  compression = "lzw",
                  lwd = 1)


degs_4reg <- read.table("../Data/DEGs_4regions_edgeR.tsv",head=T,row.names = 1,stringsAsFactors=FALSE)
degs_4reg_fc <- select(degs_4reg,logFC_Amy,logFC_Hyp,logFC_Hipp,logFC_ParTemp)

gb = grid.grabExpr(draw(Heatmap(degs_4reg_fc, col = colorpanel(75,low="#839b5c",mid="#E7E6D5",high="#74325c"), 
                                column_title = expression(paste("Log"[2],"(Fold Change)")))))


pdf("../Figures/Figure1t.pdf", w=10, h=7)
cowplot::plot_grid(v, gb, scale=c(0.85,0.95), rel_widths=c(3, 2.5), labels=c("A","B"), label_size = 20, ncol=2, align="v")
dev.off()