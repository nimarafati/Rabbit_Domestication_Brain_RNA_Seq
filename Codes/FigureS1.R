## Figure S1
library(gplots)
library(RColorBrewer)

#heatmap_amy
fc_tmp <- read.table("../Data/DEGs_Amy_edgeR.tsv",head=T,row.names = 1) %>% select(c(-logFC,-logCPM,-PValue,-FDR))
cr = cor(fc_tmp, method='spearman')
data.TMM.matrix = log2(fc_tmp+1)
data.TMM.matrix = t(scale(t(data.TMM.matrix), scale=F)) # center rows, mean substracted
gene_dist = dist(data.TMM.matrix, method='euclidean')
hc_genes = hclust(gene_dist, method='complete')
hc_samples = hclust(as.dist(1-cr), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);
partition_colors = brewer.pal(length(unique(gene_partition_assignments)), "Accent")
gene_colors = partition_colors[gene_partition_assignments]
quantBrks = quantile(data.TMM.matrix, c(0.03, 0.97))
myheatcol = colorpanel(75,low="#3C439A",mid="#E7E6D5",high="#CB4335")

pdf("../Figures/FigureS1_Amy.pdf",height=8,width=10)
heatmap.2(data.TMM.matrix, dendrogram='both', Rowv=as.dendrogram(hc_genes),
          col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, cexCol=1, lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1.5,5),lwid=c(1.5,0.2,2.5,1.5), margins=c(10,2), breaks=seq(quantBrks[1], quantBrks[2], length=76), labRow=F,main="Amygdala", key.title="", key.xlab = "Log2(TMM+1)")
dev.off()

#heatmap_hyp
fc_tmp <- read.table("../Data/DEGs_Hyp_edgeR.tsv",head=T,row.names = 1) %>% select(c(-logFC,-logCPM,-PValue,-FDR))
cr = cor(fc_tmp, method='spearman')
data.TMM.matrix = log2(fc_tmp+1)
data.TMM.matrix = t(scale(t(data.TMM.matrix), scale=F)) # center rows, mean substracted
gene_dist = dist(data.TMM.matrix, method='euclidean')
hc_genes = hclust(gene_dist, method='complete')
hc_samples = hclust(as.dist(1-cr), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);
partition_colors = brewer.pal(length(unique(gene_partition_assignments)), "Accent")
gene_colors = partition_colors[gene_partition_assignments]
quantBrks = quantile(data.TMM.matrix, c(0.03, 0.97))
myheatcol = colorpanel(75,low="#3C439A",mid="#E7E6D5",high="#CB4335")

pdf("../Figures/FigureS1_Hyp.pdf",height=8,width=10)
heatmap.2(data.TMM.matrix, dendrogram='both', Rowv=as.dendrogram(hc_genes),
          col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, cexCol=1, lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1.5,5),lwid=c(1.5,0.2,2.5,1.5), margins=c(10,2), breaks=seq(quantBrks[1], quantBrks[2], length=76), labRow=F,main="Hypothalamus", key.title="", key.xlab = "Log2(TMM+1)")
dev.off()

#heatmap_hip
fc_tmp <- read.table("../Data/DEGs_Hipp_edgeR.tsv",head=T,row.names = 1) %>% select(c(-logFC,-logCPM,-PValue,-FDR))
cr = cor(fc_tmp, method='spearman')
data.TMM.matrix = log2(fc_tmp+1)
data.TMM.matrix = t(scale(t(data.TMM.matrix), scale=F)) # center rows, mean substracted
gene_dist = dist(data.TMM.matrix, method='euclidean')
hc_genes = hclust(gene_dist, method='complete')
hc_samples = hclust(as.dist(1-cr), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);
partition_colors = brewer.pal(length(unique(gene_partition_assignments)), "Accent")
gene_colors = partition_colors[gene_partition_assignments]
quantBrks = quantile(data.TMM.matrix, c(0.03, 0.97))
myheatcol = colorpanel(75,low="#3C439A",mid="#E7E6D5",high="#CB4335")

pdf("../Figures/FigureS1_Hipp.pdf",height=8,width=10)
heatmap.2(data.TMM.matrix, dendrogram='both', Rowv=as.dendrogram(hc_genes),
          col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, cexCol=1, lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1.5,5),lwid=c(1.5,0.2,2.5,1.5), margins=c(10,2), breaks=seq(quantBrks[1], quantBrks[2], length=76), labRow=F,main="Hippocampus", key.title="", key.xlab = "Log2(TMM+1)")
dev.off()

#heatmap_par
fc_tmp <- read.table("../Data/DEGs_ParTemp_edgeR.tsv",head=T,row.names = 1) %>% select(c(-logFC,-logCPM,-PValue,-FDR))
cr = cor(fc_tmp, method='spearman')
data.TMM.matrix = log2(fc_tmp+1)
data.TMM.matrix = t(scale(t(data.TMM.matrix), scale=F)) # center rows, mean substracted
gene_dist = dist(data.TMM.matrix, method='euclidean')
hc_genes = hclust(gene_dist, method='complete')
hc_samples = hclust(as.dist(1-cr), method="complete") # cluster conditions
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);
partition_colors = brewer.pal(length(unique(gene_partition_assignments)), "Accent")
gene_colors = partition_colors[gene_partition_assignments]
quantBrks = quantile(data.TMM.matrix, c(0.03, 0.97))
myheatcol = colorpanel(75,low="#3C439A",mid="#E7E6D5",high="#CB4335")

pdf("../Figures/FigureS1_ParTemp.pdf",height=8,width=10)
heatmap.2(data.TMM.matrix, dendrogram='both', Rowv=as.dendrogram(hc_genes),
          col=myheatcol, RowSideColors=gene_colors, scale="none", density.info="none", trace="none", key=TRUE, cexCol=1, lmat=rbind(c(5,0,4,0),c(3,1,2,0)), lhei=c(1.5,5),lwid=c(1.5,0.2,2.5,1.5), margins=c(10,2), breaks=seq(quantBrks[1], quantBrks[2], length=76), labRow=F,main="Parietal temporal cortex", key.title="", key.xlab = "Log2(TMM+1)")
dev.off()
