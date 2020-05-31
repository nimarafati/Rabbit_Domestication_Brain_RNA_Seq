## Figure S5
library(tidyverse)
library(ggfortify)
library(DESeq2)
library(edgeR)

TMMcount <- read.table("../Data/Brain-R7-9-Ami-Hyp-count-s-3_TMM_FPKM.matrix",head=T,row.names = 1)
group <- rep(gl(2, 6, labels=c("Wild","Domestic"),length=12),2)
region <- list(gl(2, 3, labels=c("Hippocampus","Parietal / temporal cortex"),length=12),gl(2, 3, labels=c("Amygdala","Hypothalamus"),length=12)) %>% unlist()

sampleinfo <- data.frame(name=names(count),group=group,region=region) %>% 
  transform(region = factor(region, levels = c("Amygdala","Hypothalamus","Hippocampus","Parietal / temporal cortex")))

dgeObj <- DGEList(TMMcount) 
dgeObj <- dgeObj[rowSums(cpm(dgeObj) > 1) >= 3,]

# run PCA
logcounts <- log2(dgeObj_N$counts + 1)
pcDat <- prcomp(t(logcounts))

# plot PCA
g<-autoplot(pcDat, data = sampleinfo, col="group",fill="group", shape="region", size=5) +
  scale_shape_manual(values=c(18,17,19,15)) + 
  scale_color_manual(values=c("#945a7f","#b8cc97")) +
  scale_fill_manual(values=c("#945a7f","#b8cc97")) +
  theme_bw() + 
  theme(legend.title = element_blank())

ggsave("../Figures/FigureS5.pdf",g,width=6,height=4)
