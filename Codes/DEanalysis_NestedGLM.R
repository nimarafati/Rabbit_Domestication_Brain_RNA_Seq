library(tidyverse)
library(edgeR)

"+" = function(e1, e2){ #make function just to make it easy to connect strings
  if(is.character(c(e1, e2))){
    paste(e1, e2, sep="")
  }else{
    base::"+"(e1, e2)
  }
}

TMMcount <- read.table("../Data/Brain-R7-9-Ami-Hyp-count-s-3_TMM_FPKM.matrix",head=T,row.names = 1) #load normalized count data
data_dop <- read.table("../Data/raw_data/Dopamine_Ciliary_genes.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% filter(Feature == "Dopamine")

# exactTest
## Amygdala
count2 <- count %>% select(W2.R8.Ami, W4.R8.Ami, W6.R8.Ami, D1.R8.Ami, D2.R8.Ami, D3.R8.Ami) %>% as.matrix()
d <- DGEList(counts = count2, group = group)
d2 <- d[rowSums(cpm(d) > 1) >=3,]
d3 <- estimateTagwiseDisp(estimateCommonDisp(d2))
result <- exactTest(d3)
table1 <- as.data.frame(topTags(result, n = nrow(result)))
print(nrow(table1[table1$FDR<0.05 & abs(table1$logFC)>1,])) #number of DEGs

table1_1 <- data.frame(table1[table1$FDR<0.05 & abs(table1$logFC)>1,])
rownames(table1_1) <- NULL
table1_1$TranscriptID <- rownames(table1[table1$FDR<0.05 & abs(table1$logFC)>1,])

## Hypothalamus
count2 <- count %>% select(W2.R6.Hyp, W4.R6.Hyp, W6.R6.Hyp, D1.R6.Hyp, D2.R6.Hyp, D3.R6.Hyp) %>% as.matrix()
d <- DGEList(counts = count2, group = group)
d2 <- d[rowSums(cpm(d) > 1) >=3,]
d3 <- estimateTagwiseDisp(estimateCommonDisp(d2))
result <- exactTest(d3)
table2 <- as.data.frame(topTags(result, n = nrow(result)))
print(nrow(table2[table2$FDR<0.05 & abs(table2$logFC)>1,])) #number of DEGs

table2_1 <- data.frame(table2[table2$FDR<0.05 & abs(table2$logFC)>1,])
rownames(table2_1) <- NULL
table2_1$TranscriptID <- rownames(table2[table2$FDR<0.05 & abs(table2$logFC)>1,])

## Hippocampus
count2 <- count %>% select(P1_W2.R7, P1_W4.R7, P1_W6.R7, P3_D1.R7, P3_D2.R7, P3_D3.R7) %>% as.matrix()
d <- DGEList(counts = count2, group = group)
d2 <- d[rowSums(cpm(d) > 1) >=3,]
d3 <- estimateTagwiseDisp(estimateCommonDisp(d2))
result <- exactTest(d3)
table3 <- as.data.frame(topTags(result, n = nrow(result)))
print(nrow(table3[table3$FDR<0.05 & abs(table3$logFC)>1,])) #number of DEGs

table3_1 <- data.frame(table3[table3$FDR<0.05 & abs(table3$logFC)>1,])
rownames(table3_1) <- NULL
table3_1$TranscriptID <- rownames(table3[table3$FDR<0.05 & abs(table3$logFC)>1,])

## Parietal & temporal cortex
count2 <- count %>% select(P2_W2.R9, P2_W4.R9, P2_W6.R9, P4_D1.R9, P4_D2.R9, P4_D3.R9) %>% as.matrix()
d <- DGEList(counts = count2, group = group)
d2 <- d[rowSums(cpm(d) > 1) >=3,]
d3 <- estimateTagwiseDisp(estimateCommonDisp(d2))
result <- exactTest(d3)
table4 <- as.data.frame(topTags(result, n = nrow(result)))
print(nrow(table4[table4$FDR<0.05 & abs(table4$logFC)>1,])) #number of DEGs

table4_1 <- data.frame(table4[table4$FDR<0.05 & abs(table4$logFC)>1,])
rownames(table4_1) <- NULL
table4_1$TranscriptID <- rownames(table4[table4$FDR<0.05 & abs(table4$logFC)>1,])


# nested GLM LRT
ind <- list(rep(unlist(c("w"+1:3*2)),2),rep(unlist(c("d"+1:3)),2)) %>% unlist() %>% rep(.,2) %>% as.factor()
reg <- list(gl(2, 3, labels=c("Hipp","ParTemp"),length=12),gl(2, 3, labels=c("Hyp","Amy"),length=12)) %>% unlist()
dw <- str_split_fixed(ind, "", n=2)[,1] %>% as.factor()

design <- model.matrix(~ 0 + reg:dw)
design
colnames(design) <- c("Hipp.d", "ParTemp.d", "Hyp.d", "Amy.d", "Hipp.w", "ParTemp.w", "Hyp.w", "Amy.w")

d <- DGEList(counts = TMMcount, group = reg:dw)
d2 <- d[rowSums(cpm(d) > 1) >=3,]
d3 <- estimateGLMCommonDisp(d2, design)
d3 <- estimateGLMTrendedDisp(d3, design)
d3 <- estimateGLMTagwiseDisp(d3, design)

# GLM with likelihood ratio test
fit <- glmFit(d3, design, robust=TRUE)

con <- makeContrasts(Amy.w - Amy.d, levels=design)
qlf <- glmLRT(fit, contrast=con)
nested1 <- topTags(qlf, n = nrow(tr)) %>% data.frame() %>%
  mutate(TranscriptID = rownames(.)) %>%
  select(TranscriptID, everything()) %>%
  rename(FDR_Amy_NestedGLM = "FDR", logFC_Amy_NestedGLM = "logFC", logCPM_Amy_NestedGLM = "logCPM", PValue_Amy_NestedGLM = "PValue", LR_Amy_NestedGLM = "LR")

con <- makeContrasts(Hyp.w - Hyp.d, levels=design)
qlf <- glmLRT(fit, contrast=con)
nested2 <- topTags(qlf, n = nrow(tr)) %>% data.frame() %>%
  mutate(TranscriptID = rownames(.)) %>%
  select(TranscriptID, everything()) %>%
  rename(FDR_Hyp_NestedGLM = "FDR", logFC_Hyp_NestedGLM = "logFC", logCPM_Hyp_NestedGLM = "logCPM", PValue_Hyp_NestedGLM = "PValue", LR_Hyp_NestedGLM = "LR")

con <- makeContrasts(Hipp.w - Hipp.d, levels=design)
qlf <- glmLRT(fit, contrast=con)
nested3 <- topTags(qlf, n = nrow(tr)) %>% data.frame() %>%
  mutate(TranscriptID = rownames(.)) %>%
  select(TranscriptID, everything()) %>%
  rename(FDR_Hipp_NestedGLM = "FDR", logFC_Hipp_NestedGLM = "logFC", logCPM_Hipp_NestedGLM = "logCPM", PValue_Hipp_NestedGLM = "PValue", LR_Hipp_NestedGLM = "LR")

con <- makeContrasts(ParTemp.w - ParTemp.d, levels=design)
qlf <- glmLRT(fit, contrast=con)
nested4 <- topTags(qlf, n = nrow(tr)) %>% data.frame() %>%
  mutate(TranscriptID = rownames(.)) %>%
  select(TranscriptID, everything()) %>%
  rename(FDR_ParTemp_NestedGLM = "FDR", logFC_ParTemp_NestedGLM = "logFC", logCPM_ParTemp_NestedGLM = "logCPM", PValue_ParTemp_NestedGLM = "PValue", LR_ParTemp_NestedGLM = "LR")

edgeR_NestedGLM_allres <- nested1 %>%
  full_join(., nested2, by = "TranscriptID") %>%
  full_join(., nested3, by = "TranscriptID") %>%
  full_join(., nested4, by = "TranscriptID") %>%
  select(TranscriptID, everything())

write.table(edgeR_NestedGLM_allres, file = "../Data/Allgenes_edgeR_NestedGLM.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

#All
data_rib <- read.table("../Data/raw_data/Ribosomal_Proteins.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% mutate(Feature = "Ribosomal")
data_dop <- read.table("../Data/raw_data/Dopamine_Ciliary_genes.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% filter(Feature == "Dopamine")
data_cil <- read.table("../Data/raw_data/Dopamine_Ciliary_genes.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% filter(Feature == "Ciliary")
data_4reg <- read.table("../Data/DEGs_4regions_edgeR_name.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% mutate(Feature = "4regions")
edgeR_allres <- read.table("../Data/Allgenes_edgeR.tsv",sep="\t",head=T,stringsAsFactors=FALSE)

data_rib %>% 
  inner_join(edgeR_allres, by = "TranscriptID") %>%
  inner_join(edgeR_NestedGLM_allres, by = "TranscriptID") %>%
  select(TranscriptID, Name, Feature, FDR_Amy, FDR_Amy_NestedGLM, FDR_Hyp, FDR_Hyp_NestedGLM, FDR_Hipp, FDR_Hipp_NestedGLM, FDR_ParTemp, FDR_ParTemp_NestedGLM) %>%
  write.table(., file = "../Data/edgeR_Ribosomal_comp_NestedGLM.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

data_cil %>% 
  inner_join(edgeR_allres, by = "TranscriptID") %>% 
  inner_join(edgeR_NestedGLM_allres, by = "TranscriptID") %>%
  select(TranscriptID, Name, Feature, FDR_Amy, FDR_Amy_NestedGLM, FDR_Hyp, FDR_Hyp_NestedGLM, FDR_Hipp, FDR_Hipp_NestedGLM, FDR_ParTemp, FDR_ParTemp_NestedGLM) %>%
  write.table(., file = "../Data/edgeR_Ciliary_comp_NestedGLM.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

data_dop %>% 
  inner_join(edgeR_allres, by = "TranscriptID") %>%
  inner_join(edgeR_NestedGLM_allres, by = "TranscriptID") %>%
  select(TranscriptID, Name, Feature, FDR_Amy, FDR_Amy_NestedGLM, FDR_Hyp, FDR_Hyp_NestedGLM, FDR_Hipp, FDR_Hipp_NestedGLM, FDR_ParTemp, FDR_ParTemp_NestedGLM) %>%
  write.table(., file = "../Data/edgeR_Dopamine_comp_NestedGLM.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

data_4reg %>% 
  inner_join(edgeR_allres, by = "TranscriptID") %>%
  inner_join(edgeR_NestedGLM_allres, by = "TranscriptID") %>%
  select(TranscriptID, Name, Feature, FDR_Amy, FDR_Amy_NestedGLM, FDR_Hyp, FDR_Hyp_NestedGLM, FDR_Hipp, FDR_Hipp_NestedGLM, FDR_ParTemp, FDR_ParTemp_NestedGLM) %>%
  write.table(., file = "../Data/edgeR_4reg_comp_NestedGLM.tsv", quote=F, col.names = T, row.names = F, sep = "\t")
