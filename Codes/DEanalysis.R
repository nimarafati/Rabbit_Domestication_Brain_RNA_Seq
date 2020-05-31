library(tidyverse)
library(edgeR)
count <- read.table("../Data/Brain-R7-9-Ami-Hyp-count-s-3_TMM_FPKM.matrix",head=T,row.names = 1)  #load normalized count data
group <- factor(c("w", "w", "w", "d", "d", "d"))

overlap = function(l1, l2){ #make function to check overalap in the distribution of the CPM between two group
  x = 0
  for (i in 1:length(l1)) {
    for (j in 1:length(l2)) {
      if (l1[i] < l2[j]){
        x <- x + 1
      }
    }
  }
  if(x == 0 || x == length(l1)*length(l2)) return(1)
  else return(0)
}

##DE analysis with exactTest
###Amygdala
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
nrow(table1_1)

d_cpm <- cpm(d)
Amy_DEGs <- as.data.frame(d_cpm[table1_1$TranscriptID,])
Amy_DEGs$TranscriptID <- rownames(Amy_DEGs)
rownames(Amy_DEGs) <- NULL
wild_col <- c('W2.R8.Ami', 'W4.R8.Ami', 'W6.R8.Ami')
domestic_col <- c('D1.R8.Ami', 'D2.R8.Ami', 'D3.R8.Ami')
Amy_DEGs <- Amy_DEGs %>%
  select(TranscriptID, everything()) %>%
  mutate(SDPerMean_wild = rowSds(as.matrix(.[wild_col])) / rowMeans(.[wild_col]),
         SDPerMean_domestic = rowSds(as.matrix(.[domestic_col])) / rowMeans(.[domestic_col]))

Amy_rm <- Amy_DEGs %>% rowwise() %>% mutate(check = overlap(c(W2.R8.Ami, W4.R8.Ami, W6.R8.Ami), c(D1.R8.Ami, D2.R8.Ami, D3.R8.Ami))) %>%
  filter(check == 0) %>%
  filter(SDPerMean_wild > 1 | SDPerMean_domestic > 1) #%>% nrow()

table1_1 <- anti_join(table1_1, Amy_rm, by = "TranscriptID")

count3 <- as.data.frame(count2) %>%
  mutate(TranscriptID = rownames(count2)) %>%
  inner_join(., table1_1, by = "TranscriptID") %>%
  select(TranscriptID, everything()) %>%
  write.table(., file = "../Data/DEGs_Amy_edgeR.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

###Hypothalamus
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
nrow(table2_1)

d_cpm <- cpm(d)
Hyp_DEGs <- as.data.frame(d_cpm[table2_1$TranscriptID,])
Hyp_DEGs$TranscriptID <- rownames(Hyp_DEGs)
rownames(Hyp_DEGs) <- NULL
wild_col <- c('W2.R6.Hyp', 'W4.R6.Hyp', 'W6.R6.Hyp')
domestic_col <- c('D1.R6.Hyp', 'D2.R6.Hyp', 'D3.R6.Hyp')
Hyp_DEGs <- Hyp_DEGs %>%
  select(TranscriptID, everything()) %>%
  mutate(SDPerMean_wild = rowSds(as.matrix(.[wild_col])) / rowMeans(.[wild_col]),
         SDPerMean_domestic = rowSds(as.matrix(.[domestic_col])) / rowMeans(.[domestic_col]))

Hyp_rm <- Hyp_DEGs %>% rowwise() %>% mutate(check = overlap(c(W2.R6.Hyp, W4.R6.Hyp, W6.R6.Hyp), c(D1.R6.Hyp, D2.R6.Hyp, D3.R6.Hyp))) %>%
  filter(check == 0) %>%
  filter(SDPerMean_wild > 1 | SDPerMean_domestic > 1) #%>% nrow()

table2_1 <- anti_join(table2_1, Hyp_rm, by = "TranscriptID")

count3 <- as.data.frame(count2) %>%
  mutate(TranscriptID = rownames(count2)) %>%
  inner_join(., table2_1, by = "TranscriptID") %>%
  select(TranscriptID, everything()) %>%
  write.table(., file = "../Data/DEGs_Hyp_edgeR.tsv", quote=F, col.names = T, row.names = F, sep = "\t")


###Hippocampus
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
nrow(table3_1)

d_cpm <- cpm(d)
Hipp_DEGs <- as.data.frame(d_cpm[table3_1$TranscriptID,])
Hipp_DEGs$TranscriptID <- rownames(Hipp_DEGs)
rownames(Hipp_DEGs) <- NULL
wild_col <- c('P1_W2.R7', 'P1_W4.R7', 'P1_W6.R7')
domestic_col <- c('P3_D1.R7', 'P3_D2.R7', 'P3_D3.R7')
Hipp_DEGs <- Hipp_DEGs %>%
  select(TranscriptID, everything()) %>%
  mutate(SDPerMean_wild = rowSds(as.matrix(.[wild_col])) / rowMeans(.[wild_col]),
         SDPerMean_domestic = rowSds(as.matrix(.[domestic_col])) / rowMeans(.[domestic_col]))

Hipp_rm <- Hipp_DEGs %>% rowwise() %>% mutate(check = overlap(c(P1_W2.R7, P1_W4.R7, P1_W6.R7), c(P3_D1.R7, P3_D2.R7, P3_D3.R7))) %>%
  filter(check == 0) %>%
  filter(SDPerMean_wild > 1 | SDPerMean_domestic > 1) #%>% nrow()

table3_1 <- anti_join(table3_1, Hipp_rm, by = "TranscriptID")

count3 <- as.data.frame(count2) %>%
  mutate(TranscriptID = rownames(count2)) %>%
  inner_join(., table3_1, by = "TranscriptID") %>%
  select(TranscriptID, everything()) %>%
  write.table(., file = "../Data/DEGs_Hipp_edgeR.tsv", quote=F, col.names = T, row.names = F, sep = "\t")


###Parietal & temporal cortex
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
nrow(table4_1)

d_cpm <- cpm(d)
ParTemp_DEGs <- as.data.frame(d_cpm[table4_1$TranscriptID,])
ParTemp_DEGs$TranscriptID <- rownames(ParTemp_DEGs)
rownames(ParTemp_DEGs) <- NULL
wild_col <- c('P2_W2.R9', 'P2_W4.R9', 'P2_W6.R9')
domestic_col <- c('P4_D1.R9', 'P4_D2.R9', 'P4_D3.R9')
ParTemp_DEGs <- ParTemp_DEGs %>%
  select(TranscriptID, everything()) %>%
  mutate(SDPerMean_wild = rowSds(as.matrix(.[wild_col])) / rowMeans(.[wild_col]),
         SDPerMean_domestic = rowSds(as.matrix(.[domestic_col])) / rowMeans(.[domestic_col]))

ParTemp_rm <- ParTemp_DEGs %>% rowwise() %>% mutate(check = overlap(c(P2_W2.R9, P2_W4.R9, P2_W6.R9), c(P4_D1.R9, P4_D2.R9, P4_D3.R9))) %>%
  filter(check == 0) %>%
  filter(SDPerMean_wild > 1 | SDPerMean_domestic > 1) #%>% nrow()

table4_1 <- anti_join(table4_1, ParTemp_rm, by = "TranscriptID")

count3 <- as.data.frame(count2) %>%
  mutate(TranscriptID = rownames(count2)) %>%
  select(TranscriptID, everything()) %>%
  inner_join(., table4_1, by = "TranscriptID") %>%
  select(TranscriptID, everything()) %>%
  write.table(., file = "../Data/DEGs_ParTemp_edgeR.tsv", quote=F, col.names = T, row.names = F, sep = "\t")


##Output the results
table1_2 <- data.frame(table1)
colnames(table1_2) <- c("logFC_Amy","logCPM_Amy","Pvalue_Amy","FDR_Amy")
rownames(table1_2) <- NULL
table1_2$TranscriptID <- rownames(table1)
table1_2 %>%
  select(TranscriptID, everything()) %>%
  write.table(., file = "../Data/Allgenes_Amy_edgeR.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

table2_2 <- data.frame(table2)
colnames(table2_2) <- c("logFC_Hyp","logCPM_Hyp","Pvalue_Hyp","FDR_Hyp")
rownames(table2_2) <- NULL
table2_2$TranscriptID <- rownames(table2)
table2_2 %>%
  select(TranscriptID, everything()) %>%
  write.table(., file = "../Data/Allgenes_Hyp_edgeR.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

table3_2 <- data.frame(table3)
colnames(table3_2) <- c("logFC_Hipp","logCPM_Hipp","Pvalue_Hipp","FDR_Hipp")
rownames(table3_2) <- NULL
table3_2$TranscriptID <- rownames(table3)
table3_2 %>%
  select(TranscriptID, everything()) %>%
  write.table(., file = "../Data/Allgenes_Hipp_edgeR.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

table4_2 <- data.frame(table4)
colnames(table4_2) <- c("logFC_ParTemp","logCPM_ParTemp","Pvalue_ParTemp","FDR_ParTemp")
rownames(table4_2) <- NULL
table4_2$TranscriptID <- rownames(table4)
table4_2 %>%
  select(TranscriptID, everything()) %>%
  write.table(., file = "../Data/Allgenes_ParTemp_edgeR.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

edgeR_allres <- table1_2 %>%
  full_join(., table2_2, by = "TranscriptID") %>%
  full_join(., table3_2, by = "TranscriptID") %>%
  full_join(., table4_2, by = "TranscriptID") %>%
  select(TranscriptID, everything())

write.table(edgeR_allres, file = "../Data/Allgenes_edgeR.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

consis <- inner_join(table1_1, table2_1, by = "TranscriptID") %>%
  inner_join(table3_1, by = "TranscriptID") %>%
  inner_join(table4_1, by = "TranscriptID")
row.names(edgeR_allres) <- edgeR_allres$TranscriptID

write.table(edgeR_allres[consis$TranscriptID, ], file = "../Data/DEGs_4regions_edgeR.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

