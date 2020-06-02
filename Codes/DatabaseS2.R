library(tidyverse)
library(edgeR)
library(foreach)

"+" = function(e1, e2){ #make function just to make it easy to connect strings
  if(is.character(c(e1, e2))){
    paste(e1, e2, sep="")
  }else{
    base::"+"(e1, e2)
  }
}

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

rib <- read.table("../Data/raw_data/Ribosomal_Proteins.tsv",head=F)
names(rib) <- c("GeneID","TranscriptID","Name")

##DE analysis
count <- read.table("../Data/Brain-R7-9-Ami-Hyp-count-s-3_TMM_FPKM.matrix",head=T,row.names = 1)
gl = list(c("w","w","w","d","d","d"), c("w","w","d","w","d","d"),
          c("w","w","d","d","w","d"), c("w","w","d","d","d","w"),
          c("w","d","w","w","d","d"), c("w","d","w","d","w","d"),
          c("w","d","w","d","d","w"), c("w","d","d","w","w","d"),
          c("w","d","d","w","d","w"), c("w","d","d","d","w","w"))
write("Group_pat\tDEGs_Amy\tDEGs_Hyp\tDEGs_Hipp\tDEGs_ParTemp\tDEGs_4regions","../Data/perm_res.tsv")

foreach(i=1:10) %do% { #perform DE analysis for 10 labeling patterns
group <- factor(gl[[i]]) #obtain randomized group label

###Amygdala
count2 <- count %>% select(W2.R8.Ami, W4.R8.Ami, W6.R8.Ami, D1.R8.Ami, D2.R8.Ami, D3.R8.Ami) %>% as.matrix()
d <- DGEList(counts = count2, group = group)
d2 <- d[rowSums(cpm(d) > 1) >=3,]
d3 <- estimateTagwiseDisp(estimateCommonDisp(d2))
result <- exactTest(d3)
table1 <- as.data.frame(topTags(result, n = nrow(result)))
table1_1 <- table1[table1$FDR<0.05 & abs(table1$logFC)>1,]
rownames(table1_1) <- NULL
table1_1$TranscriptID <- rownames(table1[table1$FDR<0.05 & abs(table1$logFC)>1,])

d_cpm <- as.data.frame(cpm(d))
Amy_DEGs <- d_cpm[table1_1$TranscriptID,]
Amy_DEGs$TranscriptID <- rownames(Amy_DEGs)
rownames(Amy_DEGs) <- NULL
wild_col <- c('W2.R8.Ami', 'W4.R8.Ami', 'W6.R8.Ami')
domestic_col <- c('D1.R8.Ami', 'D2.R8.Ami', 'D3.R8.Ami')
if (nrow(Amy_DEGs) > 0){
  Amy_DEGs <- Amy_DEGs %>%
  select(TranscriptID, everything()) %>%
  mutate(SDPerMean_wild = rowSds(as.matrix(.[wild_col])) / rowMeans(.[wild_col]),
         SDPerMean_domestic = rowSds(as.matrix(.[domestic_col])) / rowMeans(.[domestic_col]))
  Amy_rm <- Amy_DEGs %>% rowwise() %>% mutate(check = overlap(c(W2.R8.Ami, W4.R8.Ami, W6.R8.Ami), c(D1.R8.Ami, D2.R8.Ami, D3.R8.Ami))) %>%
    filter(check == 0) %>%
    filter(SDPerMean_wild > 1 | SDPerMean_domestic > 1)
}else{
  Amy_rm <- Amy_DEGs
}

table1_1 <- anti_join(table1_1, Amy_rm, by = "TranscriptID")
nrow(table1_1) #number of DEGs

length(intersect(table1_1$TranscriptID, rib$TranscriptID)) #number of ribosomal proteins
n_deg_amy <- nrow(table1_1) + " (" + length(intersect(table1_1$TranscriptID,rib$TranscriptID)) + ")"

###Hypothalamus
count2 <- count %>% select(W2.R6.Hyp, W4.R6.Hyp, W6.R6.Hyp, D1.R6.Hyp, D2.R6.Hyp, D3.R6.Hyp) %>% as.matrix()
d <- DGEList(counts = count2, group = group)
d2 <- d[rowSums(cpm(d) > 1) >=3,]
d3 <- estimateTagwiseDisp(estimateCommonDisp(d2))
result <- exactTest(d3)
table2 <- as.data.frame(topTags(result, n = nrow(result)))
table2_1 <- table2[table2$FDR<0.05 & abs(table2$logFC)>1,]
rownames(table2_1) <- NULL
table2_1$TranscriptID <- rownames(table2[table2$FDR<0.05 & abs(table2$logFC)>1,])

d_cpm <- as.data.frame(cpm(d))
Hyp_DEGs <- d_cpm[table2_1$TranscriptID,]
Hyp_DEGs$TranscriptID <- rownames(Hyp_DEGs)
rownames(Hyp_DEGs) <- NULL
wild_col <- c('W2.R6.Hyp', 'W4.R6.Hyp', 'W6.R6.Hyp')
domestic_col <- c('D1.R6.Hyp', 'D2.R6.Hyp', 'D3.R6.Hyp')

if(nrow(Hyp_DEGs) > 0){
  Hyp_DEGs <- Hyp_DEGs %>%
    select(TranscriptID, everything()) %>%
    mutate(SDPerMean_wild = rowSds(as.matrix(.[wild_col])) / rowMeans(.[wild_col]),
           SDPerMean_domestic = rowSds(as.matrix(.[domestic_col])) / rowMeans(.[domestic_col]))
  Hyp_rm <- Hyp_DEGs %>% rowwise() %>% mutate(check = overlap(c(W2.R6.Hyp, W4.R6.Hyp, W6.R6.Hyp), c(D1.R6.Hyp, D2.R6.Hyp, D3.R6.Hyp))) %>%
    filter(check == 0) %>%
    filter(SDPerMean_wild > 1 | SDPerMean_domestic > 1)
}else{
  Hyp_rm <- Hyp_DEGs
}

table2_1 <- anti_join(table2_1, Hyp_rm, by = "TranscriptID")
nrow(table2_1) #number of DEGs

length(intersect(table2_1$TranscriptID, rib$TranscriptID)) #number of ribosomal proteins
n_deg_hyp <- nrow(table2_1) + " (" + length(intersect(table2_1$TranscriptID,rib$TranscriptID)) + ")"

###Hippocampus
count2 <- count %>% select(P1_W2.R7, P1_W4.R7, P1_W6.R7, P3_D1.R7, P3_D2.R7, P3_D3.R7) %>% as.matrix()
d <- DGEList(counts = count2, group = group)
d2 <- d[rowSums(cpm(d) > 1) >=3,]
d3 <- estimateTagwiseDisp(estimateCommonDisp(d2))
result <- exactTest(d3)
table3 <- as.data.frame(topTags(result, n = nrow(result)))
table3_1 <- table3[table3$FDR<0.05 & abs(table3$logFC)>1,]
table3_1$TranscriptID <- rownames(table3[table3$FDR<0.05 & abs(table3$logFC)>1,])
rownames(table3_1) <- NULL

d_cpm <- as.data.frame(cpm(d))
Hipp_DEGs <- d_cpm[table3_1$TranscriptID,]
Hipp_DEGs$TranscriptID <- rownames(Hipp_DEGs)
rownames(Hipp_DEGs) <- NULL
wild_col <- c('P1_W2.R7', 'P1_W4.R7', 'P1_W6.R7')
domestic_col <- c('P3_D1.R7', 'P3_D2.R7', 'P3_D3.R7')

if(nrow(Hipp_DEGs) > 0){
  Hipp_DEGs <- Hipp_DEGs %>%
    select(TranscriptID, everything()) %>%
    mutate(SDPerMean_wild = rowSds(as.matrix(.[wild_col])) / rowMeans(.[wild_col]),
           SDPerMean_domestic = rowSds(as.matrix(.[domestic_col])) / rowMeans(.[domestic_col]))
  
  Hipp_rm <- Hipp_DEGs %>% rowwise() %>% mutate(check = overlap(c(P1_W2.R7, P1_W4.R7, P1_W6.R7), c(P3_D1.R7, P3_D2.R7, P3_D3.R7))) %>%
    filter(check == 0) %>%
    filter(SDPerMean_wild > 1 | SDPerMean_domestic > 1)
}else{
  Hipp_rm <- Hipp_DEGs
}

table3_1 <- anti_join(table3_1, Hipp_rm, by = "TranscriptID")

nrow(table3_1) #number of DEGs
length(intersect(table3_1$TranscriptID, rib$TranscriptID)) #number of ribosomal proteins
n_deg_hipp <- nrow(table3_1) + " (" + length(intersect(table3_1$TranscriptID,rib$TranscriptID)) + ")"

###Parietal & temporal cortex
count2 <- count %>% select(P2_W2.R9, P2_W4.R9, P2_W6.R9, P4_D1.R9, P4_D2.R9, P4_D3.R9) %>% as.matrix()
d <- DGEList(counts = count2, group = group)
d2 <- d[rowSums(cpm(d) > 1) >=3,]
d3 <- estimateTagwiseDisp(estimateCommonDisp(d2))
result <- exactTest(d3)
table4 <- as.data.frame(topTags(result, n = nrow(result)))
table4_1 <- table4[table4$FDR<0.05 & abs(table4$logFC)>1,]
rownames(table4_1) <- NULL
table4_1$TranscriptID <- rownames(table4[table4$FDR<0.05 & abs(table4$logFC)>1,])

d_cpm <- as.data.frame(cpm(d))
ParTemp_DEGs <- d_cpm[table4_1$TranscriptID,]
ParTemp_DEGs$TranscriptID <- rownames(ParTemp_DEGs)
rownames(ParTemp_DEGs) <- NULL
wild_col <- c('P2_W2.R9', 'P2_W4.R9', 'P2_W6.R9')
domestic_col <- c('P4_D1.R9', 'P4_D2.R9', 'P4_D3.R9')

if(nrow(ParTemp_DEGs) > 0){
  ParTemp_DEGs <- ParTemp_DEGs %>%
    select(TranscriptID, everything()) %>%
    mutate(SDPerMean_wild = rowSds(as.matrix(.[wild_col])) / rowMeans(.[wild_col]),
           SDPerMean_domestic = rowSds(as.matrix(.[domestic_col])) / rowMeans(.[domestic_col]))
  ParTemp_rm <- ParTemp_DEGs %>% rowwise() %>% mutate(check = overlap(c(P2_W2.R9, P2_W4.R9, P2_W6.R9), c(P4_D1.R9, P4_D2.R9, P4_D3.R9))) %>%
    filter(check == 0) %>%
    filter(SDPerMean_wild > 1 | SDPerMean_domestic > 1)
}else{
  ParTemp_rm <- ParTemp_DEGs
}

table4_1 <- anti_join(table4_1, ParTemp_rm, by = "TranscriptID")
nrow(table4_1) #number of DEGs

length(intersect(table4_1$TranscriptID, rib$TranscriptID)) #number of ribosomal proteins
n_deg_partemp <- nrow(table4_1) + " (" + length(intersect(table4_1$TranscriptID,rib$TranscriptID)) + ")"

##Output the results
table1_2 <- data.frame(table1)
colnames(table1_2) <- c("logFC_Amy","logCPM_Amy","Pvalue_Amy","FDR_Amy")
rownames(table1_2) <- NULL
table1_2$TranscriptID <- rownames(table1)

table2_2 <- data.frame(table2)
colnames(table2_2) <- c("logFC_Hyp","logCPM_Hyp","Pvalue_Hyp","FDR_Hyp")
rownames(table2_2) <- NULL
table2_2$TranscriptID <- rownames(table2)

table3_2 <- data.frame(table3)
colnames(table3_2) <- c("logFC_Hipp","logCPM_Hipp","Pvalue_Hipp","FDR_Hipp")
rownames(table3_2) <- NULL
table3_2$TranscriptID <- rownames(table3)

table4_2 <- data.frame(table4)
colnames(table4_2) <- c("logFC_ParTemp","logCPM_ParTemp","Pvalue_ParTemp","FDR_ParTemp")
rownames(table4_2) <- NULL
table4_2$TranscriptID <- rownames(table4)

edgeR_allres <- table1_2 %>%
  full_join(., table2_2, by = "TranscriptID") %>%
  full_join(., table3_2, by = "TranscriptID") %>%
  full_join(., table4_2, by = "TranscriptID") %>%
  select(TranscriptID, everything())

edgeR_allres %>% filter(abs(logFC_Amy) > 1 & FDR_Amy < 0.05) %>%
  filter(abs(logFC_Hyp) > 1 & FDR_Hyp < 0.05) %>%
  filter(abs(logFC_Hipp) > 1 & FDR_Hipp < 0.05) %>%
  filter(abs(logFC_ParTemp) > 1 & FDR_ParTemp < 0.05) #%>%
#  write.table(., file = "../Data/DEanalysis_perm/DEGs_" + i + "_4regions_edgeR.tsv", quote=F, col.names = T, row.names = F, sep = "\t")

consis <- inner_join(table1_1, table2_1, by = "TranscriptID") %>%
  inner_join(table3_1, by = "TranscriptID") %>%
  inner_join(table4_1, by = "TranscriptID")
row.names(edgeR_allres) <- edgeR_allres$TranscriptID

degs_4regions <- edgeR_allres[consis$TranscriptID, ]

n_deg_4reg <- nrow(degs_4regions) + " (" + length(intersect(degs_4regions$TranscriptID,rib$TranscriptID)) + ")"
write(as.character(gl[i]) + "\t" + n_deg_amy + "\t" + n_deg_hyp + "\t" + n_deg_hipp + "\t" + n_deg_partemp + "\t" + n_deg_4reg, "../Data/perm_res.tsv", append=TRUE)

}

