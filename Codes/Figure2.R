#install.packages("GGally")
library(GGally)
library(network)
library(sna)
library(tidyverse)
library(ggrepel)
library(cowplot)
#install.packages("patchwork")
#library(patchwork)
#library(gridExtra)
#install.packages("gghighlight")
#library(gghighlight)

data_dAF <- read.table("../Data/dAF_DEGs_results.tsv",sep="\t",head=T,stringsAsFactors=FALSE)
data_dop <- read.table("../Data/raw_data/Dopamine_Ciliary_genes.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% filter(Feature == "Dopamine")
data_rib <- read.table("../Data/raw_data/Ribosomal_Proteins.tsv",head=T,stringsAsFactors=FALSE) %>% mutate(Feature1 = "Ribosomal")
mapped_amy <- read.delim("../Data/STRING_mapping_Amy.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% dplyr::select(TranscriptID = queryItem, everything())
data_amy <- read.table("../Data/DEGs_Amy_edgeR.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% 
  left_join(data_dAF, by="TranscriptID") %>%
  left_join(data_dop, by="TranscriptID") %>%
  left_join(data_rib, by="TranscriptID") %>%
  select(c("TranscriptID", "Genename", "updown100kb_dAF0.9", "Feature", "Feature1")) %>%
  mutate(sel = if_else(updown100kb_dAF0.9 > 2, "#a22041", "black"),
         col = case_when(Feature == "Dopamine" ~ "#a8bf93", #915c8b 839b5c
                         Feature1 == "Ribosomal" ~ "#595857", #ffdb4f
                         Feature == "NA" & Feature1 == "NA" ~ "#eae5e3"),
         font = case_when(Feature == "Dopamine" ~ "bold",
#                          Feature1 == "Ribosomal" ~ "bold",
                          Feature == "NA" & Feature1 == "NA" ~ "plain"),
         lab.size = case_when(Feature == "Dopamine" ~ 2.5,
#                              Feature1 == "Ribosomal" ~ 2.5,
                              Feature == "NA" & Feature1 == "NA" ~ 1.5)
  ) %>%
  mutate_all(funs(ifelse(is.na(.),1.5,.))) %>%
  mutate(col = if_else(col == 1.5, "#eae5e3", col),
         font = if_else(font == 1.5, "plain", font)) %>%
  right_join(mapped_amy, by="TranscriptID")

rownames(data_amy) <- data_amy$preferredName
amy_net <- read.delim("../Data/STRING_res_Amy.tsv", head=T, stringsAsFactors=FALSE)
amy.net <- network(amy_net[,1:2], directed = TRUE)
amy.net %v% "col" <- data_amy[network.vertex.names(amy.net),]$col
amy.net %v% "lab.col" <- data_amy[network.vertex.names(amy.net),]$sel
amy.net %v% "lab.size" <- data_amy[network.vertex.names(amy.net),]$lab.size
network::set.edge.attribute(amy.net, "weights", amy_net$combined_score*1.5)

data_cil <- read.table("../Data/raw_data/Dopamine_Ciliary_genes.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% filter(Feature == "Ciliary")
mapped_hipp <- read.delim("../Data/STRING_mapping_Hipp.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% dplyr::select(TranscriptID = queryItem, everything())
data_hipp <- read.table("../Data/DEGs_Hipp_edgeR.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% 
  left_join(data_dAF, by="TranscriptID") %>%
  left_join(data_cil, by="TranscriptID") %>%
  left_join(data_rib, by="TranscriptID") %>%
  select(c("TranscriptID", "Genename", "updown100kb_dAF0.9", "Feature", "Feature1")) %>%
  mutate(sel = if_else(updown100kb_dAF0.9 > 2, "#a22041", "black"),
         col = case_when(Feature == "Ciliary" ~ "#c4a3bf", #915c8b 
                         Feature1 == "Ribosomal" ~ "#595857", #ffdb4f
                         Feature == "NA" & Feature1 == "NA" ~ "#eae5e3"),
         font = case_when(Feature == "Ciliary" ~ "bold",
#                          Feature1 == "Ribosomal" ~ "bold",
                          Feature == "NA" & Feature1 == "NA" ~ "plain"),
         lab.size = case_when(Feature == "Ciliary" ~ 2.5,
#                              Feature1 == "Ribosomal" ~ 2.5,
                              Feature == "NA" & Feature1 == "NA" ~ 1.5)
  ) %>%
  mutate_all(funs(ifelse(is.na(.),1.5,.))) %>%
  mutate(col = if_else(col == 1.5 | Genename == "ROPN1L" | Genename == "RAB21" | Genename == "TTLL8", "#eae5e3", col),
         font = if_else(font == 1.5, "plain", font)) %>%
  right_join(mapped_hipp, by="TranscriptID")


rownames(data_hipp) <- data_hipp$preferredName
hipp_net <- read.delim("../Data/STRING_res_Hipp.tsv", head=T, stringsAsFactors=FALSE)
hipp.net <- network(hipp_net[,1:2], directed = TRUE)
hipp.net %v% "col" <- data_hipp[network.vertex.names(hipp.net),]$col
hipp.net %v% "lab.col" <- data_hipp[network.vertex.names(hipp.net),]$sel
hipp.net %v% "lab.size" <- data_hipp[network.vertex.names(hipp.net),]$lab.size
network::set.edge.attribute(hipp.net, "weights", hipp_net$combined_score*1.5)

#pdf("../Figures/Figure2.pdf", w=10, h=6)
set.seed(9)
g1 <- ggnet2(amy.net, color = "col", label.color = "lab.col", label.size = "lab.size", 
             edge.color = "#17184b", edge.size = "weights", edge.alpha = .4, node.size = 5, 
             node.alpha = .9, node.label = FALSE)  + 
  geom_text_repel(size = data_amy[network.vertex.names(amy.net),]$lab.size,
                  color = data_amy[network.vertex.names(amy.net),]$sel,
                  label = data_amy[network.vertex.names(amy.net),]$preferredName,
#                  fontface = data_amy[network.vertex.names(amy.net),]$font, 
#                  box.padding = .05, segment.size = .1,
                  nudge_x       = -.05,
                  nudge_y       = -.0001,
                  max.iter      = 10000, 
                  segment.size  = .2,
                  segment.color = "grey50",
                  segment.alpha = .8,
                  box.padding   = .05,
                  direction     = "y",
                  hjust         = .5)

g1
set.seed(22)
g2 <- ggnet2(hipp.net, color = "col", label.color = "lab.col", label.size = "lab.size", 
             edge.color = "#17184b", edge.size = "weights", edge.alpha = .4, node.size = 5, 
             node.alpha = .9, node.label = FALSE)  + 
  geom_text_repel(size = data_hipp[network.vertex.names(hipp.net),]$lab.size, 
                  color = data_hipp[network.vertex.names(hipp.net),]$sel, 
                  label = data_hipp[network.vertex.names(hipp.net),]$preferredName,
#                  fontface = data_hipp[network.vertex.names(hipp.net),]$font, 
#                  box.padding = .05, segment.size = .1,
                  
                  nudge_x       = -.05,
                  nudge_y       = -.0001,
                  max.iter      = 10000, 
                  segment.size  = .2,
                  segment.color = "grey50",
                  segment.alpha = .8,
                  box.padding   = .05,
                  direction     = "y",
                  hjust         = .5)
g2
#cowplot::plot_grid(g1, g2, scale=0.95, ncol = 2, labels = "AUTO")
#dev.off()

data_amy <- read.table("../Data/DEGs_Amy_edgeR.tsv",sep="\t",head=T,stringsAsFactors=FALSE)
data_hipp <- read.table("../Data/DEGs_Hipp_edgeR.tsv",sep="\t",head=T,stringsAsFactors=FALSE)
data_cil <- data_cil %>%
  inner_join(., data_hipp, by = "TranscriptID") %>%
  select(TranscriptID, Name, logFC, FDR)
data_tmp <- data_dop %>%
  inner_join(., data_amy, by = "TranscriptID") %>%
  select(TranscriptID, Name, logFC, FDR) %>%
  rbind(., c("Dummy", "NA", 1, 0), data_cil)

data_tmp = data_tmp %>% transform(Name=factor(Name,levels=c("ADORA2A","DRD1","DRD2","GPR88","LAMP5","PDE1B","PDE10A","PENK","PPP1R1B","RASD2","RGS9","TAC1","NA","C1orf87","CCDC37","CCDC108","CCNO","DNAH1","DNAH7","DNAH9","DNAH10","DNAH11","DNAH12","DNAI2","DYNLRB2","HYDIN","KIF19","PPP1R32","RAB21","ROPN1L","RSPH9","TTLL8","TTLL13P","ZMYND10"))) %>%
  mutate(logFC = as.numeric(logFC),
         FDR = as.numeric(FDR))
g3 <- ggplot(data_tmp,aes(x=Name,y=logFC)) +
  geom_point(aes(size=-log10(FDR),col=logFC),alpha=.9) +
  #  geom_hline(yintercept=0, col="black", linetype="dotted") + 
  scale_colour_gradient2(low = "#839b5c", high = "#74325c", mid = "#E7E6D5") +
  ylab(expression(paste("Log"[2],"(Fold Change)"))) +
  labs(size=expression(paste("-Log"[10],"(FDR)")), col=expression(paste("Log"[2],"FC")), tags = "C") +
  theme_bw(base_size = 11, base_family = "") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(face="italic", angle = 60, hjust = 1), plot.tag = element_text(hjust = 0.2, size = 15, face = "bold"), plot.margin = margin(10, 0, 0, 0))
#g3

pdf("../Figures/Figure2.pdf", w=7, h=7)
g12 <- plot_grid(g1, g2, labels=c("A","B"), label_size = 15, ncol=2, align="v") 
#g12
cowplot::plot_grid(g12, g3, rel_heights=c(3, 2.5), scale=0.95, ncol = 1)
dev.off()
