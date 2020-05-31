## Figure 3
library(tidyverse)

data_dop <- read.table("../Data/raw_data/Dopamine_Ciliary_genes.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% filter(Feature == "Dopamine")
data_cil <- read.table("../Data/raw_data/Dopamine_Ciliary_genes.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% filter(Feature == "Ciliary")
data_amy <- read.table("../Data/DEGs_Amy_edgeR.tsv",sep="\t",head=T,stringsAsFactors=FALSE)
data_hipp <- read.table("../Data/DEGs_Hipp_edgeR.tsv",sep="\t",head=T,stringsAsFactors=FALSE)
data_cil <- data_cil %>%
  left_join(., data_hipp, by = "TranscriptID") %>%
  select(TranscriptID, Name, logFC, FDR)
data_tmp <- data_dop %>%
  left_join(., data_amy, by = "TranscriptID") %>%
  select(TranscriptID, Name, logFC, FDR) %>%
  rbind(., c("Dummy", "NA", 1, 0), data_cil)

data_tmp = data_tmp %>% transform(Name=factor(Name,levels=rev(c("ADORA2A","DRD1","DRD2","GPR88","LAMP5","PDE1B","PDE10A","PENK","PPP1R1B","RASD2","RGS9","TAC1","NA","C1orf87","CCDC37","CCDC108","CCNO","DNAH1","DNAH7","DNAH9","DNAH10","DNAH11","DNAH12","DNAI2","DYNLRB2","HYDIN","KIF19","PPP1R32","RAB21","ROPN1L","RSPH9","TTLL8","TTLL13P","ZMYND10")))) %>%
  mutate(logFC = as.numeric(logFC),
         FDR = as.numeric(FDR))
g <- ggplot(data_tmp,aes(x=logFC,y=Name)) +
  geom_point(aes(size=-log10(FDR),col=logFC),alpha=.9) +
  geom_hline(yintercept=0,col="black",linetype="dotted") + 
  scale_colour_gradient2(low = "#839b5c", high = "#74325c", mid = "#E7E6D5") +
  xlab(expression(paste("Log"[2],"(Fold Change)"))) +
  labs(size=expression(paste("-Log"[10],"(FDR)")),col=expression(paste("Log"[2],"FC"))) +
  theme_bw(base_size = 11, base_family = "") +
  theme(axis.title.y=element_blank(), axis.text.y = element_text(face="italic"))

ggsave("../Figures/Figure3.pdf" , g, width=5, height=5)
