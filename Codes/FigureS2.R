## Figure S2
library(tidyverse)
library(cowplot)

data_amy<-read.table("../Data/DEGs_Amy_edgeR.tsv",head=T,stringsAsFactors=FALSE)
data_dop <- read.table("../Data/raw_data/Dopamine_Ciliary_genes.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% filter(Feature == "Dopamine")
data_amy <- data_amy %>% inner_join(., data_dop, by="TranscriptID") %>%
  pivot_longer(col=c("W2.R8.Ami", "W4.R8.Ami", "W6.R8.Ami", "D1.R8.Ami", "D2.R8.Ami", "D3.R8.Ami"), names_to="sample", values_to="tmm_fpkm") %>%
  mutate(
    class=case_when(
      sample == "W2.R8.Ami" ~ "Wild",
      sample == "W4.R8.Ami" ~ "Wild",
      sample == "W6.R8.Ami" ~ "Wild",
      sample == "D1.R8.Ami" ~ "Domestic",
      sample == "D2.R8.Ami" ~ "Domestic",
      sample == "D3.R8.Ami" ~ "Domestic"
    )
  )

g1 <- ggplot(data_amy,aes(x=class, y=log2(tmm_fpkm+1), color=class)) + 
  facet_wrap(vars(Name), ncol = 6, strip.position = "top") +
  geom_point() + 
  scale_color_manual(values=c("#b8cc97","#945a7f")) +
  ylab(expression(paste("Log"[2],"(TMM+1)"))) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.title = element_blank(),
        strip.placement = "outside", strip.background = element_blank(), strip.text = element_text(face="italic"))


data_hipp<-read.table("../Data/DEGs_Hipp_edgeR.tsv",head=T,stringsAsFactors=FALSE)
data_cil <- read.table("../Data/raw_data/Dopamine_Ciliary_genes.tsv",sep="\t",head=T,stringsAsFactors=FALSE) %>% filter(Feature == "Ciliary")
data_hipp <- data_hipp %>% inner_join(., data_cil, by="TranscriptID") %>%
  pivot_longer(col=c("P1_W2.R7", "P1_W4.R7", "P1_W6.R7", "P3_D1.R7", "P3_D2.R7", "P3_D3.R7"), names_to="sample", values_to="tmm_fpkm") %>%
  mutate(
    class=case_when(
      sample == "P1_W2.R7" ~ "Wild",
      sample == "P1_W4.R7" ~ "Wild",
      sample == "P1_W6.R7" ~ "Wild",
      sample == "P3_D1.R7" ~ "Domestic",
      sample == "P3_D2.R7" ~ "Domestic",
      sample == "P3_D3.R7" ~ "Domestic"
    )
  )

g2 <- ggplot(data_hipp,aes(x=class, y=log2(tmm_fpkm+1), color=class)) + 
  facet_wrap(vars(Name), ncol = 6, strip.position = "top") +
  geom_point() + 
  scale_color_manual(values=c("#b8cc97","#945a7f")) +
  ylab(expression(paste("Log"[2],"(TMM+1)"))) +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.title = element_blank(),
        strip.placement = "outside", strip.background = element_blank(), strip.text = element_text(face="italic"))

pdf("../Figures/FigureS2.pdf", w = 6, h = 6)
cowplot::plot_grid(g1, g2, scale=0.95, rel_heights=c(2, 4), ncol = 1, labels = "AUTO")
dev.off()
