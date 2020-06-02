## Figure S3
library(tidyverse)

data <- read.table("../Data/DEGs_edgeR_dAF.tsv", head=T, stringsAsFactors=FALSE)
data <- data %>% pivot_longer(col=c("logFC_Amy", "logFC_Hyp", "logFC_Hipp", "logFC_ParTemp"), names_to="Reg", values_to="FC") %>%
  mutate(
    Reg=case_when(
      Reg == "logFC_Amy" ~ "Amy",
      Reg == "logFC_Hyp" ~ "Hyp",
      Reg == "logFC_Hipp" ~ "Hipp",
      Reg == "logFC_ParTemp" ~ "ParTemp"
    )
  ) %>%
  transform(Reg=factor(Reg,levels=c("Amy", "Hyp", "Hipp", "ParTemp"))) %>%
  select(TranscriptID,Reg,FC,dAF)

g <- ggplot(data, aes(Reg, abs(FC), color = Reg, fill = Reg)) +
  geom_boxplot(aes(alpha=dAF), outlier.colour = NA) +
  scale_color_manual(labels=c("Amygdala","Hypothalamus", "Hippocampus", "Parietal / temporal cortex"), values=c("#35B779FF", "#31688EFF", "#440154FF", "#FDE725FF")) +
  scale_fill_manual(labels=c("Amygdala","Hypothalamus", "Hippocampus", "Parietal / temporal cortex"), values=c("#35B779FF", "#31688EFF", "#440154FF", "#FDE725FF")) +
  scale_alpha_manual(values = c(No = 0.05, Yes = 1), labels=c("Genes without high dAF in CNCR", "Genes with high dAF in CNCR")) +
  coord_cartesian(ylim = c(0,1)) +
  ylab(expression(paste("|Log"[2],"(Fold Change)|"))) +
  theme_bw() +
  theme(legend.title=element_blank(),axis.title.x = element_blank(), axis.text.x = element_blank())

ggsave("../Figures/FigureS3.pdf",g,width=6.5,height=3)

## statistics
library(car)
data_tmp <- data %>% filter(Reg == "Amy")
Anova(lm(abs(data_tmp$FC)~data_tmp$dAF))
data_tmp <- data %>% filter(Reg == "Hyp")
Anova(lm(abs(data_tmp$FC)~data_tmp$dAF))
data_tmp <- data %>% filter(Reg == "Hipp")
Anova(lm(abs(data_tmp$FC)~data_tmp$dAF))
data_tmp <- data %>% filter(Reg == "ParTemp")
Anova(lm(abs(data_tmp$FC)~data_tmp$dAF))