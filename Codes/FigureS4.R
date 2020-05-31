## Figure S4
library(tidyverse)
#library(devtools)
#install_github("trinker/plotflow")
library(plotflow)

data_dAF_TSS<-read.csv("../Data/All_dAF_TSS100kb_DEGs.txt",sep="\t",head=T)
head(data_dAF_TSS)
nT<- data_dAF_TSS %>% filter(Genes != "-") %>% nrow() #963548

data_dAF_TSS2<-read.csv("../Data/All_dAF_TSS100kb_DEGs_res.txt",sep="\t",head=T)
#sum(data_dAF_TSS2$Num[data_dAF_TSS2$Region=="Amy"])
#sum(data_dAF_TSS2$Num[data_dAF_TSS2$Region=="Amy" & data_dAF_TSS2$DEGs=="DEGs"])
data_dAF_TSS2$Exp_Num <- data_dAF_TSS2$Total_Reg_deg * (data_dAF_TSS2$Total / nT)
data_dAF_TSS2$Mvalue <- log2(data_dAF_TSS2$Num/data_dAF_TSS2$Exp_Num)
data_dAF_TSS2_degs <- data_dAF_TSS2 %>% filter(DEGs == "DEGs")
data_dAF_TSS2_nodegs <- data_dAF_TSS2 %>% filter(DEGs == "nonDEGs")

rg<-ggplot(NULL) +
  geom_line(data=data_dAF_TSS2_degs,aes(x=Bin,y=log2(Num/Exp_Num),col=Region,group=Region),size=2,alpha=.8) +
  geom_line(data=data_dAF_TSS2_nodegs,aes(x=Bin,y=log2(Num/Exp_Num),col=Region,group=Region),size=2,linetype="twodash",alpha=.8) +
  scale_color_viridis_d() +
  labs(x = "dAF bins", y = "M value") +
  theme_classic(base_size = 11, base_family = "") +
  theme(legend.title=element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))

lg <- ggplot(data_dAF_TSS2_degs[data_dAF_TSS2_degs$Region=="Amy",], aes(x=Bin, y=Total)) +
  geom_line(group=1,col="black",size=2) + 
  scale_y_continuous(labels = scales::comma) +
  labs(x = "dAF bins", y = "Number of SNPs within 100kb from TSS") +
  theme_bw(base_size = 11, base_family = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

g <- ggdual_axis(lg, rg)

ggsave("../Figures/FigureS4.pdf", g, width=5.5, height=5)
