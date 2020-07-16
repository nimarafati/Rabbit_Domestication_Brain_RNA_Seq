## Figure 3
library(tidyverse)
library(cowplot)

data_dAF<-read.table("../Data/raw_data/All_dAF_0.5over",sep="\t",head=T)
data_dAF_w<-read.table("../Data/dAF_window_10kb.txt",sep="\t",head=T)

data_dAF_x<-data_dAF[data_dAF$CHR=="chrX",]
data_dAF_w_x<-data_dAF_w[data_dAF_w$CHR=="chrX",]
g1<-ggplot(NULL) +
  geom_point(data=data_dAF_x,aes(x=POS,y=dAF),alpha=.1) +
  geom_line(data=data_dAF_w_x,aes(x=(START+END)/2,y=AVE_dAF)) +
  geom_segment(data=data_dAF_w_x, aes(x=87444396, xend=87443455, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5,col="#a22041") +
  annotate("text", x=87450000, y=.83, label="TCEAL2",col="#a22041", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=87272369, xend=87290971, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=87280000, y=.83, label="NXF2*", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=87377607, xend=87349431, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=87357000, y=.83, label="NXF2*", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=87403894, xend=87404235, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=87400000, y=.83, label="BEX5", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=87421515, xend=87422081, y=.87, yend=.87), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=87420000, y=.90, label="TCEAL6", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=87543715, xend=87555799, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=87540000, y=.83, label="ZMAT1", fontface="italic") +
  geom_hline(yintercept=quantile(data_dAF_w$AVE_dAF,probs=0.99),linetype="dashed",col="#a22041") +
  coord_cartesian(xlim = c(87250000,87550000)) +
  xlab("Chromosome X (OryCun2)") +
  ylab("dAF") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw(base_size = 11, base_family = "")

data_dAF_x<-data_dAF[data_dAF$CHR=="chr13",]
data_dAF_w_x<-data_dAF_w[data_dAF_w$CHR=="chr13",]
g2<-ggplot(NULL) +
  geom_point(data=data_dAF_x,aes(x=POS,y=dAF),alpha=.1) +
  geom_line(data=data_dAF_w_x,aes(x=(START+END)/2,y=AVE_dAF)) +
  geom_segment(data=data_dAF_w_x, aes(x=36108502, xend=36089996, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5,col="#a22041") +
  annotate("text", x=36098502, y=.83, label="NTRK1",col="#a22041", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=36072143, xend=36059702, y=.87, yend=.87), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=36066143, y=.90, label="PEAR1", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=36059409, xend=36045853, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=36041853, y=.83, label="LRRC71", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=35925920, xend=36045492, y=.87, yend=.87), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=35991492, y=.90, label="ARHGEF11", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=36109982, xend=36127839, y=.87, yend=.87), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=36118839, y=.90, label="INSRR", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=36145339, xend=36153821, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=36154339, y=.83, label="SH2D2A", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=36187459, xend=36155207, y=.87, yend=.87), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=36170459, y=.90, label="PRCC", fontface="italic") +
  geom_hline(yintercept=quantile(data_dAF_w$AVE_dAF,probs=0.99),linetype="dashed",col="#a22041") +
  coord_cartesian(xlim = c(35900000,36200000)) +
  xlab("Chromosome 13 (OryCun2)") +
  ylab("dAF") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw(base_size = 11, base_family = "")

data_dAF_x<-data_dAF[data_dAF$CHR=="chr10",]
data_dAF_w_x<-data_dAF_w[data_dAF_w$CHR=="chr10",]
g3<-ggplot(NULL) +
  geom_point(data=data_dAF_x,aes(x=POS,y=dAF),alpha=.1) +
  geom_line(data=data_dAF_w_x,aes(x=(START+END)/2,y=AVE_dAF)) +
  geom_segment(data=data_dAF_w_x, aes(x=5691866, xend=6610935, y=.9, yend=.9),size = 2, col="#a22041") +
  geom_segment(data=data_dAF_w_x, aes(x=6671117, xend=7033178, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5,col="#a22041") +
  annotate("text", x=6850000, y=.85, label="DNAH11",col="#a22041", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=6560952, xend=6645083, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=6600000, y=.85, label="SP4", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=7046368, xend=7031598, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=7040000, y=.85, label="CDCA7L", fontface="italic") +
  geom_hline(yintercept=quantile(data_dAF_w$AVE_dAF,probs=0.99),linetype="dashed",col="#a22041") +
  coord_cartesian(xlim = c(6500000,7100000)) +
  xlab("Chromosome 10 (OryCun2)") +
  ylab("dAF") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw(base_size = 11, base_family = "")

data_dAF_x<-data_dAF[data_dAF$CHR=="chr1",]
data_dAF_w_x<-data_dAF_w[data_dAF_w$CHR=="chr1",]
g4<-ggplot(NULL) +
  geom_point(data=data_dAF_x,aes(x=POS,y=dAF),alpha=.1) +
  geom_line(data=data_dAF_w_x,aes(x=(START+END)/2,y=AVE_dAF)) +
  geom_segment(data=data_dAF_w_x, aes(x=45556643, xend=45556161, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5,col="#a22041") +
  annotate("text", x=45560000, y=.83, label="RPL21",col="#a22041", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=45429244, xend=45862102, y=.87, yend=.87), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=45660000, y=.90, label="PTPRD", fontface="italic") +
  geom_segment(data=data_dAF_w_x, aes(x=46416417, xend=46420263, y=.8, yend=.8), arrow = arrow(length = unit(0.2, "cm")),size = .5) +
  annotate("text", x=46418000, y=.85, label="TMEM261", fontface="italic") +
  geom_hline(yintercept=quantile(data_dAF_w$AVE_dAF,probs=0.99),linetype="dashed",col="#a22041") +
  coord_cartesian(xlim = c(45300000,46000000)) +
  xlab("Chromosome 1 (OryCun2)") +
  ylab("dAF") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_bw(base_size = 11, base_family = "")

pdf("../Figures/Figure3.pdf",w=8,h=6)
cowplot::plot_grid(g1, g2, g3, g4, scale=0.95, ncol = 2, labels = "AUTO")
dev.off()
