############################################################
# Figure 1 
############################################################


library(Seurat)
library(ggplot2)
library(ggSCvis)
library(cowplot)
library(patchwork)
library(Matrix)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggthemes)
library(tidyverse)


# =========================
# Figure 1B
# =========================
# Load data
sce.all=readRDS("0.end.rds")
Idents(sce.all) <- sce.all$celltype
meta = sce.all@meta.data
colnames(meta)
col_df = data.frame(name = unique(meta$celltype)) %>%arrange(name)
col_df
mycol=c('#91b5a1','#347852','#d194a7','#de9d3d',"#A8A2D2",'#425785','#92b8da',
                 '#b5aa82','#705235','#c65341','#9f2b39','#61adc7')
mycol = setNames(mycol,col_df$name)
mycol

rs = sce.all@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = sce.all@meta.data$celltype) 
p <- ggplot(rs,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +
  geom_point(size = 1 , alpha =1 ) +
  scale_color_manual(values = mycol)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(),
    legend.key=element_rect(fill='white'), 
    legend.text = element_text(size=10),
    legend.key.size=unit(1,'cm') ) + 
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_segment(aes(x = min(rs$UMAP_1) , y = min(rs$UMAP_2) ,
                   xend = min(rs$UMAP_1) +3, yend = min(rs$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(rs$UMAP_1)  , y = min(rs$UMAP_2)  ,
                   xend = min(rs$UMAP_1) , yend = min(rs$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(rs$UMAP_1) +1.5, y = min(rs$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(rs$UMAP_1) -1, y = min(rs$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
p

pdf('1B.umap.pdf', width =9, height = 6)
p
dev.off()


# =========================
# Figure 1C
# =========================
# load data
sce.all=readRDS("3.sce_annotated1.rds")

Idents(sce.all) <- sce.all$ECcelltype
meta = sce.all@meta.data
colnames(meta)
col_df = data.frame(name = unique(meta$ECcelltype)) %>%arrange(name)
col_df
mycol=c('#81ad88','#f9bc53','#466eac','#FF4500', '#e1a991','#a94c0e')
mycol = setNames(mycol,col_df$name)
mycol

rs = sce.all@reductions$umap@cell.embeddings %>%  
  as.data.frame() %>% 
  cbind(cell_type = sce.all@meta.data$ECcelltype) 
rs
p <- ggplot(rs,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +
  geom_point(size = 1 , alpha =1 ) +
  scale_color_manual(values = mycol)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        axis.title = element_blank(),  
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), 
        plot.background=element_rect(fill="white"))+
  theme(
    legend.title = element_blank(),
    legend.key=element_rect(fill='white'), 
    legend.text = element_text(size=10),
    legend.key.size=unit(1,'cm') ) + 
  guides(color = guide_legend(override.aes = list(size=5)))+
  geom_segment(aes(x = min(rs$UMAP_1) , y = min(rs$UMAP_2) ,
                   xend = min(rs$UMAP_1) +3, yend = min(rs$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(rs$UMAP_1)  , y = min(rs$UMAP_2)  ,
                   xend = min(rs$UMAP_1) , yend = min(rs$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(rs$UMAP_1) +1.5, y = min(rs$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(rs$UMAP_1) -1, y = min(rs$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
p

pdf('1C.ECumap.pdf', width =9, height = 6)
p
dev.off()



# =========================
# Figure 1D
# =========================

data <- read.csv("ECORdata.csv", header = T)

p1 <- ggplot(data, aes(cid, rid))+ 
  geom_tile(aes(fill = OR), colour = NA, size = 15)+
  scale_fill_gradientn(name='OR',
                       limit = c(0,4),
                       breaks = seq(0,4),
                       na.value = "#F6C63C",
                       colours= c("#282A62","#329845","#F6C63C"))+
  theme_minimal() +
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.x = element_text(size = 14,color = 'black'),
        axis.text.y = element_text(size = 14,color = 'black')) + 
  scale_y_discrete(position = "right")
p1
ggsave(p1,filename = "1D.ECORs.pdf",height = 6, width = 5)



# =========================
# Figure 1F
# =========================
Data <- read.csv('4.AUCdata.csv')
mycol = c('#81ad88','#f9bc53','#466eac','#FF4500', '#e1a991','#a94c0e'
          
)
library(ggplot2)
p <- ggplot(Data,aes(x=cell_type,y=Score,fill=cell_type))+
  geom_violin(trim=FALSE,color="white")+
  geom_boxplot(width=0.2,position=position_dodge(0.9))+
  geom_point(aes(x=cell_type,y=mean),pch=19,position=position_dodge(0.9),size=1.5)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),
                width=0.1,
                position=position_dodge(0.9),
                color="black",
                alpha=0.7,
                size=0.5)+
  scale_fill_manual(values=mycol)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",size=20),
        axis.text.y=element_text(size=16,face="plain"),
        axis.title.y=element_text(size=20,face="plain"),
        panel.border=element_blank(),axis.line=element_line(colour="black",size=1),
        legend.text=element_text(colour="black",
                                 size=8),
        legend.title=element_text(colour="black",
                                  size=18),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  ylab("")+xlab("4.AUCell score")
p

library(ggpubr)

my_comparisons = list(c('E.04.ShearResponse','E.01.Artery'),
                      c('E.04.ShearResponse','E.02.Capillary'),
                      c('E.04.ShearResponse','E.03.Vein'),
                      c('E.04.ShearResponse','E.05.EndMT'),
                      c('E.04.ShearResponse','E.06.MHCIIHigh')
)

pdf('1E.AUCell.pdf', width = 7, height = 5)
p+stat_compare_means(comparisons = my_comparisons, label = "p.format")
dev.off()
