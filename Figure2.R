############################################################
# Figure 2 
############################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(scRNAtoolVis)
library(ggpubr)
library(ggimage)
library(Seurat)

# =========================
# Figure 2A
# =========================
df <- read.csv('Top10pathway.csv')

df$neg_log10_padj <- -log10(df$p.adjust)
df <- df %>% arrange(desc(NES))
df$Description <- factor(df$Description, levels = rev(df$Description))

p <- ggplot(df, aes(x = NES, y = Description)) +
  geom_segment(aes(x = 1.5, xend = NES, y = Description, yend = Description), 
               color = "grey70", size = 0.8) +
  geom_point(aes(color = neg_log10_padj, shape = Source), size = 4) +
  scale_color_gradient(low = "#ffc38f", high = "#e84609", 
                       limits = c(1.5, 2.0) ,
                       name = expression(-log[10](p.adjust))) +
  scale_shape_manual(values = c(16, 17, 15,18), name = "Source") +
  scale_x_continuous(limits = c(1.5, 2.1), breaks = seq(1.5, 2.1, by = 0.1)) +
  labs(
    title = "Top 10 Enriched Pathways",
    x = "Normalized Enrichment Score",
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 15)),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, color = "black", hjust = 1),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.line.y = element_line(color = "black", size = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )
p
ggsave('2A.Top10All.pdf',height = 5,width = 11)


# =========================
# Figure 2B
# =========================

seu=readRDS("3.sce_annotated1.rds")

VlnExp <- function(object,
                      group,
                      group_order,
                      features,
                      comparisons,
                      cols=NULL,
                      pieSize=NULL){
  
  #seurat obj
  meta <- object@meta.data
  if(group %in% colnames(meta)){
    Idents(object) <- group
    Idents(object) <- factor(Idents(object), levels = group_order)
  }else{
    print("Your group name does not exist, Please provide correct name in the colnames of metadata")
    
  }
  exp = FetchData(object, vars = features)
  colnames(exp) <- 'gene'
  exp$Group <- object@meta.data[,group]
  fearures_exp <- list()
  for (i in 1:length(group_order)) {
    exp1 = subset(exp, Group==group_order[i])
    colnames(exp1)[1] <- "gene"
    fearures_exp[[i]] <- exp1
    names(fearures_exp)[i] <- group_order[i]
  }
  pct_list <- list()
  for (i in 1:length(group_order)) {
    pct <- length(fearures_exp[[i]]$gene[fearures_exp[[i]]$gene>0]) / length(fearures_exp[[i]]$gene) %>% as.data.frame()
    pct_list[[i]] <- pct
    names(pct_list)[i] <- group_order[i]
  }
  exp_pct <- do.call(rbind, pct_list)
  colnames(exp_pct) <- "exp_gene"
  exp_pct$non_exp <- 1-exp_pct$exp_gene
  print(exp_pct)
  #plot pie
  if(is.null(cols)){
    colors_map = c("#FF5744","#208A42","#F98400", "#5BBCD6",'#7F3C8D' ,'#11A579', '#3969AC','#E73F74')
    cols <- colors_map[1:length(group_order)]
  }else{
    cols <- cols
  }
  plot_pie <- function(i) {
    df1 <- gather(exp_pct[i,], type, value, 1:2)
    df1$labels <- df1$value
    df1$labels <- scales::percent(df1$value,accuracy=0.1)
    df1$labels[2] <- ''
    ggplot(df1, aes(x= '', value, fill=type)) +
      geom_col(color='black') + #饼图边设置为黑色
      coord_polar(theta = 'y') +
      theme_void() + 
      theme_transparent() +
      theme(legend.position = "none")+
      geom_text(aes(label = labels), 
                position = position_stack(vjust = 0.5), size=4)+
      scale_fill_manual(values = c(cols[i],"grey80"))#修改饼图填充颜色
  }
  #data for pie
  exp_pct$pie <- lapply(1:nrow(exp_pct), plot_pie)
  exp_pct$x <- 1:length(group_order)
  exp_pct$y <- max(exp$gene)+2
  if(is.null(pieSize)){
    pieSize <- 1.5
  }else{
    pieSize <-pieSize
  }
  exp_pct$width <- pieSize
  exp_pct$height <-pieSize
 #plot Vln
  label.y.pos =  seq(max(exp$gene), max(exp$gene)+100, by = 0.3)
 p = VlnPlot(object, features = features)&
    theme_bw()&
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(color = 'black',face = "bold", size = 12),
          axis.text.y = element_text(color = 'black', face = "bold"),
          axis.title.y = element_text(color = 'black', face = "bold", size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color="black",size = 1.2, linetype="solid"),#修改边框大小
          panel.spacing = unit(0.12, "cm"),
          plot.title = element_text(hjust = 0.5, face = "bold.italic"),#字体加粗斜体
          legend.position = 'none')&
    stat_compare_means(method="wilcox.test",hide.ns = F, #显著性检验
                       comparisons = comparisons,
                       label="p.format",
                       bracket.size=0.8,
                       tip.length=0,
                       size=5,
                       vjust = 0.6,
                       label.y = label.y.pos)&
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))& #修改y轴，让其适应显著性检验标注
    scale_fill_manual(values = cols)&
    geom_subview(aes(x=x, y=y, 
                     subview=pie,width = width,
                     height = height), data=exp_pct)&
    ylim(0,max(exp$gene)+3)
  return(p)
}

pdf("2B.PIEZO1pieVlinpvalue.pdf", width = 5, height = 7)
VlnExp(object = seu, group="E4orno",
          group_order=c("E.04.ShearResponse","others"),
          features="PIEZO1",comparisons= list(c("E.04.ShearResponse","others")))
dev.off()



# =========================
# Figure 2C
# =========================
sce_cca=readRDS("3.sce_annotated1.rds")
markers = c('PIEZO1')

pdf("Fig2C.pdf", width =4.5, height = 5)
jjDotPlot(object=sce_cca,
          gene=markers,
          id='group',
          ytree = F)

dev.off()