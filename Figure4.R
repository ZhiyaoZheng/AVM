############################################################
# Figure 4
############################################################

library(Seurat)
library(ggplot2)
library(ggforce)
library(dplyr)


# =========================
# Figure 4A
# =========================
data <- readRDS('3.sce_annotated1.rds')
table(data$ECcelltype)
Idents(data) <- data$ECcelltype
E4 <- data[, Idents(data) %in% c( 'E.04.ShearResponse' )]

genes <- c('ICAM1','PECAM1','SELE','VCAM1','SELP','ICAM2','MADCAM1')
features <- list("genes" = genes)

pdf('4A.adhesion doplot.pdf', width = 6, height = 2.5)
DotPlot(object = E4, features=features,
        col.min=0,col.max=1.5
)&
  theme_bw()& 
  scale_size(range=c(0,7))&
  theme(axis.title = element_blank(),
        axis.text.x = element_text(color = 'black',angle = 90, hjust = 1, vjust = 0.5,face = "bold"),
        axis.text.y = element_text(color = 'black', face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(), 
        strip.text = element_blank(),
        plot.margin=unit(c(1, 1, 1, 1),'cm'),
        panel.border = element_rect(color="black",size = 1.2, linetype="solid"),
        panel.spacing = unit(0.12, "cm"),
        legend.frame = element_rect(colour = "black"),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.5, "cm"),
        legend.title = element_text(color = 'black', face = "bold", size=9))& 
  scale_color_gradientn(colours = colorRampPalette(c('#b6d7e5','#cd7196'))(100))& 
  labs(tag = "                                       Leukocyte-endothelial adhesion genes")& 
  theme(plot.tag.position = c(0.3, 1.05), 
        plot.tag = element_text(size = 12,face = "bold"))& 
  guides(size=guide_legend(title="Proportion of\nexpressing cells"), 
         colour=guide_colorbar(title="Average\nexpression"))
dev.off()


# =========================
# Figure 4F
# =========================
dt <- read.csv('data.csv')
dt <- dt %>%
  mutate(
    ymax = cumsum(ratio),
    ymin = c(0, head(ymax, n=-1)),
    label_position = (ymax + ymin) / 2,
    label = paste0(ratio, "%")
  )
p <- ggplot(dt) +
  geom_rect(aes(ymax = ymax, ymin = ymin, 
                xmax = 4, xmin = 2, fill = Type),
            color = "white", size = 2, linejoin = "round") +
  geom_text(aes(x = 3.5, y = label_position, 
                label = label),
            color = "black", size =5, fontface = "bold") +
  geom_text(aes(x = 0, y = 0, label = "Immune Cell-Targeted Pathways"),
            size = 4, fontface = "bold", color = "#333333") +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  scale_fill_manual(values = c( "steelblue","#559d3c","#ed8431" ,"#d1362b", "#624094",'#2abfa3' )) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 11, face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )
p
ggsave('4F.GO Targeted Pathways.pdf')


# =========================
# Figure 4H
# =========================
p <- ggplot(plot_data, aes(x = count, y = target)) +
  geom_segment(
    aes(x = 0, xend = count, y = target, yend = target),
    color = "#4E79A7",
    size = 2,  
    alpha = 0.6
  ) +
  geom_point(aes(size = count), color = "#4E79A7", alpha = 0.7) +
  geom_text(
    aes(x = count, y = target, label = count),
    hjust = -1.7, 
    size = 3.5,    
    color = "black", 
    family = "sans"  
  ) +
  scale_x_continuous(limits = c(0, 35)) + 
  scale_y_discrete(limits = rev) +
  scale_size_continuous(range = c(3, 10)) +
  labs(
    x = "Number of Interactions (Count)",
    y = "T Cell Subset",
    title = "Increased Interaction Frequency in AVM", 
    size = "Count"
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black")
  )

pdf('4H.CellchatNumber.pdf', width = 6, height = 6)
p
dev.off()



# =========================
# Figure 4I
# =========================
dat.plot <- read.csv('immunecelltypefreq.csv')
compraisions.list=list(c("EP","AVM"))
g.colSet = list(group = c("EP" = "#80ACF9","AVM" = "#E99C93"))
p.list.perMcls <- lapply(unique(sort(dat.plot$Celltype)),function(mcls){
  
  text.size = 10
  text.angle = 45
  text.hjust = 1
  legend.position = "none"
  p <- ggboxplot(dat.plot[Celltype== mcls ,],x="Group",y="freq",
                 color = "Group", legend="none",title=mcls,
                 xlab="",ylab="percentage within Immune cell",
                 add = "jitter",outlier.shape=NA) +
    scale_color_manual(values=g.colSet$group) +
    stat_compare_means(label="p.format",comparisons=compraisions.list) +
    coord_cartesian(clip="off") +
    theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
          axis.ticks = element_line(color = "black"),
          axis.title = element_text(size = text.size,color ="black"), 
          axis.text = element_text(size=text.size,color = "black"),
          axis.text.x = element_text(angle = text.angle, hjust = text.hjust ),
          legend.position = legend.position,
          legend.text = element_text(size= text.size),
          legend.title= element_text(size= text.size),
          strip.background = element_rect(color="black",size= 1, linetype="solid") 
    )
  return(p)
})
p.prop = wrap_plots(p.list.perMcls,ncol = 5)
p.prop
ggsave(p.prop,filename = "4I.TNKCellIDstatis.pdf",width = 30,height = 22,units = "cm")

# =========================
# Figure 4N
# =========================
env_result <- read.csv('env_result.csv')
pdf("4N.CrossL.pdf", width = 5, height = 5)  #
plot(env_result, 
     main = "AVM Ripley's Cross L-function",
     xlab = "Distance r (µm)",
     ylab = 'L[E.04.ShearResponse,CD8T.02.Teff_GZMB](r)',
     legend = FALSE,
     col = c("red", "black", "transparent", "transparent"),
     lty = c(1, 2, 0, 0),
     lwd = c(2, 1, 0, 0),
     shadecol = "gray80") 
legend("topleft",
       legend = c("Observed L(r)", "Theoretical L(r)"),
       col = c("red","black"),
       lty = c(1, 3, 1),
       lwd = c(2, 1, 8),
       bg = "white")
text(x = max(env_result$r) * 0.7, 
     y = max(env_result$obs) * 0.95, 
     labels = paste("p =", format(p_value, digits = 3, scientific = FALSE)),
     cex = 1.2,  
     font = 2,   
     col = "black")
dev.off()
