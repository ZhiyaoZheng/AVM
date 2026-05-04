############################################################
# Figure 5
############################################################

library(Seurat)
library(ggplot2)
library(ggforce)
library(dplyr)
library(ggpubr)


# =========================
# Figure 5A
# =========================
data <- read.csv("AllORdata.csv", header = T)
target_cells <- c('Endothelial cell', 'Smooth muscle cell', 'Fibroblast', 'Pericyte')
filtered_data <- data[rid %in% target_cells]
head(filtered_data)
p1 <- ggplot(filtered_data, aes(cid, rid))+ 
  geom_tile(aes(fill = OR), colour = NA, size = 15)+
  scale_fill_gradientn(name='OR',
                       limit = c(0,3),
                       breaks = seq(0,3),
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
ggsave(p1,filename = "5A.ORsVascular.pdf",height = 6, width = 6)



# =========================
# Figure 5B
# =========================
data <- read.csv('data.csv')
p <- ggplot(data, aes(x = PeriinEC, y =Teff, color = Sample)) +
  geom_point(size = 5,color = "#507fbd") +
  geom_smooth(method = "lm", color = "black", fill = "#99b3c6",
              size = 1, se = T) + 
  labs(title = "Rupture AVM",
       x = "Pericyte coverage",
       y = "Teff_GZMB percentage in All cell") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 12), 
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 12) 
  ) +
  annotate("text", 
           x = -Inf,  
           y = Inf,  
           label = paste("R =", correlation_rounded,
                         "\nP =", p_value_rounded),
           color = "black", 
           size = 5, 
           hjust = -0.1,  
           vjust = 1.5)  

p
pdf('5B.pericyte-teff.pdf', width = 5, height = 5)
p
dev.off()


# =========================
# Figure 5D
# =========================
pdf("5D.AVMpericyteCrossL.pdf", width = 5, height = 5)  #
plot(env_result, 
     main = "AVM Ripley's Cross L-function",
     xlab = "Distance r (µm)",
     ylab = 'L[Pericyte,CD8T.02.Teff_GZMB](r)',
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

