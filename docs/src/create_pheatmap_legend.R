library(tidyverse)
library(dichromat)
library(RColorBrewer)

r <- tibble(y = rep(1, 101), x = seq(from = 0, to = 1, by = 0.01))
heatmap_legend_palette <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(101)
legend <- ggplot(r, aes(y=y, x=factor(x), fill = factor(x))) +
  geom_tile(height = 0.5) +
  scale_fill_manual(values = heatmap_legend_palette) +
  scale_x_discrete(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1")) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 45, face = "bold"),
        plot.margin = margin(5,5,5,30),
        legend.position = "none")

ggsave("/Users/kayla/norm_for_RNAseq_coexp/plots/heatmap_legend.png", 
       plot = legend, width = 30, height = 2)