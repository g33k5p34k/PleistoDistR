# By: Ethan Gyllenhaal (egyllenhaal@unm.edu)
# Last updated: 12Apr2022
#
# R script for making box plots for Fiji msprime output

library(ggplot2)
library(gridExtra)
library(wesanderson)

setwd('/path/to/fiji_msprime')

# load combined output data
data <- read.csv("combined_output.csv")

# set color values
color_values = wes_palette("Zissou1", 2, type="continuous")

# factor by dispersal value, color and fill by static vs dynamic
fst <- ggplot(data, aes(x= factor(dispersal), y = fst, color = type, fill = type))
# make boxplot with blank axes to add later, make solid lines and transleucent center
fst_plot <- fst + geom_boxplot(show.legend=FALSE, notch=FALSE, lwd = 0.3, outlier.size = 0.7) + 
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ggtitle("Viti Levu-Kadavu FST") +
  scale_color_manual(values = color_values) + scale_fill_manual(values = alpha(color_values, 0.3))

piV <- ggplot(data, aes(x= factor(dispersal), y = vitiPi, color = type, fill = type))
piV_plot <- piV + geom_boxplot(show.legend=FALSE, notch=FALSE, lwd = 0.3, outlier.size = 0.7) + 
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ylim(0, 0.19) + ggtitle("Viti Levu Pi") +
  scale_color_manual(values = color_values) + scale_fill_manual(values = alpha(color_values, 0.3))

piK <- ggplot(data, aes(x= factor(dispersal), y = kadPi, color = type, fill = type))
piK_plot <- piK + geom_boxplot(show.legend=FALSE, notch=FALSE, lwd = 0.3, outlier.size = 0.7) + 
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ylim(0, 0.19) + ggtitle("Kadavu Pi") +
  scale_color_manual(values = color_values) + scale_fill_manual(values = alpha(color_values, 0.3))

# grid arrange with bottom axis label, save png and pdf for later
msp_grid <- grid.arrange(fst_plot, piV_plot, piK_plot, nrow=1, bottom="Mean Dispersal (km)")
ggsave(file = "prelim_box_msp.pdf", units = "in", width = 10, height = 4, dpi=300, msp_grid)
ggsave(file = "prelim_box_msp.png", units = "in", width = 10, height = 4, dpi=300, msp_grid)
