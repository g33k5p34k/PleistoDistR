# By: Ethan Gyllenhaal (egyllenhaal@unm.edu)
# Last updated: 12Apr2022
#
# R script for making box plots for Samao SLiM output

library(ggplot2)
library(gridExtra)
library(wesanderson)

setwd('/path/to/samoa_slim')

# load combined output data
fst_data <- read.csv("combined_output_fst.csv")
pi_data <- read.csv("combined_output_pi.csv")

## FST plotting

# set color values
color_values = wes_palette("Zissou1", 2, type="continuous")

# factor by dispersal value, color and fill by static vs dynamic
fst_SU <- ggplot(fst_data, aes(x = factor(Dispersal), y = savaii.upolu, color = Type, fill = Type))
# make boxplot with blank axes to add later, make solid lines and transleucent center
fst_plot_SU <- fst_SU + geom_boxplot(show.legend=FALSE, notch=FALSE, lwd = 0.3, outlier.size = 0.7) +  
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ylim(0, .36) + ggtitle("Savaii-Upolu FST") +
  scale_colour_manual(values = color_values) + scale_fill_manual(values = alpha(color_values, 0.3))

fst_TuU <- ggplot(fst_data, aes(x = factor(Dispersal), y = tutuila.upolu, color = Type, fill = Type))
fst_plot_TuU <- fst_TuU + geom_boxplot(show.legend=FALSE, notch=FALSE, lwd = 0.3, outlier.size = 0.7) +  
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ylim(0, .36) + ggtitle("Tutuila-Upolu FST") +
  scale_colour_manual(values = color_values) + scale_fill_manual(values = alpha(color_values, 0.3))

fst_TaU <- ggplot(fst_data, aes(x = factor(Dispersal), y = tau.upolu, color = Type, fill = Type))
fst_plot_TaU <- fst_TaU + geom_boxplot(show.legend=FALSE, notch=FALSE, lwd = 0.3, outlier.size = 0.7) +  
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  ylim(0, .36) + ggtitle("Tau-Upolu FST")+
  scale_colour_manual(values = color_values) + scale_fill_manual(values = alpha(color_values, 0.3))

## Pi plotting

pi_upolu <- ggplot(pi_data, aes(x = factor(Dispersal), y = upolu, color = Type, fill = Type))
pi_plot_upolu <- pi_upolu + geom_boxplot(show.legend=FALSE, notch=FALSE, lwd = 0.3, outlier.size = 0.7) + 
  theme_bw(base_size = 10) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  ylim(0,0.006) + ggtitle("Upolu Pi") + 
  scale_colour_manual(values = color_values) + scale_fill_manual(values = alpha(color_values, 0.3))


pi_tutuila <- ggplot(pi_data, aes(x = factor(Dispersal), y = tutuila, color = Type, fill = Type))
pi_plot_tutuila <- pi_tutuila + geom_boxplot(show.legend=FALSE, notch=FALSE, lwd = 0.3, outlier.size = 0.7) +  
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  ylim(0,0.006) + ggtitle("Tutuila Pi") +
  scale_colour_manual(values = color_values) + scale_fill_manual(values = alpha(color_values, 0.3))


pi_tau <- ggplot(pi_data, aes(x = factor(Dispersal), y = tau, color = Type, fill = Type))
pi_plot_tau <- pi_tau + geom_boxplot(show.legend=FALSE, notch=FALSE, lwd = 0.3, outlier.size = 0.7) + 
  theme_bw(base_size = 10) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  ylim(0,0.006) + ggtitle("Tau Pi") +
  scale_colour_manual(values = color_values) + scale_fill_manual(values = alpha(color_values, 0.3))

# grid arrange with bottom axis label, save png and pdf for later
slim_grid <- grid.arrange(fst_plot_SU, fst_plot_TuU, fst_plot_TaU, pi_plot_upolu, pi_plot_tutuila, pi_plot_tau, nrow=2)
ggsave(file = "prelim_box_slim.pdf", units = "in", width = 10, height = 7.5, dpi=300, slim_grid)
ggsave(file = "prelim_box_slim.png", units = "in", width = 10, height = 7.5, dpi=300, slim_grid)
