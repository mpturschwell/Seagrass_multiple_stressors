library(patchwork)

load(file = "Scripts/Data/population_figs_for_paper.RDA")
load(file = "Scripts/Data/consumer-resource_figs_for_paper.RDA")


Fig3 <- f3A /f3B
Fig3 <- Fig3 + plot_annotation(tag_levels = 'A')
Fig3
ggsave(path = "Plots", filename = paste0(mytime, "_Figure3.tiff"), Fig3 , width = 8, height = 8, units = c("in"), dpi = 300)


Figure4 <- (f4A + f4C) / (f4B + f4D)
Figure4 <- Figure4 + plot_annotation(tag_levels = 'A')
Figure4
ggsave(path = "Plots", filename = paste0(mytime, "_Figure4.tiff"), Figure4 , width = 8, height = 8, units = c("in"), dpi = 300)
