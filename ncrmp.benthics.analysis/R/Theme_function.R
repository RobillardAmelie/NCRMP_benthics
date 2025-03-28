theme_Publication <- function(base_size=14, base_family = 'Arial') {
  library(ggthemes)
  
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            panel.spacing.x = unit(8, "mm"),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = "grey50", size = 1),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(),
            axis.line.x = element_line(colour="grey50", size = 1),
            axis.line.y = element_line(colour="grey50", size = 1),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_line(colour="#f0f0f0"),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.margin = margin(0,0,0,0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_blank(),
            strip.text = element_text(face="plain")
    ))
  
}
