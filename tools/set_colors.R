library(ggsci)
library(ggplot2)
library(dplyr)

# pal_npg("nrc")(10)[-8]
# [1] "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" [7] "#91D1C2FF" "#7E6148FF" "#B09C85FF

benchmark_method = c("Tangram","Cell2location","Stereoscope", "DestVI","Seurat","SPOTlight","SpatialDWLS","CARD","SpatialDC")
method_color = c("#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF","#91D1C2FF", "#7E6148FF", "#B09C85FF", "#E64B35FF")
scales::show_col(method_color) # overview of the colors


p1 = plot_df %>% ggplot(aes(x=Noise,y=val_ccc,color=Method)) + 
    geom_line() + geom_point() + facet_wrap(~CellType) + 
    scale_x_continuous(limits=c(20,80),breaks=seq(25,75,25),label=seq(25,75,25),expand = c(0,0.2)) +
    theme_bw()+ 
    theme(plot.background = element_blank(),legend.position = "bottom",
            legend.background = element_blank(),#panel.grid = element_blank(),axis.text.x = element_blank(),
            text = element_text(size = 10)#,axis.ticks.x.bottom = element_blank()
    )+            
    scale_color_manual(values = method_color) + # using this to set color manually
    scale_y_continuous(limits=c(-0.1,1),breaks=seq(0,1,0.25),label=seq(0,1,0.25),expand = c(0,0.02))+
    labs(x = "External Gaussian Noise", y = "CCC") 
