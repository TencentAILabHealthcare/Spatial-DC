library(dplyr)

setwd("/aaa/zihanwu/yyyli2/project1_spatial_deconv/experiment_records/20240722/exp21/HumanTonsil2023/01_cell_perc/03_benchmark_methods/exp_ref46/SpatialDC")

rm(list=ls())
data_df = read.csv("SpatialDC.csv",header=T,row.names=1)

corr_df = cor(data_df, method="pearson")

# set break range for pheatmap 
range_all <- range(c(-1, 1))    
my_palette <- colorRampPalette(c("navy", "white","firebrick3"))(n=100)
breaks = seq(range_all[1], range_all[2], length.out = 101)


p1=pheatmap::pheatmap(corr_df, color=my_palette, breaks=breaks)

pdf("corr_spatialDC.pdf")
print(p1)
dev.off()
