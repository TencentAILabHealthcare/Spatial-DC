# Analysis of co-localization with predicted proteomic profiles
# Yuan
# 20241213
# Fig. 5h

library(dplyr)

rm(list=ls())
setwd("")

data_df = read.csv("SpatialDC/reconstruct.csv",header=T,row.names=1)

data_df_filter = data_df %>% t() %>% as.data.frame()
data_df_filter = data_df_filter[which(rowSums(data_df_filter) > 0),]

data_df_filter = data_df_filter %>% t() %>% as.data.frame()

corr_df = cor(data_df_filter, method="pearson")

# set break range for pheatmap 
range_all <- range(c(-1, 1))    
my_palette <- colorRampPalette(c("navy", "white","firebrick3"))(n=100)
breaks = seq(range_all[1], range_all[2], length.out = 101)

p1=pheatmap::pheatmap(corr_df, color=my_palette, breaks=breaks)

pdf("figures/panel_corr_spatialDC.pdf",height=5,width=5.5)
print(p1)
dev.off()

write.csv(corr_df[p1$tree_row$labels[p1$tree_row$order],p1$tree_col$labels[p1$tree_col$order]],"analysis/panel_co_loc.csv")
