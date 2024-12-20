# Converting from AnnData to Seurat via h5Seurat
# 20241119
# Yuan
library(Seurat)
library(SeuratData)
library(SeuratDisk)

# ------------------------------
# h5ad to h5seurat
files = list.files("xxxxx",pattern = "*.h5ad",full.names = T)
for (file in files){
        print(file)
        Convert(source = file, dest="h5seurat", overwrite=T,meta=F)
}


# -----------------------------
# Seurat Obj to h5Seurat
# h5Seurat to h5ad
cite_seurat = readRDS("xxx.rds")
rna_counts <- cite_seurat@assays$ADT@counts
cite_seurat@meta.data$celltype = cite_seurat@meta.data$annotation_figure_1
seurat_obj <- CreateSeuratObject(counts = rna_counts, meta.data=cite_seurat@meta.data) # create Seurat Obj.

# Convert
SaveH5Seurat(seurat_obj, "20220215_tonsil_atlas_cite_seurat_obj.h5Seurat")
Convert("20220215_tonsil_atlas_cite_seurat_obj.h5Seurat", dest = "h5ad")


