library(Seurat)
library(dplyr)
library(SeuratDisk)

# --------------------------------------------
start_time = Sys.time()

setwd("")

sc_file_path = "01_data/reference/scp2021_1003_Reference.h5seurat"
sc_meta_path = "01_data/reference/scp2021_1003_Reference_meta.csv"
celltype_key = "celltype"

sp_data_dir = "01_data/simulations/synthetic_cellnum_noise_Lung5/"

method = "Seurat" # change - 2
output_dir = paste("03_benchmark_methods/exp_ct3_Lung5",method,sep="/")

if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

sc_rna <- LoadH5Seurat(sc_file_path,meta=F)
meta=read.csv(sc_meta_path,header=T,row.names=1,sep=",")
sc_rna@meta.data = meta

# sc_rna <- SCTransform(sc_rna)
sc_rna = NormalizeData(sc_rna, normalization.method = "LogNormalize", scale.factor = 10000)
sc_rna <- FindVariableFeatures(sc_rna, selection.method = "vst", nfeatures = 1000)
# all.genes <- rownames(pbmc)
sc_rna <- ScaleData(sc_rna)

run_Seurat = function(reference, query, celltype_key=celltype_key, output_file_path=output_file_path){
    anchors <- FindTransferAnchors(reference=reference, query = query, normalization.method = 'LogNormalize',k.score=10,dims=1:10,npcs=10)#,k.score=10,npcs=5,dims=1:5)#,k.score=10,npcs=5,dims = 1:5)
    predictions <- TransferData(anchorset = anchors, refdata = reference@meta.data[,celltype_key],k.weight=20)#,dims = 1:5,k.weight=5)# if error, anno  # ,k.weight=20, related to anchors
   
    celltype_names = unique(reference@meta.data[,celltype_key])

    names(predictions) = c("id",celltype_names,"score")

    temp_save_path = paste0(output_file_path, ".csv")
    write.csv(predictions,file=temp_save_path)
}

seeds <- c(0)
cells <- c(10, 15, 20) 
noises <- c(0, 0.25, 0.5, 0.75, 1.0)

for (seed in seeds) {
  for (cell in cells) {
    for (noise in noises) {
      
      print(paste("seed:", seed, "_cells:", cell, "_noise", as.integer(noise * 100)))
            
      output_file_path <- paste(output_dir, "/", method, "_seed", seed, "_cells", cell, "_noise", as.integer(noise * 100), sep = "")
      bulk_file_path <- paste(sp_data_dir, "/Simu_seed", seed, "_cells", cell, "_noise", as.integer(noise * 100), ".h5seurat", sep = "")

      tryCatch({
        spatial <- LoadH5Seurat(bulk_file_path,meta=F)        
        spatial <- NormalizeData(spatial, normalization.method = "LogNormalize", scale.factor = 10000)
        spatial <- FindVariableFeatures(spatial, selection.method = "vst", nfeatures = 1000)        
        spatial <- ScaleData(spatial)

        run_Seurat(reference=sc_rna, query=spatial, 
          celltype_key=celltype_key, output_file_path=output_file_path)

      },error = function(e) {
        message("Error occurred:",conditionMessage(e))
      
      })

    }
  }
}

end_time <- Sys.time()
print("Total time consumption: seconds")
print(end_time - start_time)
