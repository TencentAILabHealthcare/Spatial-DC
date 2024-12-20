library(Seurat)
library(dplyr)
library(SeuratDisk)
library(SPOTlight)


# ------------------------------------------------
start_time = Sys.time()

setwd("")

reference_data_dir = "01_data/reference/synthetic_noise_v2"
sp_data_dir = "01_data/simulations/synthetic_cellnum_noise_v2"

method = "SPOTlight"
output_dir = paste("03_output/exp_conditions_v2",method,sep="/")

if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

run_SPOTlight = function(reference, query, mgs, output_file_path=output_file_path){
    
    predictions <- SPOTlight(
      x = reference,
      y = query,      
      mgs = mgs,
      gene_id = "gene",
      group_id = "cluster",
      weight_id = "avg_log2FC",slot="data"      
        )

    temp_save_path = paste0(output_file_path, ".csv")
    write.csv( predictions$mat,file=temp_save_path)    
    
}

seeds <- c(0)
cells <- c(10, 15, 20) 
noises <- c(0, 0.25, 0.5, 0.75, 1.0)

celltype_key = "celltype"
for (reference_noise in noises){
  sc_file_path = paste0(reference_data_dir,"/","scp2021_1003_Reference_noise",as.integer(reference_noise * 100),".h5seurat")
  sc_meta_path = paste0(reference_data_dir,"/","scp2021_1003_Reference_noise",as.integer(reference_noise * 100),"_meta.csv")
  sc_rna <- LoadH5Seurat(sc_file_path,meta=F)

  meta=read.csv(sc_meta_path,header=T,row.names=1,sep=",")
  sc_rna@meta.data = meta
  Idents(object = sc_rna) <- meta[,celltype_key]

  # sc_rna <- SCTransform(sc_rna)
  sc_rna = NormalizeData(sc_rna, normalization.method = "LogNormalize", scale.factor = 10000)
  sc_rna <- FindVariableFeatures(sc_rna, selection.method = "vst", nfeatures = 1000)
  # all.genes <- rownames(pbmc)
  sc_rna <- ScaleData(sc_rna)

  cluster_markers_all <- Seurat::FindAllMarkers(object = sc_rna, logfc.threshold=0.25,test.use="t",
                                                assay = "RNA",
                                                slot = "data",
                                                verbose = TRUE, 
                                                only.pos = TRUE)
  for (seed in seeds) {
    for (cell in cells) {
      for (noise in noises) {
        
        print(paste("seed:", seed, "_cells:", cell, "_noise", as.integer(noise * 100)))
              
        output_file_path <- paste(output_dir, "/", method,"_reference_noise",as.integer(reference_noise * 100), "_seed", seed, "_cells", cell, "_noise", as.integer(noise * 100), sep = "")
        bulk_file_path <- paste(sp_data_dir, "/Simu_seed", seed, "_cells", cell, "_noise", as.integer(noise * 100), ".h5seurat", sep = "")

        tryCatch({
          spatial <- LoadH5Seurat(bulk_file_path,meta=F)        
          spatial <- NormalizeData(spatial, normalization.method = "LogNormalize", scale.factor = 10000)
          spatial <- FindVariableFeatures(spatial, selection.method = "vst", nfeatures = 1000)        
          spatial <- ScaleData(spatial)

          run_SPOTlight(sc_rna, query=spatial, mgs=cluster_markers_all, output_file_path=output_file_path)

        },error = function(e) {
          message("Error occurred:",conditionMessage(e))
        
        })

      }
    }
  }
  }


end_time <- Sys.time()
print("Total time consumption: seconds")
print(end_time - start_time)
