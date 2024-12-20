library(CARD)
library(Seurat)
library(dplyr)
library(SeuratDisk)
library(MuSiC)

# ----------------------------
start_time = Sys.time()
setwd("")
reference_data_dir = ""
sp_data_dir = ""

method = "CARD"
output_dir = paste("03_output/exp_v1",method,sep="/")

if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

run_CARD = function(sc_exp, sc_meta, spatial_count, spatial_location, output_file_path=output_file_path){
    
      CARD_obj = createCARDObject(
        sc_count = sc_exp,
        sc_meta = sc_meta,
        spatial_count = spatial_count,
        spatial_location = spatial_location,
        ct.varname = "cellType",
        ct.select = unique(sc_meta$cellType),
        sample.varname = "sampleInfo",
        minCountGene = 0,
        minCountSpot = 0) 

      CARD_obj@spatial_location = spatial_location
      CARD_obj@spatial_countMat = spatial_count
      
      CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
      
      temp_save_path = paste0(output_file_path, ".csv")
      write.csv( CARD_obj@Proportion_CARD,file=temp_save_path)
       
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
  meta$cellID = row.names(meta)
  meta$cellType = meta$celltype
  meta$sampleInfo = "sample1"
  meta = meta[,c("cellID","cellType","sampleInfo")]

  sc_rna@meta.data = meta  
  sc_exp = sc_rna@assays$RNA@counts

  for (seed in seeds) {
    for (cell in cells) {
      for (noise in noises) {
        
        print(paste("seed:", seed, "_cells:", cell, "_noise", as.integer(noise * 100)))
              
        output_file_path <- paste(output_dir, "/", method,"_reference_noise",as.integer(reference_noise * 100), "_seed", seed, "_cells", cell, "_noise", as.integer(noise * 100), sep = "")
        bulk_file_path <- paste(sp_data_dir, "/Simu_seed", seed, "_cells", cell, "_noise", as.integer(noise * 100), ".h5seurat", sep = "")
        sp_loc_path <- paste(sp_data_dir, "/Simu_seed", seed, "_cells", cell, "_noise", as.integer(noise * 100), "_loc.csv", sep = "")
      
        sp_loc = read.csv(sp_loc_path,header=T,row.names=1)
        tryCatch({
          spatial <- LoadH5Seurat(bulk_file_path,meta=F)        
          spatial_count <- spatial@assays$RNA@counts# %>% as.data.frame()
          row.names(sp_loc) = colnames(spatial_count)   

          run_CARD(sc_exp, meta,  spatial_count, sp_loc, output_file_path=output_file_path)

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
