library(Giotto)
library(dplyr)
instrs = createGiottoInstructions(save_plot = FALSE,show_plot = FALSE,python_path="/usr/bin/python3")


# -----------------------------------------------------
start_time = Sys.time()

setwd("")

reference_data_dir = "01_data/reference/synthetic_noise_v2"
sp_data_dir = "01_data/simulations/synthetic_cellnum_noise_v2"
# output dir
method = "SpatialDWLS" # change - 2
output_dir = paste("03_output/exp_conditions_v2",method,sep="/")

if (!file.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


run_SpatialDWLS = function(sign_matrix, spatial_locs, spatialRNA,celltype_key=celltype_key, output_file_path=output_file_path){
    
    ##Generate Giotto objects and cluster spots
    st_data <- createGiottoObject(raw_exprs = t(spatialRNA),
                                spatial_locs = spatial_locs,
                                instructions = instrs)

    st_data <- filterGiotto(gobject = st_data,
                            expression_threshold = 0,
                            gene_det_in_min_cells = 0,
                            min_det_genes_per_cell = 0,
                            expression_values = c('raw'),
                            verbose = T)

    st_data <- normalizeGiotto(gobject = st_data)
    st_data <- calculateHVG(gobject = st_data)
    gene_metadata = fDataDT(st_data) 
    st_featgenes = gene_metadata[hvg == 'yes']$gene_ID 

    st_data <- runPCA(gobject = st_data, genes_to_use = st_featgenes, scale_unit = F)
    st_data <- createNearestNetwork(gobject = st_data, dimensions_to_use = 1:10, k = 10)
    st_data <- doLeidenCluster(gobject = st_data, resolution = 0.4, n_iterations = 1000)
    st_data <- runDWLSDeconv(gobject = st_data, sign_matrix = sign_matrix)
    
    temp_save_path = paste0(output_file_path, ".csv")
    write.csv(st_data@spatial_enrichment$DWLS,file=temp_save_path,row.names=FALSE)
}

seeds <- c(0)
cells <- c(10, 15, 20) 
noises <- c(0, 0.25, 0.5, 0.75, 1.0)

celltype_key = "celltype"
for (reference_noise in noises){
  sc_file_path = paste0(reference_data_dir,"/","scp2021_1003_Reference_noise",as.integer(reference_noise * 100),"_expression.csv")
  sc_meta_path = paste0(reference_data_dir,"/","scp2021_1003_Reference_noise",as.integer(reference_noise * 100),"_meta.csv")  

  sc_count = read.csv(sc_file_path, header=T,row.names=1)

  sc_meta = read.csv(sc_meta_path, header=T,row.names=1)
  sc_data <- createGiottoObject(raw_exprs = sc_count,instructions = instrs)
  sc_data <- normalizeGiotto(gobject = sc_data, scalefactor = 6000, verbose = T)
  sc_data <- calculateHVG(gobject = sc_data)
  gene_metadata = fDataDT(sc_data) 
  sc_featgenes = gene_metadata[hvg == 'yes']$gene_ID 

  sc_data <- runPCA(gobject = sc_data, genes_to_use = sc_featgenes, scale_unit = F)

  #calculate Sig for deconvolution, This step use DEG function implemented in Giotto
  sc_data@cell_metadata$leiden_clus <- sc_meta[,celltype_key]  # change this
  scran_markers_subclusters = findMarkers_one_vs_all(gobject = sc_data,
                                                      method = 'gini',logFC=0.25,pval=0.05,
                                                      expression_values = 'normalized',
                                                      cluster_column = 'leiden_clus')
  Sig_scran <- unique(scran_markers_subclusters$genes[which(scran_markers_subclusters$comb_rank <= 1000)])#1155 signature gene

  #Calculate median expression value of signature genes in each cell type
  norm_exp <- 2^(sc_data@norm_expr)-1
  id <- sc_data@cell_metadata$leiden_clus 
  ExprSubset <- norm_exp[Sig_scran,] 
  Sig_exp <- NULL 
  for (i in unique(id)){
  Sig_exp <- cbind(Sig_exp,(apply(ExprSubset,1,function(y) mean(y[which(id==i)]))))
  }
  colnames(Sig_exp) <- unique(id) 

  for (seed in seeds) {
    for (cell in cells) {
      for (noise in noises) {
        
        print(paste("seed:", seed, "_cells:", cell, "_noise", as.integer(noise * 100)))
              
        output_file_path <- paste(output_dir, "/", method,"_reference_noise",as.integer(reference_noise * 100), "_seed", seed, "_cells", cell, "_noise", as.integer(noise * 100), sep = "")
        bulk_file_path <- paste(sp_data_dir, "/Simu_seed", seed, "_cells", cell, "_noise", as.integer(noise * 100), "_expression.csv", sep = "")
        sp_loc_path <- paste(sp_data_dir, "/Simu_seed", seed, "_cells", cell, "_noise", as.integer(noise * 100), "_loc.csv", sep = "")
      
        sp_loc = read.csv(sp_loc_path,header=T,row.names=1)
        tryCatch({
          spatial_count <- read.csv(bulk_file_path,header=T,row.names=1)                    
          names(spatial_count) = gsub("\\.", "_", names(spatial_count))

          run_SpatialDWLS(sign_matrix=Sig_exp, spatial_locs=sp_loc,spatialRNA=spatial_count, output_file_path=output_file_path)

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
