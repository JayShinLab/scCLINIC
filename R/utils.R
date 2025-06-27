########################################################################################################################################################Global Function
## Standard Seurat Clustering
#' scCLINIC subcluster clustering function
#'
#' @param obj seurat object
#' @param res resolution
#' @param npc npc, default is 30 (Seurat Default)
#' @param VF find variable features, TRUE = enable, FALSE = disable
#' @keywords scCLINIC subcluster clustering function
#' @export
#' @examples
#' standard_seurat_clustering(obj)
#'
standard_seurat_clustering <- function(obj, res = 0.1, npc = 30,VF = TRUE, Verbose = FALSE, SCT = FALSE, Global_Top20 = Global_Top20){
  message("Clustering started.")
  if (ncol(obj) < npc){
    npc = ncol(obj)
  }

  if (SCT){
    if(VF) {
      obj <- SCTransform(obj, vst.flavor = "v2", conserve.memory = TRUE, verbose = Verbose)
    } else {
      obj <- SCTransform(obj, vst.flavor = "v2", conserve.memory = TRUE, verbose = Verbose, min_cells = 0)
      VariableFeatures(obj) <- Global_Top20$gene
    }
    } else {
      obj <- NormalizeData(obj, verbose = Verbose)
      if (VF){
        obj <- FindVariableFeatures(obj, verbose = Verbose)
      }
    obj <- ScaleData(obj, verbose = Verbose)
    }
  obj <- RunPCA(obj,npcs = npc, verbose = Verbose)
  #remove scale.data slot to save memory
  if  (SCT) {
    obj[["SCT"]]$scale.data <- NULL
    } else {
    obj[["RNA"]]$scale.data <- NULL
    }
  obj <- FindNeighbors(obj, dims = 1:npc, reduction = "pca", verbose = Verbose)
  obj <- FindClusters(obj, resolution = as.numeric(res), verbose = Verbose)
  obj <- RunUMAP(obj, dims = 1:npc, verbose = Verbose)
  message("Clustering completed.")
  return(obj)
}

## Function to create new folder for output
#' create folder function
#'
#' @param folder_path folder directory to create
#' @export
#' @examples
#' create_folder_if_not_exists(folder_path)
#'
create_folder_if_not_exists <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    dir.create(folder_path, recursive = TRUE)
  }
}

#' Load test data
#'
#' @return The test data
#' @export
load_data <- function() {
  file_path <- system.file("data", "scCLINIC_Demo.rds", package = "scCLINIC")
  readRDS(file_path)
}

rm_tmp <- function() {
    dir_delete(paste0(Output,Name,"_Step1"))
    dir_delete(paste0(Output,Name,"_Step2/",OverlapRatio,"_recluster"))
  }
