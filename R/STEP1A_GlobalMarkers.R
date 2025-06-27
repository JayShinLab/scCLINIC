######################################################################################################################################################STEP1A_GlobalMarkers
## Step 1A-calculate global markers (before merging)
#' Step 1A-calculate global markers (before merging)
#' @param Name Name of Run
#' @param Input Input Directory of Seurat Object, or seurat object
#' @param Outdir Output Directory of scCLINIC Results
#' @param resol Step1A: FindClusters resolution for initial clustering in Step1, default using Seurat FindClusters default resolution: 0.8
STEP1A_GlobalMarkers <- function(Input,Outdir,Name,resol=0.8,verbose = TRUE, SCT = FALSE){
  message("Step1A started.")
  if (is.character(Input) && grepl("\\.RDS$", Input, ignore.case = TRUE)) {
    # Load the .RDS file
    obj <- readRDS(Input)
    } else {
    obj <- Input
  }

  if(SCT){
    res <- paste0("SCT_snn_res.",resol)
  } else {
    res <- paste0("RNA_snn_res.",resol)
  }
  
  obj <- obj[["RNA"]]$counts #only extract raw counts
  obj <- CreateSeuratObject(obj)
  obj$nCount_RNA <- SeuratObject::.CalcN(obj[["RNA"]]$counts)$nCount
  obj$nFeature_RNA <- SeuratObject::.CalcN(obj[["RNA"]]$counts)$nFeature
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-") #percentage Mito

  folder_path_Step1 <- paste0(Outdir,Name,"_Step1/")
  create_folder_if_not_exists(folder_path_Step1)

  message("Clustering cells.")
  obj <- standard_seurat_clustering(obj, SCT = SCT, Verbose = verbose)
  obj <- FindClusters(obj, resolution = resol, n.start = 10) #Initial clustering resolution based on resol parameter
  #obj@meta.data[[res]] <- paste0("M",obj@meta.data[[res]])
  p1 <- DimPlot(obj, reduction = "umap",group.by = res, raster=FALSE, label = T, repel = T)
  p2 <- FeaturePlot(obj,features = "nCount_RNA",reduction = 'umap')
  ggsave(filename = paste0(folder_path_Step1,"1a_resol_",resol,".png"), p1+p2, height = 10, width = 20, dpi = 300)
  message("Clustering completed.")
  
  #Calculate global markers
  message("Finding markers")
  Idents(object = obj) <- res

  if(inherits(obj[["RNA"]]$counts,what = "IterableMatrix")){ #check for BPcells object, removes temp transpose matrices to prevent IO error
    Global.markers <- data.frame()
    for (i in levels(Idents(obj))) {
      message(paste0("Cluster ",i))
      marker<- FindMarkers(obj, ident.1 = i, min.pct = 0.25, logfc.threshold = 0.25, verbose = verbose) # Run FindMarkers one cluster at a time
      marker$cluster = i
      marker$gene = rownames(marker)
      Global.markers <- rbind(Global.markers,marker)
      unlink(list.files(tempdir(), full.names = TRUE), recursive = TRUE) #Remove tmp files before running next cluster
      }
    Global.markers <- subset(Global.markers, p_val < 0.01) # return.thresh < 0.01
    rownames(Global.markers) <- make.unique(names = as.character(x = Global.markers$gene))
    } else {
    Global.markers <- FindAllMarkers(obj, min.pct = 0.25, logfc.threshold = 0.25, verbose = verbose)
    }
  
  write.table(Global.markers,file = paste0(folder_path_Step1,"1a_resol_",resol,"_Global_Marker.csv"), sep = ",")
  message("Step1A completed.")
  return(obj)
}
