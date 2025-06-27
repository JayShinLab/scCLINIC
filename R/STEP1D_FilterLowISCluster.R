
#' STEP1D Filter Low IS Cluster
#' @param Name Name of Run
#' @param obj Seurat Object return by Step1C
#' @param Outdir Output Directory of scCLINIC Results
#' @param resol Step1A: FindClusters resolution for initial clustering in Step1, default using Seurat FindClusters default resolution: 0.8
#' @param OverlapRatio Step1B: The overlap ratio to carry forward to Step1C, default using 0.5, decrease to have coarse resolution. User manual annotation metadata colname, if CELLANNOTATION = TRUE,
#' @param ISThreshold Step1C: The threshold for Indentity Score (IS) to classify low quality cells
#' @param CELLANNOTATION default is FALSE, TRUE if using user cluster annotation
STEP1D_FilterLowISCluster <- function(obj, Outdir, Name, resol = 0.8, OverlapRatio = 0.5, ISThreshold = 0,CELLANNOTATION = FALSE, SCT = FALSE){
  #Filter low IS cluster, based on ISThreshold, by default = 0 (IS -ve is low quality cells)
  message("Step1D started.")
  folder_path_Step1 <- paste0(Outdir,Name,"_Step1/")

  if (CELLANNOTATION){
    message("Using user-annotated clusters.")
    OverlapRatio <- "annotation_index" #User manual cellannotation
  }else{
    OverlapRatio <- paste0("Overlap_Ratio_",OverlapRatio)
    message(paste0("Using ",OverlapRatio))
  }

  message("Filtering low quality cell clusters.")
  lowqualityqc <- read.table(paste0(folder_path_Step1,"1c_","resol_",resol,"_",OverlapRatio,"_IdentityScore.csv"), sep=",", header = T)
  lowqualityqc$cluster <- as.character(lowqualityqc$cluster)
  goodqualitycluster <- lowqualityqc[lowqualityqc$IdentityScore > ISThreshold,]$cluster

  expr <- FetchData(obj, vars = OverlapRatio)
  subGlobal <- obj[,which(expr[,1] %in% goodqualitycluster)] #subset ONLY those good quality major clusters
  p1 <- DimPlot(subGlobal, reduction = "umap",group.by = OverlapRatio,raster=FALSE, pt.size = 0.25)
  p2 <- FeaturePlot(subGlobal,features = "nCount_RNA",reduction = 'umap')
  ggsave(filename = paste0(folder_path_Step1,"1d_","resol_",resol,"_",OverlapRatio,".png"), p1+p2, height = 10, width = 20, dpi = 300)

  #Save RDS
  # if (ncol(subGlobal) < ncol(obj)){
  #   saveRDS(subGlobal,file = paste0(folder_path_Step1,"1d_","resol_",resol,"_",OverlapRatio,".rds"))
  # }else{
  #   if (file_exists(paste0(folder_path_Step1,"1b_","resol_",resol,".rds"))){
  #     createLink(link=paste0(folder_path_Step1,"1d_","resol_",resol,"_",OverlapRatio,".rds"), target=paste0(folder_path_Step1,"1b_","resol_",resol,".rds"), overwrite=TRUE)
  #   }
  #   }
  saveRDS(subGlobal,file = paste0(folder_path_Step1,"1d_","resol_",resol,"_",OverlapRatio,".rds"))

  #Calculate Global Marker after removing low quality cells
  if (length(goodqualitycluster) != length(lowqualityqc$cluster)){
    Idents(object = subGlobal) <- OverlapRatio

    if(inherits(subGlobal[["RNA"]]$counts,what = "IterableMatrix")){ #check for BPcells object, removes temp transpose matrices to prevent IO error
    Global.markers <- data.frame()
    for (i in levels(Idents(subGlobal))) {
      marker<- FindMarkers(subGlobal, ident.1 = i, min.pct = 0.25, logfc.threshold = 0.25, verbose = verbose) # Run FindMarkers one cluster at a time
      marker$cluster = i
      marker$gene = rownames(marker)
      Global.markers <- rbind(Global.markers,marker)
      unlink(list.files(tempdir(), full.names = TRUE), recursive = TRUE) #Remove tmp files before running next cluster
      }
    Global.markers <- subset(Global.markers, p_val < 0.01) # return.thresh < 0.01
    rownames(Global.markers) <- make.unique(names = as.character(x = Global.markers$gene))
    } else {
    Global.markers <- FindAllMarkers(subGlobal, min.pct = 0.25, logfc.threshold = 0.25)}
  }else{
    Global.markers <- read.table(paste0(folder_path_Step1,"1c_","resol_",resol,"_",OverlapRatio,"_GlobalMarker.csv"), sep = ",", header = T, row.names = 1)
  }
  write.table(Global.markers,file = paste0(folder_path_Step1,"1d_","resol_",resol,"_",OverlapRatio,"_GlobalMarker.csv"),sep = ",")
  return(subGlobal)
  message("Step1D completed.")
}
