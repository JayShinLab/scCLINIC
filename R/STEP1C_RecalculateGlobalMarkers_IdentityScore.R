#Step1c-recalculate markers and calculate Identity Score (IS)
#Find marker for scCLINIC resolution
#' Step1c-recalculate markers and calculate Identity Score (IS)
#' @param Name Name of Run
#' @param obj Seurat Object return by Step1B
#' @param Outdir Output Directory of scCLINIC Results
#' @param resol Step1A: FindClusters resolution for initial clustering in Step1, default using Seurat FindClusters default resolution: 0.8
#' @param OverlapRatio Step1B: The overlap ratio to carry forward to Step1C, default using 0.5, decrease to have coarse resolution. User manual annotation metadata colname, if CELLANNOTATION = TRUE,
#' @param gene_n Number of global markers to used
#' @param CELLANNOTATION default is FALSE, TRUE if using user cluster annotation
STEP1C_RecalculateGlobalMarkers_IdentityScore <- function(obj, Outdir, Name, resol = 0.8, OverlapRatio = 0.5, gene_n = 150, CELLANNOTATION = FALSE, verbose = TRUE, SCT = FALSE){
  message("Step1C started.")
  folder_path_Step1 <- paste0(Outdir,Name,"_Step1/")
  create_folder_if_not_exists(folder_path_Step1)

  if (CELLANNOTATION){
    message("Using user-annotated clusters.")
    obj@meta.data$annotation_index <- as.numeric(factor(obj@meta.data[,OverlapRatio]))#change the manual annotation to index, eg. CellType0 CellType1 CellType2 CellType3 to 1 2 3 4
    OverlapRatio <- "annotation_index" #User manual cellannotation
    # Check if object has UMAP
    if (!("umap" %in% names(obj@reductions))){
          if (!("pca" %in% names(obj@reductions))){
              if (is.null(obj[[DefaultAssay(obj)]]["scale.data"])){
                  obj <- standard_seurat_clustering(obj, Verbose = verbose, SCT = SCT)
              } else {
                  obj <- RunPCA(obj)
                  obj <- RunUMAP(obj,dims=1:30)
              }
              } else {
                  obj <- RunUMAP(obj,dims=1:30)
              }
          }
  }else{
    OverlapRatio <- paste0("Overlap_Ratio_",OverlapRatio)
    message(paste0("Using ",OverlapRatio))
  }

  clus.names <- obj@meta.data[[OverlapRatio]]
  if (any(!(grepl("^[a-zA-Z]|^\\.[^0-9]", clus.names)))) {
    clus.names <- ifelse(
      !(grepl("^[a-zA-Z]|^\\.[^0-9]", clus.names)),
      paste0("M", clus.names),
      clus.names
    )
    obj@meta.data[[OverlapRatio]] <- clus.names
    obj@meta.data[[OverlapRatio]] <- factor(obj@meta.data[[OverlapRatio]]) #Preserve levels
    message = paste0("Appending `M` to cluster names to ensure valid variable names")
  }
  Idents(object = obj) <- OverlapRatio
  message("Finding markers.")

  if(inherits(obj[["RNA"]]$counts,what = "IterableMatrix")){ #check for BPcells object, removes temp transpose matrices to prevent IO error
    Global.markers <- data.frame()
    for (i in levels(Idents(obj))) {
      message(paste0("Cluster ",i))
      marker<- FindMarkers(obj, ident.1 = i, min.pct = 0.25, logfc.threshold = 0.25, verbose = verbose) # Run FindMarkers one cluster at a time
      marker <- subset(marker, p_val < 0.01) # return.thresh < 0.01
      if (nrow(marker) == 0) {
        stop(paste0("Error: No significant markers (p_val < 0.01) found for cluster ", i, ".")) # Check for clusters with no marker genes
        }
      marker$cluster = i
      marker$gene = rownames(marker)
      Global.markers <- rbind(Global.markers,marker)
      unlink(list.files(tempdir(), full.names = TRUE), recursive = TRUE) #Remove tmp files before running next cluster
      }
    rownames(Global.markers) <- make.unique(names = as.character(x = Global.markers$gene))
    } else {
    Global.markers <- FindAllMarkers(obj, min.pct = 0.25, logfc.threshold = 0.25)
    clusters_with_no_markers <- setdiff(levels(Idents(obj)), unique(Global.markers$cluster))
    if (length(clusters_with_no_markers) > 0) {
        stop(paste0("Error: No significant markers (p_val < 0.01) found for the following cluster(s): ", # Check for clusters with no marker genes
              paste(clusters_with_no_markers, collapse = ", "), "."))
    }}
  
  write.table(Global.markers,file = paste0(folder_path_Step1,"1c_","resol_",resol,"_",OverlapRatio,"_GlobalMarker.csv"),sep = ",")
  
  obj$all <- 1
  avegloballst <- AverageExpression(obj,group.by = OverlapRatio,slot = "counts")#average expression of each scCLINIC major clusters
  aveallcell <- AverageExpression(obj,group.by = "all",slot = "counts")
  ncell_table <- table(obj@meta.data[,OverlapRatio])

  message("Calculating Identity Score.")
  #Global.markers <- Global.markers[Global.markers$p_val_adj <= 0.05,] note: by default not using pval_adj to select top markers
  Global.markers <- Global.markers[Global.markers$avg_log2FC > 0,]
  Global.markers$ES <- Global.markers$avg_log2FC*(Global.markers$pct.1-Global.markers$pct.2)
  Global.markers$DPE <- Global.markers$pct.1-Global.markers$pct.2
  Global.markers <- Global.markers[order(Global.markers$ES, decreasing = TRUE), ]

  Global.markers %>%
    group_by(cluster) %>%
    top_n(n = gene_n, wt = ES) %>%
    ungroup()-> Global_Top20 #Select Top N global markers for each major clusters based on avg_log2FC*(pct.1-pct.2)

  Global.markers <- Global_Top20

  #Extract information for the Top N global markers for each major clusters, including average expression, max expression, min expression, number of cells,
  #Note: local = major cluster x, ref = all cells except major cluster x, all = all cells
  Global.markers$aveexp_local <- 0
  Global.markers$aveexp_ref <- 0
  Global.markers$maxexp_ref <- 0
  Global.markers$minexp_ref <- 0
  Global.markers$ncell_local <- 0
  Global.markers$aveexp_all <- 0
  for (row in seq(nrow(Global.markers))){
    Global.markers[row,]$aveexp_local <- avegloballst$RNA[Global.markers[row,]$gene,as.character(Global.markers[row,]$cluster)]
    refexp <- avegloballst$RNA[Global.markers[row,]$gene,]#each global markers' average expression of each scCLINIC major clusters
    refexp <- refexp[-which(names(refexp) == as.character(Global.markers[row,]$cluster))]
    Global.markers[row,]$ncell_local <- ncell_table[as.character(Global.markers[row,]$cluster)]
    Global.markers[row,]$aveexp_ref <- mean(refexp)
    Global.markers[row,]$maxexp_ref <- max(refexp)
    Global.markers[row,]$minexp_ref <- min(refexp)
    Global.markers[row,]$aveexp_all <- aveallcell$RNA[Global.markers[row,]$gene,]
  }

  #Summarized the global marker information based on scCLINIC major clusters
  Global.markers$count <- 1
  lowqualityqc <- Global.markers %>%
    group_by(cluster) %>%
    summarize(
      IdentityScore = mean(DPE),
      No_Gene = sum(count),
      No_Cell = mean(ncell_local),
      avg_log2FC = mean(avg_log2FC),
      pct.1 = mean(pct.1),
      pct.2 = mean(pct.2),
      Avg_Exp = mean(aveexp_local),
      Avg_Exp_Others = mean(aveexp_ref)
    ) %>%
    ungroup()

  write.table(lowqualityqc,file = paste0(folder_path_Step1,"1c_","resol_",resol,"_",OverlapRatio,"_IdentityScore.csv"),sep=",",row.names = FALSE)
  return(obj)
  message("Step1C completed.")
}
