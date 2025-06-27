
#' STEP2A Subcluster Major Cluster
#' @param Name Name of Run
#' @param obj Seurat Object return by Step1D
#' @param Outdir Output Directory of scCLINIC Results
#' @param resol Step1A: FindClusters resolution for initial clustering in Step1, default using Seurat FindClusters default resolution: 0.8
#' @param OverlapRatio Step1B: The overlap ratio to carry forward to Step1C, default using 0.5, decrease to have coarse resolution. User manual annotation metadata colname, if CELLANNOTATION = TRUE,
#' @param gene_n Number of global markers to used
#' @param CELLANNOTATION default is FALSE, TRUE if using user cluster annotation
STEP2A_Subcluster <- function(obj, Outdir, Name, resol = 0.8, OverlapRatio = 0.5, gene_n = 150, CELLANNOTATION = FALSE, verbose = FALSE, SCT = FALSE){
  ###No Return
  message("Step2A started.")
  folder_path_Step1 <- paste0(Outdir,Name,"_Step1/")

  ##Step 2###################################################################Creating Output Folders
  if (CELLANNOTATION){
    message("Using user-annotated clusters.")
    OverlapRatio <- "annotation_index" #User manual cellannotation
  }else{
    OverlapRatio <- paste0("Overlap_Ratio_",OverlapRatio)
    message(paste0("Using ",OverlapRatio))
  }

  folder_path_Step2 <- paste0(Outdir,Name,"_Step2/")
  create_folder_if_not_exists(folder_path_Step2)

  folder_path_Step2_L1R <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/")
  create_folder_if_not_exists(folder_path_Step2_L1R)

  folder_path_Step2_L1R_Marker <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/Marker/")
  create_folder_if_not_exists(folder_path_Step2_L1R_Marker)

  folder_path_Step2_Output <- paste0(Outdir,Name,"_Step2/Output_",OverlapRatio,"/")
  create_folder_if_not_exists(folder_path_Step2_Output)
  ###############################################################################
  ##Subset & Recluster - selected markers for the variable features during reclustering, positive markers with pct1 > pct2, pct1 > 50%, pct2 < 66.66%, top 150 markers (each major clusters) based on avglog2FC
  Global.markers <- read.table(paste0(folder_path_Step1,"1d_","resol_",resol,"_",OverlapRatio,"_GlobalMarker.csv"), sep = ",", header = T, row.names = 1)
  #Global.markers <- Global.markers[Global.markers$p_val_adj <= 0.05,]
  #Select global markers for reclustering (top150 marker per clusters)
  Global.markers <- Global.markers[Global.markers$avg_log2FC > 0,]
  Global.markers <- Global.markers[Global.markers$pct.1 > Global.markers$pct.2,]
  Global.markers$ES <- Global.markers$avg_log2FC*(Global.markers$pct.1-Global.markers$pct.2)
  Global.markers <- Global.markers[order(Global.markers$ES, decreasing = TRUE), ]

  Global.markers %>%
    group_by(cluster) %>%
    top_n(n = gene_n, wt = ES) %>%
    ungroup()-> Global_Top20 #select top marker for each major clusters based on ES

  qcsclst <- sort(unique(obj@meta.data[,OverlapRatio]))

  #Subset each major cluster and recluster based on selected global markers as variable features
  for (i in sort(qcsclst)){
    message(paste0("Subclustering ",i))
    expr <- FetchData(object = obj, vars = OverlapRatio)
    cell_recluster <- obj[, which(x = expr == i)]#Subset major cluster
    cell_recluster <- DietSeurat(cell_recluster, counts = TRUE, data = TRUE, scale.data = FALSE)
    tryCatch({
      VariableFeatures(cell_recluster) <- Global_Top20$gene#Variable Features = selected top global markers
      cell_recluster <- standard_seurat_clustering(cell_recluster,res=1.5,VF=FALSE,Verbose=verbose, SCT = SCT)##Recluster with seurat resolution 1.5, set VF to FALSE
      cell_recluster$scCLINIC_subcluster <- paste0("S",cell_recluster$seurat_clusters)

      p3<- DimPlot(cell_recluster, reduction = "umap",group.by = "scCLINIC_subcluster",raster=FALSE,pt.size = 0.1) + ggtitle(paste0(i," subclusters"))
      ggsave(paste0(folder_path_Step2_L1R,OverlapRatio,"_cluster_",i,".png"),p3, height = 5, width = 5, dpi = 300)

      saveRDS(cell_recluster,file = paste0(folder_path_Step2_L1R,OverlapRatio,"_cluster_",i,".rds"))

      ###

    }, error = function(e) {
      warning(paste("scCLINIC cannot subcluster cluster ", i, ": ", conditionMessage(e)))#Print errors if major cluster not able to re clustered, while code proceed
    })
  }


  cluster_to_consider <- list()#The list of major clusters which had been successfully re clustered
  for (i in sort(qcsclst)){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".rds")
    if (file_name %in% list.files(folder_path_Step2_L1R)){
      cluster_to_consider <- unlist(c(cluster_to_consider,i))
    }
  }

  #Calculate the avglogFC, pct.1 and pct.2 of the top150 global marker genes for each subclusters using FindAllMarkers
  #Top 150 global markers (each major clusters) selection based on avglog2FC*(pct.1-pct.2)
  n_topgene <- gene_n#################################ONLY DEPENDENT
  Global.markers <- read.table( paste0(folder_path_Step1,"1d_","resol_",resol,"_",OverlapRatio,"_GlobalMarker.csv"), sep = ",", header = T, row.names = 1)
  #Global.markers <- Global.markers[Global.markers$p_val_adj <= 0.05,]
  Global.markers <- Global.markers[Global.markers$avg_log2FC > 0,]
  Global.markers <- Global.markers[Global.markers$pct.1 > Global.markers$pct.2,]
  Global.markers$ES <- Global.markers$avg_log2FC*(Global.markers$pct.1-Global.markers$pct.2)
  Global.markers <- Global.markers[order(Global.markers$ES, decreasing = TRUE), ]
  Global.markers <- Global.markers[Global.markers$cluster %in% cluster_to_consider,]
  Global.markers %>%
    group_by(cluster) %>%
    top_n(n = n_topgene, wt = ES) %>%
    ungroup()-> Global_Top20 #select top marker for each major clusters based on ES

  #For each major cluster, FindAllMarkers within the top150 global markers for each subclusters
  for (i in sort(qcsclst)){
    message(paste0("Finding markers for ",i))
    file_name <- paste0(OverlapRatio,"_cluster_",i,".rds")
    if (file_name %in% list.files(folder_path_Step2_L1R)){
      cell_recluster <- readRDS(paste0(folder_path_Step2_L1R,file_name))
      Idents(object = cell_recluster) <- "scCLINIC_subcluster"

      if(inherits(cell_recluster[["RNA"]]$counts,what = "IterableMatrix")){ #check for BPcells object, removes temp transpose matrices to prevent IO error
        sc.markers <- data.frame()
        for (i in levels(Idents(cell_recluster))) {
          marker<- FindMarkers(cell_recluster, features = unique(Global_Top20$gene),ident.1 = i,
                               min.cells.group = 0,
                               min.cells.feature = 0,
                               min.pct = 0,
                               logfc.threshold = 0, 
                               return.thresh = Inf,
                               only.pos = FALSE) #findallmarker only calculate all top markers of major clusters. NOTE: These parameter values calculate all markers input in "features"
          marker$cluster = i
          marker$gene = rownames(marker)
          sc.markers <- rbind(sc.markers,marker)
          unlink(list.files(tempdir(), full.names = TRUE), recursive = TRUE) #Remove tmp files before running next cluster
        }
        rownames(sc.markers) <- make.unique(names = as.character(x = sc.markers$gene))
      } else {
      sc.markers <- FindAllMarkers(cell_recluster, features = unique(Global_Top20$gene),
                                   min.cells.group = 0,
                                   min.cells.feature = 0,
                                   min.pct = 0,
                                   logfc.threshold = 0,
                                   return.thresh = Inf,
                                   only.pos = FALSE)#findallmarker only calculate all top markers of major clusters. NOTE: These parameter values calculate all markers input in "features"
      }
      message(paste0("Major Cluster Markers: ",length(unique(Global_Top20$gene))))
      message(paste0("Subcluster ", i," Markers: ",nrow(sc.markers)))
      write.table(sc.markers,file = paste0(folder_path_Step2_L1R_Marker,sub("\\.rds$", "", file_name),".csv"),sep=",")
    }
  }
  message("Step2A completed.")
}
