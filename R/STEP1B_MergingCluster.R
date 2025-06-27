#Step1b-merging cluster
#' Step1b-merging cluster
#' @param Name Name of Run
#' @param obj Seurat Object return by Step1A
#' @param Outdir Output Directory of scCLINIC Results
#' @param resol Step1A: FindClusters resolution for initial clustering in Step1, default using Seurat FindClusters default resolution: 0.8
#' @param overlapRatioList Step1B: The overlap ratio to try in Step1B (recommended range: 0.1 to 0.5)
#' @param OverlapRatio Step1B: The overlap ratio to carry forward to Step1C, default using 0.5, decrease to have coarse resolution. User manual annotation metadata colname, if CELLANNOTATION = TRUE,
#' @param gene_n Number of global markers to used
STEP1B_MergingCluster <- function(obj, Outdir, Name, resol = 0.8, overlapRatioList = c(0.1,0.25,0.5,0.75,0.9), gene_n = 150, SCT = FALSE){
  message("Step1B started.")
  folder_path_Step1 <- paste0(Outdir,Name,"_Step1/")

  if(SCT){
    res <- paste0("SCT_snn_res.",resol)
  } else {
    res <- paste0("RNA_snn_res.",resol)
  }

  #Global.markers <- readRDS(paste0(folder_path_Step1,"1a_resol_",resol,"_Global_Marker.rds"))
  Global.markers <- read.table(paste0(folder_path_Step1,"1a_resol_",resol,"_Global_Marker.csv"), sep= ",", header = T, row.names = 1)
  Global.markers <- Global.markers[Global.markers$avg_log2FC > 0,]
  Global.markers$ES <- Global.markers$avg_log2FC*(Global.markers$pct.1-Global.markers$pct.2)
  Global.markers <- Global.markers[order(Global.markers$ES, decreasing = TRUE), ]
  for (overlapratio in overlapRatioList){
    message(paste0("Merging clusters for Overlap Ratio ",overlapratio))
    Global.markers %>%
      group_by(cluster) %>%
      top_n(n = gene_n, wt = ES) %>%
      ungroup()-> Global_Top20 #Select Top N global markers for each major clusters based on avg_log2FC*(pct.1-pct.2)

    #Find the overlapped global markers between major clusters
    duplicates <- Global_Top20[duplicated(Global_Top20$gene) | duplicated(Global_Top20$gene, fromLast = TRUE), ]

    result <- duplicates %>%
      group_by(gene) %>%
      summarize(
        p_val = toString(p_val),
        avg_log2FC = toString(avg_log2FC),
        pct.1 = toString(pct.1),
        pct.2 = toString(pct.2),
        p_val_adj = toString(p_val_adj),
        cluster = toString(cluster),
        ES = toString(ES)
      ) %>%
      ungroup()

    # Count the occurrences of the overlapped global markers in each major clusters
    cluster_pair_counts <- duplicates %>%
      group_by(cluster) %>%
      summarise(count = n()) %>%
      filter(count > 0)
    #"result" is the information of each overlapped global markers
    result$cluster <- strsplit(result$cluster, ', ')
    result$cluster <- lapply(result$cluster, as.integer)
    overla <- list()
    uniqla <- list()

    # Count the occurrences of the overlapped global markers in each major clusters
    for (i in cluster_pair_counts$cluster){
      overlaps <- list()
      for (j in result$cluster){#for each overlapped global markers (j) in the major cluster (i), append all source of artifacts in "overlaps"
        if (i %in% j){
          overlaps <- c(overlaps, j)
        }
      }
      ele <- unlist(overlaps)
      ele <- ele[ele != i]
      overla <- c(overla,list(ele))
      element_counts <- table(ele)# "element_counts" records each source of artifacts occurrences of major cluster (i) overlapped global markers
      # EXAMPLE
      # element_counts for one of the major cluster (i)
      # 0  1  2  3  4  5  6  7  8  9 10 11 13 14 15 17 18 19 20 21 22 23 -> source of artifacts
      # 1  2 31  2 30  1  1 18  1  3  2  8  3  1 15  7  8  2  8 24 11  4 -> occurrences
      elements_appearing_four_or_more <- as.numeric(names(element_counts[element_counts >= overlapratio*nrow(Global_Top20[Global_Top20$cluster == i,])]))#Check which source of artifacts have occurrences higher than overlapratio*(no. of top gene_n global markers for major cluster (i))
      uniqla <- c(uniqla,list(unique(elements_appearing_four_or_more)))
    }
    cluster_pair_counts$overlap <- overla #all the occurrences of respective source of artifacts of each major clusters
    cluster_pair_counts$unioverlap <- uniqla #the source of artifacts of each major clusters which higher than overlapratio*(no. of top gene_n global markers for major cluster)

    clusterpair <- cluster_pair_counts
    clusterpair$cluster <- as.matrix(clusterpair)[,"cluster"]
    clusterpair$cluster <- as.integer(clusterpair$cluster)
    clusterpair$network <- lapply(1:nrow(clusterpair), function(x) list())
    for (row in seq(nrow(clusterpair))){
      clusterpair[row,]$network <- list(c(clusterpair[row,]$cluster, unlist(clusterpair[row,]$unioverlap)))
    }

    clusterpair$overlap <- as.character(cluster_pair_counts$overlap)
    clusterpair$unioverlap <- as.character(cluster_pair_counts$unioverlap)

    OverlapRatio <- paste0("Overlap_Ratio_",overlapratio)
    #write.table(clusterpair[,1:4],paste0(folder_path_Step1,"1b_OverlapInfo_","resol_",resol,"_",OverlapRatio,".csv"),sep = ",")

    #Merging major clusters which higher than overlapratio*(no. of top gene_n global markers for major cluster)
    #Eg. 3, 5, 11 major clusters overlapped, 1, 11 major cluster overlapped, -> network = 1, 3, 5, 11 major clusters overlapped
    #https://stackoverflow.com/questions/47322126/merging-list-with-common-elements
    # network<- unique(sapply(clusterpair$network, function(x)
    #   unique(unlist(clusterpair$network[sapply(clusterpair$network, function(y)
    #     any(x %in% y))])))) -> restructure...because terminal fail to run

    network <- unique(sapply(clusterpair$network, function(x) {
      unique(unlist(clusterpair$network[sapply(clusterpair$network, function(y) {
        any(x %in% y)
      })]))
    }))

    network <- lapply(network, sort)

    #Replace seurat resolution with scCLINIC resolution (merging resolution): eg. 1, 3, 5, 11 will now only become one cluster (Major Cluster 1)
    obj$HighRes <- obj@meta.data[res]
    obj@meta.data[,OverlapRatio] <- obj$HighRes
    for (net in network){
      obj@meta.data[rownames(obj@meta.data[obj@meta.data$HighRes %in% net, ]),][,OverlapRatio] <- net[1]#took the first major cluster in the merging list to represent the list of merged major clusters
    }

    #Reassign major cluster ID
    overlap_table <- table(obj@meta.data[OverlapRatio])
    # Remove the entries with a count of zero
    non_zero_overlap <- overlap_table[overlap_table != 0]
    # Get the original values and their counts
    original_values <- as.numeric(names(non_zero_overlap))
    counts <- as.numeric(non_zero_overlap)
    # Create a mapping from original values to new values starting from 1
    new_values <- seq_along(original_values)
    # Create a named vector for the mapping
    value_mapping <- setNames(new_values, original_values)
    # Reassign the major cluster ID based on the mapping
    obj@meta.data[,OverlapRatio] <- value_mapping[as.character(obj@meta.data[,OverlapRatio])]

    p1 <- DimPlot(obj, reduction = "umap",group.by = OverlapRatio,raster=FALSE, pt.size = 0.25, sizes.highlight = 0.25)
    p2 <- VlnPlot(obj,features = "nCount_RNA",group.by = OverlapRatio)
    ggsave(filename = paste0(folder_path_Step1,"1b_OverlapInfo_","resol_",resol,"_",OverlapRatio,".png"), p1+p2, height = 10, width = 20, dpi = 300)
    message(paste0("Overlap Ratio ",overlapratio," merging completed."))
  }

  saveRDS(obj, paste0(folder_path_Step1,"1b_","resol_",resol,".rds"))
  return(obj)
  message("Step1b completed.")
}
