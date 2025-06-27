
#' STEP2B Calculate Contamination Score
#' @param Name Name of Run
#' @param obj Seurat Object return by Step1C
#' @param Outdir Output Directory of scCLINIC Results
#' @param resol Step1A: FindClusters resolution for initial clustering in Step1, default using Seurat FindClusters default resolution: 0.8
#' @param OverlapRatio Step1B: The overlap ratio to carry forward to Step1C, default using 0.5, decrease to have coarse resolution. User manual annotation metadata colname, if CELLANNOTATION = TRUE,
#' @param gene_n Number of global markers to used
#' @param CELLANNOTATION default is FALSE, TRUE if using user cluster annotation
STEP2B_ContaminationScore<-function(obj,Outdir,Name,resol=0.8,OverlapRatio=0.5,gene_n=150,CELLANNOTATION = FALSE, SCT = FALSE){
  message("Step2B started.")
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
  folder_path_Step2_L1R <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/")
  folder_path_Step2_L1R_Marker <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/Marker/")
  folder_path_Step2_Output <- paste0(Outdir,Name,"_Step2/Output_",OverlapRatio,"/")

  ###scCLINIC
  totalcells <- ncol(obj)
  avegloballst <- AverageExpression(obj,group.by = OverlapRatio,slot = "counts")#avg expression of each major clusters
  gcell_table <- table(obj@meta.data[,OverlapRatio])#number of cells of each major clusters
  qcsclst <- sort(unique(obj@meta.data[,OverlapRatio]))#list of major clusters seurat ID
  cluster_to_consider <- list()
  for (i in sort(qcsclst)){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".csv")
    if (file_name %in% list.files(folder_path_Step2_L1R_Marker)){
      cluster_to_consider <- unlist(c(cluster_to_consider,i))
    }
  }

  #Load the top150 global markers
  n_topgene <- gene_n#################################ONLY DEPENDENT
  Global.markers <- read.table(paste0(folder_path_Step1,"1d_","resol_",resol,"_",OverlapRatio,"_GlobalMarker.csv"), sep= ",", header = T, row.names = 1)
  #Global.markers <- Global.markers[Global.markers$p_val_adj <= 0.05,]
  Global.markers <- Global.markers[Global.markers$avg_log2FC > 0,]
  Global.markers <- Global.markers[Global.markers$pct.1 > Global.markers$pct.2,]
  Global.markers$ES <- Global.markers$avg_log2FC*(Global.markers$pct.1-Global.markers$pct.2)
  Global.markers <- Global.markers[order(Global.markers$ES, decreasing = TRUE), ]
  Global.markers <- Global.markers[Global.markers$cluster %in% cluster_to_consider,]

  Global.markers %>%
    group_by(cluster) %>%
    top_n(n = n_topgene, wt = ES) %>%
    ungroup()-> Global_Top20

  #Identify and calculate the Identity Score of the Artifacts in each Subclusters and calculate the scCLINIC score for each subclusters
  overlaplst <- list()

  for (i in cluster_to_consider){#For each major cluster (i)
    message(paste0("Calculating Enrichment Score for ",i))
    file_name <- paste0(OverlapRatio,"_cluster_",i,".rds")
    if (paste0(sub("\\.rds$", "", file_name),".csv") %in% list.files(folder_path_Step2_L1R_Marker)){ #Load each subcluster markers
      sc.markers <- read.table(paste0(folder_path_Step2_L1R_Marker,sub("\\.rds$", "", file_name),".csv"), sep = ",", header = T, row.names = 1)
      subc <- readRDS(paste0(folder_path_Step2_L1R,OverlapRatio,"_cluster_",i,".rds"))
      avegenelst <- AverageExpression(object = subc, group.by = "scCLINIC_subcluster",slot = "counts")#Avg marker expression of each subclusters
      ncell_table <- table(subc@meta.data$"scCLINIC_subcluster")#Number of cells of each subclusters

      #Prefilter the findallmarker result, only keep positive markers with pct1 > pct2, ensure no negative ES
      if (length(sc.markers)!=0){
        sc.markers$ES <- sc.markers$avg_log2FC*(sc.markers$pct.1-sc.markers$pct.2) #Calculate ES score for each markers
        sc.markers$ES <- ifelse(sc.markers$avg_log2FC < 0, 0, sc.markers$ES) #Set marker ES to 0, if marker avglog2FC < 0
        sc.markers$ES <- ifelse(sc.markers$pct.1 < sc.markers$pct.2, 0, sc.markers$ES) #Set marker ES to 0, if marker pct.1 - pct.2 < 0

        sc.markers <- sc.markers[order(sc.markers$ES, decreasing = TRUE), ]

        localsclst <- unique(sc.markers$cluster)
        for (j in sort(localsclst)){#For each subcluster (j) in major cluster (i)
          OwnGene <- Global_Top20[Global_Top20$cluster %in% i,]#Global marker which are present in major cluster i
          Global_sub_Top20 <- anti_join(Global_Top20, OwnGene, by = "gene")# Global marker genes which are not present in major cluster i

          sc_sub_Top20 <- sc.markers[sc.markers$cluster %in% j,]# Subcluster marker (j)

          overlap <- sc_sub_Top20[sc_sub_Top20$gene %in% Global_sub_Top20$gene,]#Subcluster marker (j) which also present in other major clusters

          if (nrow(overlap)>0){
            overlap$globalcluster <- i #global cluster = major cluster i
            contam_global <- Global_sub_Top20[Global_sub_Top20$gene %in% sc_sub_Top20$gene,] #extract the global marker information, including avglog2FC, pct1, pct2
            merged_df <- merge(overlap, contam_global, by = "gene", suffixes = c("_local", "_ref")) #merge both information, _local = subclusters marker information, _ref = global marker information
            merged_df$ncells <- ncell_table[[j]] #number of cells in subcluster j
            merged_df$refcells <- gcell_table[[i]] #number of cells in that major cluster where the global marker present
            if (length(merged_df$gene) == 1){
              merged_df$aveexp_local <- avegenelst$RNA[merged_df$gene,as.character(merged_df$cluster_local)] #avg expression of that global marker in that major cluster where the global marker present
              merged_df$aveexp_ref <- avegloballst$RNA[merged_df$gene,as.character(merged_df$cluster_ref)]  #avg expression of that subcluster marker in subcluster j
            }
            else{
              merged_df$aveexp_local <- avegenelst$RNA[merged_df$gene,as.character(j)]
              merged_df$aveexp_ref <- diag(as.matrix(avegloballst$RNA[merged_df$gene,as.character(merged_df$cluster_ref)]))#as.character(merged_df$cluster_ref)
            }
            overlaplst <- rbind(overlaplst,merged_df)
          }
        }
      }}}
  message("Calculating Contamination Scores.")
  #Notes: Subcluster (j) markers which not present in major cluster (i) aka artifacts
  overlaplst.filtered <- overlaplst
  overlaplst.filtered$cluster_ref <- overlaplst.filtered$cluster_ref
  overlaplst.filtered$cluster_local <- overlaplst.filtered$cluster_local
  overlaplst.filtered$globalcluster <- overlaplst.filtered$globalcluster
  overlaplst.filtered$dp1 <- overlaplst.filtered$ES_local  #Identity Score of the artifact aka ES_local

  write.table(overlaplst.filtered,file = paste0(folder_path_Step2_Output,"ArtifactsInfo.csv"),sep = ",")

  referenceoverlap <- overlaplst.filtered[, c("dp1","cluster_local" , "globalcluster", "cluster_ref")]
  overlaplst.filtered <-  overlaplst.filtered[, c("gene","dp1","p_val_local", "avg_log2FC_local", "pct.1_local", "pct.2_local", "p_val_adj_local", "cluster_local" ,   "ES_local", "globalcluster", "ncells", "aveexp_local")]

  #Remove duplicated artifacts contributed by different major clusters, eg. Artifact A present across major cluster 1 2 3, in overlaplst.filtered appeared three times...
  dup_rows <- duplicated(overlaplst.filtered)
  overlaplst.filtered <- overlaplst.filtered[!dup_rows, ]

  #Summarize the artifacts per subclusters, calculate the scCLINIC score for each subclusters
  contamdfall <- overlaplst.filtered %>%
    group_by(cluster_local, globalcluster) %>%
    summarise(count = n(),mean(dp1),mean(ncells))#sum the Identity Score for each artifacts in the subcluster aka Subcluster scCLINIC Score, count the number of unique artifacts in the subcluster

  contamdfall$cluster_local <- as.matrix(contamdfall)[,"cluster_local"]
  #contamdfall$cluster_local <- as.integer(contamdfall$cluster_local)
  #contamdfall$globalcluster <- as.integer(contamdfall$globalcluster)#Summary including the Subclusters ID (globalcluster & cluster_local) and the respective scCLINIC Score (mean(dp1))

  #Classified subclusters to either contaminated or non-contaminated subclusters based on elbow point
  contamdfall$scCLINICScore <- contamdfall$`mean(dp1)`
  contamdfall_backup <- contamdfall
  contamdfall <- contamdfall[order(contamdfall$"scCLINICScore", decreasing = TRUE), ]

  #write.table(contamdfall,file = paste0(folder_path_Step2_Output,"ArtifactsInfo.csv"),sep = ",")

  x1 <- as.numeric(rownames(contamdfall))#Subclusters instances
  y1 <- contamdfall$"scCLINICScore"#scCLINIC Score of respective subclusters
  z1 <- paste0(contamdfall$globalcluster,"_",contamdfall$cluster_local)#Respective subclusters ID, eg. M1S2 = Subcluster 2 in Major Cluster 1
  ncell1 <- contamdfall$`mean(ncells)`#no. of Cells of respective subclusters

  ycopy <- y1
  # Step 3: Calculate the first and second derivatives
  # y1 = y1#sgolayfilt(y1, p = 3, n = 5) #Optional Smoothing
  # first_derivative <- diff(y1) / diff(x1)
  # second_derivative <- diff(first_derivative) / diff(x1[-length(x1)])
  # # Step 5: Identify local maxima of the second derivative
  # local_maxima <- findpeaks(second_derivative, minpeakheight = 0.01, minpeakdistance = 5)
  # maxima_x <- sort(x1[local_maxima[,2]],decreasing = F)
  # maxima_y <- sort(y1[local_maxima[,2]] ,decreasing = T)

  #Maxima
  local_maxima <- findpeaks(diff(diff(y1)), npeaks = 7,minpeakdistance = 5, sortstr = TRUE)
  maxima_x <- sort(local_maxima[,2], decreasing = F)
  maxima_y <- sort(y1[local_maxima[,2]] ,decreasing = T)

  #write.table(local_maxima,file = paste0(folder_path_Step2_Output,"second_derivative.csv"),sep = ",")

  #Optimal maxima
  optimal_maxima <- findpeaks(diff(diff(diff(y1))), npeaks = 7,minpeakdistance = 5, sortstr = TRUE)
  opmaxima_x <- optimal_maxima[,2]+1
  #opmaxima_y <- sort(y1[optimal_maxima[,2]+1] ,decreasing = T)

  # Function to find the closest point in maxima_x for a given value
  find_closest <- function(value, maxima_x) {
    # Compute the absolute differences
    differences <- abs(maxima_x - value)
    # Find the index of the minimum difference
    closest_index <- which.min(differences)
    # Return the closest value
    return(maxima_x[closest_index])
  }

  # Apply the function to each element in opmaxima_x
  closest_points <- sapply(opmaxima_x, find_closest, maxima_x = maxima_x)

  #write.table(optimal_maxima,file = paste0(folder_path_Step2_Output,"third_derivative.csv"),sep = ",")

  #Plot scCLINIC Score vs Percentage Droplets (%)
  y1 <- ycopy #Plot using raw score not smoothing (if smoothing is applied)
  data <- data.frame(x = x1, y = y1,z = z1,ncelllst = cumsum(ncell1)/totalcells*100)#ncelllst = accumulative percentage droplets of each subclusters
  # Calculate the scaling factor
  range_y1 <- range(data$y)
  range_ncell1_accumulated <- range(data$ncelllst)
  scaling_factor <- diff(range_y1) / diff(range_ncell1_accumulated)

  # Base plot
  p <- ggplot(data, aes(x = x1, y = y1, label=z1)) +
    geom_point() +
    geom_line(aes(y = y), col = "#53a0b9") +
    geom_line(aes(y = ncelllst * scaling_factor), col = "#d55046") +
    geom_text(aes(label = z), angle = 90, hjust = 0.5, vjust = -0.5) +
    labs(title = paste0("Contamination Score (no. of cells = ", totalcells, ")"),
         x = "Subclusters (ID)", y = "Contamination Score") +
    scale_y_continuous(
      name = "Contamination Score",
      sec.axis = sec_axis(~./scaling_factor, name = "Percentage Droplets (%)")
    ) +
    theme_minimal() +
    theme(
      axis.line.y = element_line(color = "#53a0b9"),
      axis.text.y = element_text(color = "#53a0b9"),
      axis.line.y.right = element_line(color = "#d55046"),
      axis.text.y.right = element_text(color = "#d55046"),
      axis.line = element_line(color = "black"),
      axis.text.x = element_text(color = "black")
    )+
    annotate("text", x =  x1[length(x1)], y = y1[1], label = length(maxima_x)+1, hjust = 1, col = "black")


  # Create list to store additional layers
  additional_layers <- list()

  nelbowidx <- 0
  # Loop to create vertical and horizontal lines and annotations
  for (nelbow in maxima_x) {
    nelbowidx <- nelbowidx + 1
    nelbow_y <- y1[nelbow]
    additional_layers <- c(
      additional_layers,
      geom_vline(xintercept = nelbow, linetype = "dashed", color = "black"),
      #geom_hline(yintercept = nelbow_y, linetype = "dashed", color = "#53a0b9"),
      #geom_hline(yintercept = data$ncelllst[nelbow] * scaling_factor, linetype = "dashed", color = "#d55046"),
      # geom_segment(x = nelbow, y = data$ncelllst[nelbow] * scaling_factor,
      #                  xend = x1[length(x1)], yend = data$ncelllst[nelbow] * scaling_factor,
      #              linetype = "dashed", color = "#d55046"),
      #annotate("text", x = 0, y = nelbow_y, label = round(nelbow_y, 3), vjust = -1, col = "#53a0b9"),
      annotate("text", x = nelbow, y = data$ncelllst[nelbow] * scaling_factor, label = floor(data$ncelllst[nelbow]*1000)/1000, vjust = 1,hjust=-0.1, col = "#d55046"),
      annotate("text", x = nelbow, y = y1[1], label = nelbowidx, hjust = 1.5, col = "black")
    )
  }

  # Add the additional layers to the plot
  p <- p + additional_layers +
    geom_vline(xintercept = closest_points[1], linetype = "dashed", color = "#d55046")+
    geom_vline(xintercept = closest_points[2], linetype = "dashed", color = "#7c3b36")

  ggsave(filename = paste0(folder_path_Step2_Output,"scCLINICScoreKneePlot",".png"), p, height = 10, width = 10, dpi = 300)

  #Summarize the artifacts per subclusters per source of artifacts, calculate the scCLINIC score for each subclusters and each source of artifacts
  summary_clusterref <- referenceoverlap %>% #sum the Identity Score for each artifacts in the subcluster per source of artifacts aka Subcluster scCLINIC Score per source of artifacts, count the number of unique artifacts in the subcluster per source of artifacts
    group_by(cluster_local, globalcluster,cluster_ref) %>%
    summarise(count = n(),mean(dp1))

  summary_clusterref$cluster_local <- as.matrix(summary_clusterref)[,"cluster_local"]
  #summary_clusterref$cluster_local <- as.integer(summary_clusterref$cluster_local)
  summary_clusterref$cluster_ref <- as.matrix(summary_clusterref)[,"cluster_ref"]
  #summary_clusterref$cluster_ref <- as.integer(summary_clusterref$cluster_ref)
  #summary_clusterref$globalcluster <- as.integer(summary_clusterref$globalcluster)
  summary_clusterref$DiffDP1 <- summary_clusterref$`mean(dp1)`

  #Extend the source of artifacts information to each contaminated subclusters
  GlobalPlusCellSpecific <- inner_join(contamdfall, summary_clusterref,
                                       by = c("cluster_local", "globalcluster"))
  GlobalPlusCellSpecific <- GlobalPlusCellSpecific[,c("cluster_local", "globalcluster","scCLINICScore", "cluster_ref","DiffDP1")]


  ################################
  #Store scCLINIC Result in Metadata of the original seurat object
  #scCLINIC Score (scCLINICScore), Major Cluster ID (L0), Subcluster ID (L1),
  #scCLINIC level (scCLINIC_Level)
  message("Saving results in Seurat Object.")
  obj@meta.data$"scCLINICScore" <- 0
  obj@meta.data$"L0" <- NA
  obj@meta.data$"L1" <- NA
  obj@meta.data$"scCLINIC_Level" <- length(maxima_y)+1

  contaminatedinfo <- GlobalPlusCellSpecific
  contaminatedinfo <- contaminatedinfo[order(contaminatedinfo$"DiffDP1", decreasing = TRUE), ]

  refinfo <- contaminatedinfo %>% #Summarized the source of artifacts for each subclusters
    group_by(cluster_local, globalcluster) %>%
    summarise(cluster_ref = list(cluster_ref))

  contamV <- contamdfall_backup#Including the scCLINIC score for each subclusters (for both contaminated and non-contaminated subclusters)
  for (i in sort(qcsclst)){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".rds")
    if (file_name %in% list.files(folder_path_Step2_L1R)){
      recluster <- readRDS(paste0(folder_path_Step2_L1R,file_name))

      obj$L0 <- ifelse(colnames(obj) %in% colnames(recluster), i ,obj@meta.data$L0)#Store L0

      for (clusterid in unique(recluster$"scCLINIC_subcluster")){
        obj$L1 <- ifelse(colnames(obj) %in% rownames(recluster@meta.data[recluster@meta.data$"scCLINIC_subcluster" == clusterid,]), clusterid ,obj@meta.data$L1)#Store L1
      }

      subclusterlst <- contamV[contamV$globalcluster == i,]$cluster_local
      if (length(subclusterlst)>0){
        for (subclusterid in subclusterlst){
          recluster.filterout <- rownames(recluster@meta.data[recluster@meta.data$scCLINIC_subcluster == subclusterid,]) #cell name of each contaminated subclusters
          curDP1_sum <- contamV[contamV$globalcluster == i & contamV$cluster_local == subclusterid,]$scCLINICScore #curDP1_sum is the scCLINIC Score for the current subcluster (subclusterid)
          obj@meta.data[recluster.filterout,]$scCLINICScore <- rep(curDP1_sum, length(recluster.filterout))#Store scCLINIC Score for each subclusters
        }
      }
    }
  }

  #Store the scCLINIC_Level for each subclusters
  contam_level = length(maxima_y)
  for (maxima_y_cur in sort(maxima_y, decreasing = FALSE)){
    obj$scCLINIC_Level <- ifelse(obj$scCLINICScore >= maxima_y_cur,contam_level, obj$scCLINIC_Level)
    contam_level = contam_level - 1
  }
  obj@meta.data$"scCLINIC_ClusterID" <- paste0(obj@meta.data$"L0","_",obj@meta.data$"L1")#Store scCLINIC_ClusterID

  #Plotting and visualize the scCLINIC results on the UMAP,
  palette15 <- c(
    "#7c3b36",
    "#71b84e",
    "#7c4ccb",
    "#c4944a",
    "#cf4ba6",
    "#597b4a",
    "#d55046",
    "#53a0b9",
    "#4e336c",
    "#b488bb",
    "#8d8d0f",
    "#3b8d89",
    "#e51c23",
    "#1e88e5",
    "#ff9800"
  )
  obj$scCLINIC_Level <- factor(obj$scCLINIC_Level, levels = sort(unique(obj$scCLINIC_Level),decreasing=F))
  col <- setNames(palette15, levels(obj$scCLINIC_Level))

  p4 <- DimPlot(obj, reduction = "umap",group.by = "scCLINIC_Level",raster=FALSE,pt.size = 0.1,cols = col)#rev(viridis(length(unique(obj$scCLINIC_Level)))))
  p6 <- DimPlot(obj, reduction = "umap",group.by = OverlapRatio,raster=FALSE,pt.size = 0.1)
  p8 <- FeaturePlot(obj, reduction = "umap",features = "scCLINICScore",pt.size = 0.1)
  ggsave(filename = paste0(folder_path_Step2_Output,"scCLINICResult.png"), p6+p4+p8, height = 5, width = 15, dpi = 300)

  ##Store score in Subcluster (sync with updated seurat object), OPTIONAL
  for (i in cluster_to_consider){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".rds")
    if (file_name %in% list.files(folder_path_Step2_L1R)){
      recluster <- readRDS(paste0(folder_path_Step2_L1R,file_name))
      recluster@meta.data <- obj@meta.data[rownames(recluster@meta.data),]#copy-paste the metadata of updated seurat object to subcluster's metadata
      #Plotting and visualize the scCLINIC results on the subcluster's UMAP,
      p4 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_Level",raster=FALSE,  sizes.highlight = 0.1,cols=col)
      p7 <- DimPlot(recluster, reduction = "umap",group.by = "scCLINIC_ClusterID",raster=FALSE,  sizes.highlight = 0.1)
      p8 <- FeaturePlot(recluster, reduction = "umap",features = "scCLINICScore")
      ggsave(filename = paste0(folder_path_Step2_Output,i,"_scCLINICResult.png"), p7+p4+p8, height = 5, width = 15, dpi = 300)
    }
  }

  ###Step4 Interpretation
  #Load the contamination information
  file_name <- paste0("ArtifactsInfo.csv")
  curcluster <- read.table(paste0(folder_path_Step2_Output,file_name), sep= ",", header = T, row.names = 1)
  source_of_artifact <- unique(curcluster$cluster_ref)

  #
  artifact_info_lst <- t(data.frame(row.names = c(colnames(curcluster),"multiplehit","globalduplicate","nduplicate")))
  for (source_of_artifact_i in source_of_artifact){#For each source of artifacts (source_of_artifact_i)
    artifact_info <- curcluster[curcluster$cluster_ref== source_of_artifact_i,]#Subset each artifacts from the same source (source_of_artifact_i)
    artifact_info$multiplehit <- duplicated(artifact_info$gene) | duplicated(artifact_info$gene, fromLast = T)# Identify artifacts from the same source which present multiple times (multiplehit = TRUE)

    #Investigate the artifacts from the same source which present multiple times
    dupligenetable <- artifact_info[artifact_info$multiplehit==T,]
    artifact_info$globalduplicate <- 1
    artifact_info$nduplicate <- 1
    for (gene_dupli in unique(dupligenetable$gene)){
      ndupli <- length(dupligenetable[dupligenetable$gene == gene_dupli,]$gene)#ndupli: Number of times the artifacts from the same source are present.
      globaldupli <- length(unique(dupligenetable[dupligenetable$gene == gene_dupli,]$globalcluster))#globaldupli: The major clusters of contaminated subclusters where the artifacts are present.
      #Store both information in artifact_info
      artifact_info$nduplicate <- ifelse(artifact_info$gene == gene_dupli, ndupli ,artifact_info$nduplicate)
      artifact_info$globalduplicate <- ifelse(artifact_info$gene == gene_dupli, globaldupli ,artifact_info$globalduplicate)
    }
    artifact_info_lst <- rbind(artifact_info_lst, artifact_info)#rbind each source of artifacts' artifact_info
  }

  artifact_info <- artifact_info_lst

  #Change the name of the columns
  colnames(artifact_info) <- c(

    "Artifact_Gene" ,   "p_val"  ,    "avg_log2FC" ,"pct.1"  ,    "pct.2"  ,
    "p_val_adj",  "Subcluster" ,   "Enrichment_Score"    ,     "Major_Cluster" ,   "SoA_p_val"  ,
    "SoA_avg_log2FC" ,  "SoA_pct.1"  ,      "SoA_pct.2"   ,     "SoA_p_val_adj" ,   "Source_of_Artifacts_SoA" ,
    "SoA_Enrichment_Score"         ,  "No. cell"      ,     "SoA_No. cell"      ,   "Avg_Exp"   ,  "SoA_Avg_Exp"   ,
    "dp1"          ,    "multiplehit"   ,   "globalduplicate" , "nduplicate"
  )

  #Remove certain columns
  artifact_info <- artifact_info[, !colnames(artifact_info) %in% c("dp1", "multiplehit", "globalduplicate", "nduplicate")]

  # Reorder the columns in artifact_info
  new_order <- c(

    "Artifact_Gene",

    "Major_Cluster", "Subcluster", "No. cell", "Enrichment_Score",  "avg_log2FC", "pct.1", "pct.2", "p_val","p_val_adj", "Avg_Exp",

    "Source_of_Artifacts_SoA","SoA_No. cell", "SoA_Enrichment_Score","SoA_avg_log2FC", "SoA_pct.1", "SoA_pct.2", "SoA_p_val","SoA_p_val_adj","SoA_Avg_Exp"

    )
  artifact_info <- artifact_info[, new_order, drop = FALSE]

  file_name <- paste0("ArtifactsInfo.csv")
  write.table(artifact_info,file = paste0(folder_path_Step2_Output,file_name),sep = ",")

  ###Add the scCLINIC Score for each source of artifacts into the seurat object metadata
  #scCLINIC Subcluster ID
  artifact_info$MajorSub <- paste0(artifact_info$Major_Cluster,"_",artifact_info$Subcluster)
  #Summarize ES score and their source of major cluster for each subclusters
  result <- artifact_info %>%
    group_by(MajorSub) %>%
    summarize(
      cluser_reflst = list(Source_of_Artifacts_SoA),
      dp1lst = list(Enrichment_Score)
    ) %>%
    ungroup()

  # Function to calculate the average ES score (dp1lst) for each source of artifacts (cluster_reflst)
  tabulate_Cluster_Contribution_Score <- function(cluster_reflst, dp1lst) {
    components <- unlist(cluster_reflst)
    dp1_values <- unlist(dp1lst)
    CCS_Matrix <- tapply(dp1_values, components, mean, na.rm = TRUE)
    return(CCS_Matrix)
  }

  # For each subclusters (each row in result), calculate the average ES score for each source of artifacts
  tabulate_CCS <- mapply(tabulate_Cluster_Contribution_Score, result$cluser_reflst, result$dp1lst, SIMPLIFY = FALSE)

  # List of all source of artifacts which contaminated major cluster X
  lst_of_source_of_artifacts <- unique(unlist(lapply(tabulate_CCS, names)))

  # Create a matrix, each row represent one subclusters and each column present each source of artifacts, to store the CCS for each major clusters
  CCS_Matrix <- matrix(NA, nrow = length(tabulate_CCS), ncol = length(lst_of_source_of_artifacts), dimnames = list(NULL, lst_of_source_of_artifacts))

  # store the CCS for each major clusters in the matrix
  for (i in seq_along(tabulate_CCS)) {
    CCS_Matrix[i, names(tabulate_CCS[[i]])] <- tabulate_CCS[[i]]
  }

  # Replace NA with 0
  CCS_Matrix[is.na(CCS_Matrix)] <- 0

  CCS_Matrix <- as.data.frame(CCS_Matrix)
  CCS_Matrix$MajorSub <- result$MajorSub #Named each rows with their Subclusters ID
  #Rowname = MajorSub and remove the column
  rownames(CCS_Matrix) <- CCS_Matrix$MajorSub
  CCS_Matrix <- CCS_Matrix[, !colnames(CCS_Matrix) %in% "MajorSub"]
  #Add each source of artifacts into metadata as seperate column
  for (colname in colnames(CCS_Matrix)){
    obj@meta.data[[colname]]<-0
  }
  #Assign each source of artifacts' scCLINIC score into metadata for each subclusters
  for (rowccs in rownames(CCS_Matrix)){
    for (colname in colnames(CCS_Matrix)){
      obj@meta.data[[colname]] <- ifelse(obj$scCLINIC_ClusterID == rowccs, CCS_Matrix[rowccs,colname] ,obj@meta.data[[colname]])#Store L0
    }
  }
  #Tabulated subclusters information, including their scCLINIC_ClusterID, scCLINIC_Level, PercentageDroplets, ElbowPoint, n_artifacts, n_cells, scCLINIC score, CCS

  scoreinfo <- obj@meta.data[,c(c("scCLINIC_ClusterID","scCLINIC_Level"),colnames(CCS_Matrix))]
  scoreinfo <- scoreinfo[!duplicated(scoreinfo), ]
  scoreinfo <- scoreinfo[complete.cases(scoreinfo), ]
  rownames(scoreinfo) <- NULL

  contamdfall$scCLINIC_ClusterID <- paste0(contamdfall$globalcluster,"_",contamdfall$cluster_local)
  contamdfall$PercentageDroplets <- contamdfall$`mean(ncells)`/totalcells*100
  contamdfall <- merge(contamdfall, scoreinfo, by = "scCLINIC_ClusterID", all = TRUE)
  contamdfall$ElbowPoint <- 0
  contam_level <- length(maxima_y)
  for (maxima_y_cur in sort(maxima_y,decreasing = F)){
    contamdfall$ElbowPoint <- ifelse(contamdfall$scCLINIC_Level == contam_level,maxima_y_cur,contamdfall$ElbowPoint)
    contam_level <- contam_level -1
  }
  contamdfall <- transform(contamdfall, n_artifacts = count, n_cells = `mean(ncells)`)
  contamdfall <- subset(contamdfall, select = -c(count  , mean.dp1., mean.ncells.))

  #Remove certain columns
  contamdfall <- contamdfall[, !colnames(contamdfall) %in% c("cluster_local", "globalcluster","n_artifacts")]

  #Change the name of the columns
  colnames(contamdfall)[colnames(contamdfall) == "n_cells"] <- "No. cell"
  colnames(contamdfall)[colnames(contamdfall) == "PercentageDroplets"] <- "No. cell (%)"
  colnames(contamdfall)[colnames(contamdfall) == "ElbowPoint"] <- "scCLINIC_Elbow_Point"
  # Identify columns with "^M" in their names
  cols_with_M <- grep("\\M", colnames(contamdfall), value = TRUE)
  # Assign new column names to the data frame
  colnames(contamdfall)[colnames(contamdfall) %in% cols_with_M] <- paste0("CCS","_", cols_with_M)

  # Reorder the columns in artifact_info
  new_order <- c(
    "scCLINIC_ClusterID" ,  "scCLINICScore" , "No. cell (%)" ,"No. cell",  "scCLINIC_Level", "scCLINIC_Elbow_Point",grep("\\CCS", colnames(contamdfall), value = TRUE)
  )

  contamdfall <- contamdfall[, new_order, drop = FALSE]

  write.table(contamdfall,file = paste0(folder_path_Step2_Output,"scCLINICScore.csv"),sep = ",",row.names = FALSE)

  #Tidy up Metadata (Remove unnecessary columns)
  obj_metadata <- obj@meta.data
  metadata_saved <- c(OverlapRatio,"scCLINICScore","scCLINIC_Level","scCLINIC_ClusterID",grep("^M",colnames(obj_metadata),value = T))
  obj_metadata <- obj_metadata[,metadata_saved]
  # Get the column names
  colnames <- colnames(obj_metadata)
  # Add "CCS_" in front of column names that start with "M"
  colnames(obj_metadata) <- ifelse(grepl("^M", colnames), paste0("CCS_", colnames), colnames)
  obj@meta.data <- obj_metadata

  #Save the updated seurat object (including scCLINIC results)
  saveRDS(obj,file = paste0(folder_path_Step2_Output,"scCLINICResult.rds"))

  message("Step2B completed.")

  return(obj)
}
