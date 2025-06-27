#' Plot Contamination Pattern
#' @param Name Name of Run
#' @param obj Seurat Object return by Step1D
#' @param Outdir Output Directory of scCLINIC Results
#' @param OverlapRatio Step1B: The overlap ratio to carry forward to Step1C, default using 0.5, decrease to have coarse resolution. User manual annotation metadata colname, if CELLANNOTATION = TRUE,
#' @param CELLANNOTATION default is FALSE, TRUE if using user cluster annotation
PlotContaminationPattern <- function(obj,Outdir,Name,OverlapRatio=0.5,CELLANNOTATION = FALSE, verbose = TRUE){
  ###No Return
  message("Plot Contamination Pattern.")

  if (CELLANNOTATION){
    message("Using user-annotated clusters.")
    OverlapRatio <- "annotation_index" #User manual cellannotation
  }else{
    OverlapRatio <- paste0("Overlap_Ratio_",OverlapRatio)
    message(paste0("Using ",OverlapRatio))
  }

  folder_path_Step2_Output <- paste0(Outdir,Name,"_Step2/Output_",OverlapRatio,"/")

  #Load artifacts information and scCLINIC result rds
  contamgeneinfo <- read.csv(paste0(folder_path_Step2_Output,"ArtifactsInfo.csv"),na.strings = c("", "NA"))

  #scCLINIC Subcluster ID
  contamgeneinfo$MajorSub <- paste0(contamgeneinfo$Major_Cluster,"_",contamgeneinfo$Subcluster)
  #Summarize ES score and their source of major cluster for each subclusters
  result <- contamgeneinfo %>%
    group_by(MajorSub) %>%
    summarize(
      cluser_reflst = list(Source_of_Artifacts_SoA),
      dp1lst = list(Enrichment_Score)
    ) %>%
    ungroup()

  # Sort the data frame based on numeric row names
  contamgeneinfo <- contamgeneinfo[order(as.numeric(rownames(contamgeneinfo))), ]
  #Remove Artifacts with ES = 0
  contamgeneinfo <- contamgeneinfo[contamgeneinfo$Enrichment_Score > 0, ]
  contamgeneinfo <- contamgeneinfo[, !colnames(contamgeneinfo) %in% c("MajorSub")]
  #Update ArtifactsInfo.csv
  write.table(contamgeneinfo,file = paste0(folder_path_Step2_Output,"ArtifactsInfo.csv"),sep = ",", row.names = FALSE)


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

  # Convert to dataframe for plotting
  CCS_Matrix <- as.data.frame(CCS_Matrix)
  CCS_Matrix$MajorSub <- result$MajorSub #Named each rows with their Subclusters ID
  data <- CCS_Matrix %>%
    pivot_longer(cols = -MajorSub, names_to = "Component", values_to = "value") #dependency tidyr

  for (clusterx in unique(obj@meta.data[,OverlapRatio])){
    # Filter out other clusters, only keep cluster X which wish to display
    filtered_data <- data %>%
      dplyr::filter(grepl(paste0("^", clusterx, "(?![0-9])"), MajorSub, perl = TRUE))

    if (nrow(filtered_data) != 0){
      # Heatmap
      # Convert dataframe
      heatmap_data <- dcast(filtered_data, MajorSub ~ Component, value.var = "value")
      heatmap_data <- heatmap_data[, !colnames(heatmap_data) %in% clusterx]#remove itself as source of artifacts

      rownames(heatmap_data) <- heatmap_data$MajorSub #Named rownames as the subcluster ID (MajorSub)
      heatmap_data$MajorSub <- NULL #remove MajorSub columns

      # Summarize scCLINIC Score (CS) of each subclusters
      CS_table <- obj@meta.data %>%
        group_by(scCLINIC_ClusterID) %>%
        summarize(
          scCLINICScore = mean(scCLINICScore), #scCLINICScore of each cells within each subclusters are same value...
        ) %>%
        ungroup()

      #Convert to plotting dataframe
      heatmap_data <- as.data.frame(heatmap_data)  # Convert heatmap_data to a dataframe
      heatmap_data <- heatmap_data %>%
        mutate(Score = CS_table$scCLINICScore[match(rownames(heatmap_data), CS_table$scCLINIC_ClusterID)]) #Score

      # Convert data to matrix
      heatmap_matrix <- as.matrix(heatmap_data)

      # Reorder and prepare the data
      multicolplot <- data.frame((heatmap_matrix))

      rownames_sorted <- rownames(multicolplot)[order(as.numeric(gsub(".*S", "", rownames(multicolplot))), decreasing = TRUE)]
      multicolplot <- multicolplot[rownames_sorted, ]

      colnames_sorted <- colnames(multicolplot)[order(as.numeric(gsub("M", "", colnames(multicolplot))), decreasing = FALSE)]
      multicolplot <- multicolplot[, colnames_sorted]

      multicolplot$Category <- as.character(rownames(multicolplot))
      multicolplot$Category <- factor(multicolplot$Category, levels = rownames_sorted)

      x_limits <- c(0, ceiling(max(multicolplot[, -ncol(multicolplot)], na.rm = TRUE) * 100) / 100)
      custom_breaks <- seq(x_limits[1], x_limits[2], length.out = 2)
      text_size <- 8

      # Create individual plots without y-axis text
      plot_list <- list()
      for (i in seq_along(colnames_sorted)) {
        clus <- colnames_sorted[i]
        p <- ggplot(multicolplot, aes_string(x = "Category", y = clus)) +
          geom_bar(stat = "identity", fill = viridis(length(colnames_sorted))[i]) +
          labs(title = paste(clus), y = "", x = "") +
          theme_minimal() +
          theme(
            panel.grid = element_line(color = "grey", size = 0.5),
            axis.line = element_line(color = "black"),
            axis.text = element_text(size = text_size),
            axis.title = element_text(size = text_size),
            plot.title = element_text(size = text_size),
            #axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            plot.margin = unit(c(0,0,0,0), "cm")
          ) +
          coord_flip() +
          scale_y_continuous(breaks = custom_breaks, limits = x_limits)

        plot_list[[i]] <- p
      }

      # Create a plot for y-axis text labels
      y_axis_plot <- ggplot(multicolplot, aes_string(x = "Category", y = 1)) +
        geom_blank() +
        #geom_bar(stat = "identity") +
        labs(title = "", y = "", x = "") +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.text = element_text(size = text_size),
          axis.title = element_text(size = text_size),
          plot.title = element_text(size = text_size),
          axis.text.x = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm")
        ) +
        coord_flip() +
        scale_y_continuous(breaks = custom_breaks, limits = x_limits)

      # Combine the y-axis labels plot with the other plots
      #combined_plots <- cowplot::plot_grid(
      #  y_axis_plot,
      #  cowplot::plot_grid(plotlist = plot_list, ncol = length(plot_list)),
      #  rel_widths = c(1, 15)
      #)
      combined_plots <- y_axis_plot + 
      patchwork::wrap_plots(plot_list,nrow=1) + 
      patchwork::plot_annotation(title = "Artifact Source",theme = theme(plot.title = element_text(size=text_size,hjust = 0.5))) + 
      patchwork::plot_layout(nrow=1,widths = c(1,10000))
      # Save the combined plot
      ggsave(filename = paste0(folder_path_Step2_Output, clusterx, "_scCLINICScore.png"), 
             combined_plots, height = 0.2*length(rownames_sorted), width = 1.2*length(colnames_sorted), dpi = 300, limitsize = FALSE)
    }
  }

  ###Plot Source of Artifacts Contamination Patterns
  #Plotting and visualize the scCLINIC score of each source of artifacts on the major cluster's UMAP,
  p1 <- FeaturePlot(obj,features = paste0("CCS_",colnames(CCS_Matrix)[colnames(CCS_Matrix) != "MajorSub"]))
  ggsave(filename = paste0(folder_path_Step2_Output,"ContaminationPattern_SourceOfArtifacts",".png"),p1,  height = 20, width = 20, dpi = 300)

  folder_path_Step2_L1R <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/")
  folder_path_Step2_L1R_Marker <- paste0(Outdir,Name,"_Step2/",OverlapRatio,"_recluster/Marker/")

  qcsclst <- sort(unique(obj@meta.data[,OverlapRatio]))#list of major clusters seurat ID
  cluster_to_consider <- list()
  for (i in sort(qcsclst)){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".csv")
    if (file_name %in% list.files(folder_path_Step2_L1R_Marker)){
      cluster_to_consider <- unlist(c(cluster_to_consider,i))
    }
  }

  for (i in cluster_to_consider){
    file_name <- paste0(OverlapRatio,"_cluster_",i,".rds")
    if (file_name %in% list.files(folder_path_Step2_L1R)){
      recluster <- readRDS(paste0(folder_path_Step2_L1R,file_name))
      recluster@meta.data <- obj@meta.data[rownames(recluster@meta.data),]#copy-paste the metadata of updated seurat object to subcluster's metadata
      #Plotting and visualize the scCLINIC score of each source of artifacts on the subcluster's UMAP,
      p1 <- FeaturePlot(recluster, features = paste0("CCS_",colnames(CCS_Matrix)[!colnames(CCS_Matrix) %in% c("MajorSub", i)]))
      ggsave(filename = paste0(folder_path_Step2_Output,i,"_SourceofArtifacts.png"), p1, height = 10, width = 10, dpi = 300)
    }
  }
}
