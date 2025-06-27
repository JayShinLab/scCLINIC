#' This function allows you to run scCLINIC in one step.
#' @param Name Name of Run
#' @param Input Input Directory of Seurat Object, or seurat object
#' @param Outdir Output Directory of scCLINIC Results
#' @param resol Step1A: FindClusters resolution for initial clustering in Step1, default using Seurat FindClusters default resolution: 0.8
#' @param overlapRatioList Step1B: The overlap ratio to try in Step1B (recommended range: 0.1 to 0.5)
#' @param OverlapRatio Step1B: The overlap ratio to carry forward to Step1C, default using 0.5, decrease to have coarse resolution. User manual annotation metadata colname, if CELLANNOTATION = TRUE,
#' @param ISThreshold Step1C: The threshold for Identity Score (IS) to classify low quality cells
#' @param gene_n Number of global markers to used
#' @param CELLANNOTATION default is FALSE, TRUE if using user cluster annotation
#' @param rm_tmp default is FALSE, TRUE to remove temporary files
#' @keywords scCLINIC automated all in one function
#' @export
#' @examples
#' scCLINIC_Main_Function("kidney",obj,"~/user/output/")
#' scCLINIC_Main_Function("kidney","~/user/kidney.rds","~/user/output/")
#' scCLINIC_Main_Function("kidney","~/user/kidney.rds","~/user/output/", OverlapRatio = "CT.Park",CELLANNOTATION = TRUE)
scCLINIC_Main_Function <-function(Name,Input,Output,resol=0.8,overlapRatioList=c(0.1,0.25,0.5,0.75,0.9),OverlapRatio=0.5,ISThreshold=0,gene_n=150,CELLANNOTATION = FALSE,rm_tmp = FALSE, SCT = FALSE){

 if (CELLANNOTATION){
    if (is.character(Input) && grepl("\\.RDS$", Input, ignore.case = TRUE)) {
      # Load the .RDS file
      message("Read data started.")
      obj <- readRDS(Input)
    } else {
      obj <- Input
    }
    obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol="Manual",OverlapRatio,gene_n,CELLANNOTATION = TRUE, SCT = SCT)
    obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol="Manual",OverlapRatio,ISThreshold,CELLANNOTATION = TRUE, SCT = SCT)
    STEP2A_Subcluster(obj,Output,Name,resol="Manual",OverlapRatio,gene_n,CELLANNOTATION = TRUE, SCT = SCT)
    obj <- STEP2B_ContaminationScore(obj,Output,Name,resol="Manual",OverlapRatio,gene_n,CELLANNOTATION = TRUE, SCT = SCT)
    PlotContaminationPattern(obj,Output,Name,OverlapRatio,CELLANNOTATION = TRUE)
  }else{
    obj <- STEP1A_GlobalMarkers(Input,Output,Name,resol, SCT = SCT)
    obj <- STEP1B_MergingCluster(obj,Output,Name,resol,overlapRatioList,gene_n, SCT = SCT)
    obj <- STEP1C_RecalculateGlobalMarkers_IdentityScore(obj,Output,Name,resol,OverlapRatio,gene_n, SCT = SCT)
    obj <- STEP1D_FilterLowISCluster(obj,Output,Name,resol,OverlapRatio,ISThreshold, SCT = SCT)
    STEP2A_Subcluster(obj,Output,Name,resol,OverlapRatio,gene_n, SCT = SCT)
    obj <- STEP2B_ContaminationScore(obj,Output,Name,resol,OverlapRatio,gene_n, SCT = SCT)
    PlotContaminationPattern(obj,Output,Name,OverlapRatio)
  }
  if (CELLANNOTATION){
    OverlapRatio <- "annotation_index" #User manual cellannotation
  }else{
    OverlapRatio <- paste0("Overlap_Ratio_",OverlapRatio)
  }
  if(rm_tmp){
    dir_delete(paste0(Output,Name,"_Step1"))
    dir_delete(paste0(Output,Name,"_Step2/",OverlapRatio,"_recluster"))
  }
}
