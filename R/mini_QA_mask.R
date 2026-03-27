#' QA mask for minitree
#'
#' A function that takes the digital terrain model (DTM) as input and returns recommended QA masks
#'
#'
#' @param DTM_fp Filepath to the DTM. A DTM of 1m ground resolution was used in Thapa et al. (2026)
#' @param out_dir Output directory.
#' @return A List with two objects. (1) sf object showing the treetops (POINTS) and (2) sf object showing the segmented crowns (POLYGONS)
#' @import gdalUtilities
#' @import terra
#' @export

mini_QA_mask<-function(DTM_fp,outdir){
  
  gdaldem(mode = "slope",
          input_dem = DTM_fp,
          output_map = paste0(out_dir,"/Slope_for_QA.tif"))
  
  gdaldem(mode = "slope",
          input_dem = paste0(out_dir,"/Slope_for_QA.tif"),
          output_map = paste0(out_dir,"/Slope_of_slope_for_QA.tif"))
  
  
  
}