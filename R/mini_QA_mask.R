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
#' @import dplyr
#' @import sf
#' @export

mini_QA_mask<-function(DTM_fp,outdir,mini_segment_sf_list,slope_thres = 45, sos_thres = 84.5, save_mask = F){
  
  
  # 1. Define file paths for the intermediate slope products
  slope_path <- file.path(outdir, "Slope_for_QA.tif")
  slope_of_slope_path <- file.path(outdir, "Slope_of_slope_for_QA.tif")
  mask_path <- file.path(outdir, "QA_Mask.tif")
  
  # 2. Create slope and slope of slope rasters
  gdaldem(mode = "slope",
          input_dem = DTM_fp,
          output_map = slope_path)
  
  gdaldem(mode = "slope",
          input_dem = slope_path,
          output_map = slope_of_slope_path)
  
  #3. Load the rasters
  
  s <- rast(slope_path)
  ss <- rast(slope_of_slope_path)
  
  #4. Add QA column
  
  treetops<-mini_segment_sf_list[[1]]
  
  slope_ext<-terra::extract(s,treetops)
  sos_ext<-terra::extract(ss,treetops)
  
  treetops<-treetops %>% 
    dplyr::mutate(slope = slope_ext$Slope_for_QA,
                  slope_of_slope = sos_ext$Slope_of_slope_for_QA,
                  .after = Z) %>% 
    dplyr::mutate(QA_pass = ifelse(slope<=slope_thres & slope_of_slope<=sos_thres,"pass","fail")) %>% 
    dplyr::select(-c(slope,slope_of_slope))
  
  crowns<-mini_segment_sf_list[[2]]
  QA<-sf::st_drop_geometry(treetops) %>% 
    dplyr::select(treeID,QA_pass)
  
  crowns<-left_join(crowns,QA,by="treeID")
  
  mini_segment_sf_list<-list(seg_treetops = treetops,
                             seg_crowns = crowns)
  
  return(mini_segment_sf_list)
  
  if(save_mask==T){
    
    # 4. Create the mask
    # Logical: Keep pixels where (Slope <= 45) AND (Slope_of_Slope <= 84.5)
    # Areas failing this (extreme slopes/noise) become NA
    qa_mask <- (s <= slope_thres) & (ss <= sos_thres)
    
    # Convert FALSE (0) to NA so only the "Good" areas remain as 1
    qa_mask[qa_mask == 0] <- NA
    
    # 5. Save the final mask
    writeRaster(qa_mask, mask_path, overwrite = TRUE)
    
  }
  unlink(slope_path)
  unlink(slope_of_slope_path)
  
}