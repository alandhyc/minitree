#' Segmenting small trees
#'
#' A function to create segment trees using LiDAR point clouds as an input. Specifically designed to be sensitive to small regenerating trees.
#'
#'
#' @param shm_r A spatRaster of canopy heights created by the lidarSHM pipeline, created by increasing max_shrub_ht to maximum height of trees
#' @return A List with two objects. (1) sf object showing the treetops (POINTS) and (2) sf object showing the segmented crowns (POLYGONS)
#' @import terra
#' @import lidR
#' @import sf
#' @import dplyr
#' @export

mini_segment<-function(shm_r){
  
  #Step 1. Smoothing the CHM created by the SHM algorithm
  
  kernel<-matrix(1,3,3)
  smoothed<-terra::focal(shm_r, w=kernel, fun = mean, na.rm = TRUE)
  smoothed<-terra::ifel(smoothed<0,0,smoothed)
  
  #Step 2: Create shrub-adjusted CHM that removes shrub height
  local_shrub_height <- terra::focal(shm_r, w = matrix(1, 19, 19), fun = mean, na.rm = TRUE)
  
  valid_low_veg <- shm_r <= 5 & !is.na(shm_r) & !is.na(local_shrub_height)
  adjusted <- shm_r - local_shrub_height
  adjusted <- terra::app(adjusted, function(x) pmax(x, 0))
  relative_chm <- shm_r  
  relative_chm <- terra::lapp(
    c(shm_r, adjusted, valid_low_veg),
    fun = function(c, a, m) ifelse(m, a, c)
  )
  
  #Step 2: Use smoothed CHM to detect treetops
  
  #Function for locate_trees
  
  f <- function(x) {
    y <- 3*(-(exp(-0.15*(x-1))-1))+2
    return(y)
  }
  
  #Locate trees
  
  if(terra::ncell(smoothed)==0) return(NULL)
  
  treetops <- lidR::locate_trees(smoothed, lmf(ws = f, hmin = 0))
  
  #Step 3: Extract heights from shrub adjusted CHM for filtering
  
  if(nrow(treetops)==0) return(NULL)
  
  treetops$treeID<-1:nrow(treetops)
  
  relative_z<-extract(relative_chm,treetops)
  names(relative_z)<-c("treeID","Z")
  
  treetops<-treetops %>% 
    dplyr::select(-Z) %>% 
    left_join(relative_z,by = "treeID")
  
  treetops<-treetops[which(treetops$Z>0.5),]
  
  if(nrow(treetops)==0) return(NULL)
  
  #Step 4: Grow the crowns using the smoothed CHM
  
  smoothed<-terra::toMemory(smoothed)
  
  crowns <- lidR::dalponte2016(
    smoothed,
    treetops,
    th_tree = 0.5, #0.5
    th_seed = 0.4, #0.4
    th_cr = 0.55, #0.55
    max_cr = 100, #100
    ID = "treeID"
  )()
  
  crown_seg<-terra::as.polygons(crowns)
  crown_seg<-sf::st_as_sf(crown_seg)
  
  names(crown_seg)[1]<-"treeID"
  
  crown_seg<-crown_seg %>% 
    dplyr::left_join(sf::st_drop_geometry(treetops),by = "treeID")
  
  #Step 5: Get original heights for the trees (not the ones without shrubs)
  
  #Currently the Z values are from the ralative CHM (shrub removed)
  #Convert it back to values from the CHM
  
  accurate_z<-terra::extract(shm_r,treetops)
  
  accurate_z$ID<-treetops$treeID[accurate_z$ID]
  names(accurate_z)<-c("treeID","adjusted_z")
  
  crown_seg<-crown_seg %>% 
    dplyr::left_join(accurate_z,by = "treeID")
  
  crown_seg<-crown_seg %>% 
    dplyr::mutate(area = as.numeric(st_area(crown_seg)),.after = "adjusted_z") %>% 
    dplyr::select(-Z) %>% 
    dplyr::rename(Z = adjusted_z)
  
  seg<-list(
    seg_treetops = treetops,
    seg_crowns = crown_seg
  )
  
  return(seg)
  
}