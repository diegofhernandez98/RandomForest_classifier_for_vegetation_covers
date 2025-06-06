library(lidR)
library(ggplot2)
library(raster)
library(terra)
library(sf)
library(here)

#reading las
target_crs <- "EPSG:7850"
tdir <- here::here("Z:\\DEC\\ForestManagementPlan_2014_2023\\DATA\\Working\\SiteResistivity_Balingup\\Outputs\\LAS_FILES")

#https://stackoverflow.com/questions/20799989/loop-over-sub-folder-in-r
sub.folders <- list.dirs(tdir, recursive=TRUE)[-1]
folder_names <- basename(sub.folders)

#Reading buffer of transect
tdir2 <- here("shp\\")
f.list2 <-list.files(tdir2, pattern = "buffer.*\\.shp$", full.names = TRUE)

buffer_list <- list()
for (i in seq_along(f.list2)) {
  buffer_zone <- st_read(f.list2[i]) %>%
    st_transform(crs = target_crs)
  
  name_site <- buffer_zone$Site

  #filtering the buffers by folder name  
  if (name_site %in% folder_names) {
    buffer_list[[name_site]] <- buffer_zone
    
  }
}

# RUN ONLY WHEN CREATING NEW LAS
las_list <- list()

for (i in seq_along(sub.folders)) {
  las_fold <- sub.folders[[i]]
  fold_name <- basename(las_fold)
  las_files <- list.files(las_fold, pattern = "\\.las$", full.names = TRUE)
  ctg <- readLAScatalog(las_files)
  # las_check(ctg)
  study_area <- buffer_list[[i]]
  buffer_geom <- st_geometry(study_area)

  las_clipped <- clip_roi(ctg, buffer_geom)
  las_list[[fold_name]] <- las_clipped
  saveRDS(las_clipped,
          here("LAS_Processing", "las", paste0("LAS_",fold_name,".rds")))
}
####################################################################################


las_dir <- here::here("Z:\\DEC\\ForestManagementPlan_2014_2023\\DATA\\Working\\SiteResistivity_Balingup\\LAS_Processing\\las")
rds_files <- list.files(las_dir, pattern = "\\.rds$", full.names = TRUE)



las_list2 <- list()
for (file in rds_files) {
  name <- tools::file_path_sans_ext(basename(file))
  las_list2[[name]] <- readRDS(file)
}

writeLAS(las_list2[[2]], here("LAS_Processing","las","lasfile_n02.las"))

####################################################################################
#dsm creation
ndsm_list <- list()

for (i in seq_along(las_list2)) {
  las <- las_list2[[i]]
  las_fold <- sub.folders[[i]]
  fold_name <- basename(las_fold)
  dsm <- rasterize_canopy(las, res = 0.1, algorithm = dsmtin())
  dtm <- rasterize_terrain(las, res = 0.1, algorithm = tin())
  ndsm <- dsm - dtm
  
  ndsm_list[[i]] <- ndsm
  writeRaster(ndsm, here("LAS_Processing","Transects",paste0("ndsm_transect_",fold_name,".tif")),overwrite = TRUE)
}


folder_names
####################################################################################
#dsm AND dtm creation

for (i in seq_along(las_list2)) {
  las <- las_list2[[i]]
  las_fold <- sub.folders[[i]]
  fold_name <- basename(las_fold)
  dsm <- rasterize_canopy(las, res = 0.1, algorithm = dsmtin())
  dtm <- rasterize_terrain(las, res = 0.1, algorithm = tin())
  plot(dtm)
  # writeRaster(dtm, here("LAS_Processing","Transects",paste0("DTM_transect_",fold_name,".tif")),overwrite = TRUE)
  # writeRaster(dsm, here("LAS_Processing","Transects",paste0("DSM_transect_",fold_name,".tif")),overwrite = TRUE)
  
}




####################################################################################
#dsm models
chm <- rasterize_canopy(las_list[[1]], 0.1, p2r(subcircle = 0.2), pkg = "terra")
# #writeRaster(chm, here("LAS_Processing", paste0("CHM", name_site, ".tif")),
#             overwrite = TRUE)


dsm1 <- rasterize_canopy(las_list[[1]], res = 0.1, algorithm = dsmtin())
# #writeRaster(dsm1 , here("LAS_Processing", paste0("dsm1", name_site, ".tif")),
#             overwrite = TRUE)

dsm2 <- rasterize_canopy(las_list[[1]], res = 0.1, algorithm = p2r())

# #writeRaster(dsm2 , here("LAS_Processing", paste0("dsm2", name_site, ".tif")),
#             overwrite = TRUE)
####################################################################################
dtm1 <- rasterize_terrain(las_list2[[2]], res = 0.1, algorithm = tin())
plot(dtm1)


dtm2 = rasterize_terrain(las_list2[[2]], algorithm = knnidw(k = 6L, p = 2))
plot(dtm2)


