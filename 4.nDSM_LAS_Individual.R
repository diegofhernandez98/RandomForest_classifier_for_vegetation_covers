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


# RUN ONLY WHEN CREATING NEW LAS

las_fold <- sub.folders[[14]]
fold_name <- basename(las_fold)
las_files <- list.files(las_fold, pattern = "\\.las$", full.names = TRUE)
ctg <- readLAScatalog(las_files)
las <- readLAS(ctg)



dsm <- rasterize_canopy(las, res = 0.1, algorithm = dsmtin())
plot(dsm)

writeRaster(dsm, here("LAS_Processing/Transects", paste0("DSM_transect_PTrial.tif")),
            overwrite = TRUE)


dtm <- rasterize_terrain(las, res = 0.1, algorithm = tin())
plot(dtm)

ndsm <- dsm - dtm
plot(ndsm)
writeRaster(ndsm, here("LAS_Processing/Transects", paste0("ndsm_transect_PTrial.tif")),
            overwrite = TRUE)
