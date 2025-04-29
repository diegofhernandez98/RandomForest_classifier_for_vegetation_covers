library(sf)
library(tidyverse)
library(here)
library(raster)
library(terra)
library(dplyr)

tdir <- here("outputs/Processing/Balingup_L1_P1/Processing/Agisoft/P1_Collie2/DSM_PRODUCTS")

shp_dir <- here("outputs/Processing/Balingup_L1_P1/shp")


raster_files <- list.files(tdir, pattern = "\\.tif$", full.names = TRUE)
shp_files <- list.files(shp_dir, pattern = "\\.shp$", full.names = TRUE)


dsm_raster <- rast(file.path(tdir, "dsm.tif"))

dtm_raster <- rast(file.path(tdir, "dtm.tif"))

ortho_raster<- rast(file.path(tdir,"Orthomosaic_dsm.tif"))

#shp_obj <- st_read(shp_files)

shp_obj <- st_read(file.path(shp_dir, "Transects2.shp"))


filtered_data <- shp_obj[shp_obj$Site == "Collie02", ] %>%
  sf::st_buffer(30)

filtered_data <- filtered_data %>%
  dplyr::select(Site)

st_write(filtered_data, here("shp/bufferCollie1.shp"), delete_layer = TRUE)
plot(filtered_data)

name_site<- unlist(filtered_data$Site)
##############################################################################################
#ndsm making run first 
# NDSM <- dsm_raster - dtm_raster
# plot(NDSM)
##############################################################################################
#ndsm making forced 
dsm_raster_cropped <- crop(dsm_raster, filtered_data)
dtm_raster_cropped <- crop(dtm_raster, filtered_data)

aligned_dtm <- resample(dtm_raster_cropped, dsm_raster_cropped)
NDSM <- dsm_raster_cropped  - aligned_dtm



NDSM <- dsm_raster_cropped  - dtm_raster_cropped
plot(NDSM)
##############################################################################################
#CROPPING ndsm
cropped<- crop(NDSM,filtered_data)
masked<- terra::mask(cropped,filtered_data)
plot(masked)
writeRaster(masked, here::here(paste0("Transects/ndsm_transect_",name_site,".tif")), overwrite=TRUE)

#CROPPING orthomosaic
cropped2<- crop(ortho_raster,filtered_data)
masked2<- terra::mask(cropped2,filtered_data)
plot(masked2)

writeRaster(masked2, here::here(paste0("Transects/ortho_transect_",name_site,".tif")), overwrite=TRUE)
# writeRaster(masked, 
#             here("transects/ndsm_transect.tif"), overwrite=TRUE)

##############################################################################################

tdir <- "Z:/DEC/ForestManagementPlan_2014_2023/DATA/Working/TuartDecline_Mar2025/RPA_Imagery/Processing/Agisoft/DSM_PRODUCTS"

shp_dir <- "Z:/DEC/ForestManagementPlan_2014_2023/DATA/Working/TuartDecline_Mar2025/shp"

raster_files <- list.files(tdir, pattern = "\\.tif$", full.names = TRUE)
shp_files <- list.files(shp_dir, pattern = "\\clip.shp$", full.names = TRUE)


shp_obj <- st_read(shp_files)

filtered_data <-shp_obj


writeRaster(NDSM, "ndsm.tif", overwrite=TRUE)
writeRaster(cropped2, "Ortho_transect.tif", overwrite=TRUE)

