library(sf)
library(tidyverse)
library(here)
library(raster)
library(ranger)
library(caret)
library(terra)

#Set the training directory
tdir <- here("training\\")

target_crs <- "EPSG:7850"

#List all shapefiles in the training directory
f.list <- list.files(tdir, pattern = "shp$", full.names = FALSE)

#Initialize an empty data frame for training data
trn.all <- st_read(paste0(tdir, f.list[1])) %>% 
  mutate(class = str_split_fixed(f.list[1], "\\.", 2)[,1])
trn.all <- trn.all[0,]

#Read and combine all training shapefiles
for (i in 1:length(f.list)){
  trn <- st_read(paste0(tdir, f.list[i])) %>% 
    mutate(class = str_split_fixed(f.list[i], "\\.", 2)[,1])
  if (nrow(trn) > 50){
    trn <- trn[1:50, ]
  }
  trn.all <- rbind(trn.all, trn)
}

#Plot the combined training data

plot(trn.all[,2])


#Load the DSM and ortho images
nsm <- rast(here("Transects/ndsm_transect_N01.tif"))
plot(nsm)
nsm_reprojected <- project(nsm, target_crs, method = "near")
img <- brick(here("Transects/ortho_transect_N01.tif"))[[c(1:3)]] 
plot(img)

#Ensure CRS match
if (st_crs(trn.all) != crs(nsm_reprojected)) {
  trn.all <- st_transform(trn.all, crs(nsm_reprojected))
}

extracted_values<- raster::extract(nsm_reprojected, st_coordinates(trn.all))

# Ensure the extracted values are a simple vector
extracted_values <- unlist(extracted_values)

trn.all$nsm <- extracted_values

head(trn.all)
st_write(trn.all, here("points_shapefile.shp"), delete_layer = TRUE)




