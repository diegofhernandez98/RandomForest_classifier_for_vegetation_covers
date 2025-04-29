library(sf)
library(tidyverse)
library(here)
library(raster)
library(ranger)
library(caret)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)



# Set the training directory
tdir <- here("training\\")


target_crs <- "EPSG:7850"

# List all shapefiles in the training directory
f.list <- list.files(tdir, pattern = "shp$", full.names = FALSE)

# Initialize an empty data frame for training data
trn.all <- st_read(paste0(tdir, f.list[1])) %>% 
  mutate(class = str_split_fixed(f.list[1], "\\.", 2)[,1])
trn.all <- trn.all[0,]


# Read and combine all training shapefiles
for (i in 1:length(f.list)){
  trn <- st_read(paste0(tdir, f.list[i])) %>% 
    mutate(class = str_split_fixed(f.list[i], "\\.", 2)[,1])
  if (nrow(trn) > 50){
    trn <- trn[1:50, ]
  }
  trn.all <- rbind(trn.all, trn)
}

# Plot the combined training data
plot(trn.all[,2])

# Load the DSM and ortho images
nsm <- rast(here("Transects/ndsm_transect_N02.tif"))
plot(nsm)
nsm_reprojected <- project(nsm, target_crs, method = "near")
plot(nsm)
img <- brick(here("Transects/ortho_transect_N02.tif"))[[c(1:3)]] 


ext <- extent(img)
template_raster <- raster(nrows=100, ncols=100, xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4])
crs(template_raster) <- crs(img)
res(template_raster) <- 0.10


resamp_ras <- resample(img, template_raster, method="bilinear")
#writeRaster(resampled_raster, here("t.tif"), overwrite = TRUE)

##reading mode
ranger_model <- readRDS("Models/rf.coverClas6.02")

# Print the confusion matrix
print(confusionMatrix(ranger_model))

# Calculate GLI for the image and stack with other layers
GLI <- ((2 * resamp_ras[[2]] - resamp_ras[[3]] - resamp_ras[[1]]) / (2 * resamp_ras[[2]] + resamp_ras[[3]] + resamp_ras[[1]]))
RLI <- ((2 * resamp_ras[[3]] - resamp_ras[[2]] - resamp_ras[[1]]) / (2 * resamp_ras[[3]] + resamp_ras[[2]] + resamp_ras[[1]]))
img <- stack(resamp_ras, GLI, RLI)

# Resample DSM to match the image resolution and stack with other layers
nsm <- raster(nsm_reprojected)
nsm2 <- resample(nsm, img)
stk <- stack(img, nsm2)
names(stk) <- c("blue", "green", "red", "GLI", "RLI","nsm")

# Predict using the updated RasterStack
rst <- raster::predict(stk, ranger_model)
# Plot and save the prediction result
plot(rst)
writeRaster(rst, here("Rasters_results/predict_transect_N02.tif"), overwrite = TRUE)


