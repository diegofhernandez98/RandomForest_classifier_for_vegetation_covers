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
  mutate(class = str_split_fixed(f.list[1], "\\.", 2)[,1])%>%
  st_transform(crs = target_crs)

trn.all <- trn.all[0,]

# Read and combine all training shapefiles
for (i in 1:length(f.list)){
  trn <- st_read(paste0(tdir, f.list[i])) %>% 
    mutate(class = str_split_fixed(f.list[i], "\\.", 2)[,1]) %>%
    st_transform(crs = target_crs)
  
  trn.all <- rbind(trn.all, trn) 
  }
  
# Plot the combined training data
#plot(trn.all[,2])


#Reading buffer of transect
buffer_zone <- st_read("shp/bufferP21.shp") %>%
  st_transform(crs = target_crs)
name_site<- unlist(buffer_zone$Site)
#plot(buffer_zone)

#filtering trn.all
selection <- st_intersection(trn.all, buffer_zone)  %>%
  dplyr::select(class)

#plot(selection)

trn.all <- selection

plot(trn.all[,2])

# Load the DSM and ortho images
nsm <- rast(here("Transects/ndsm_transect_P21.tif"))
#plot(nsm)
nsm_reprojected <- project(nsm, target_crs, method = "near")
#plot(nsm)
img <- brick(here("Transects/ortho_transect_P21.tif"))[[c(1:3)]] 

ext <- extent(img)
template_raster <- raster(nrows=100, ncols=100, xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4])
crs(template_raster) <- crs(img)
res(template_raster) <- 0.10


resamp_ras <- resample(img, template_raster, method="bilinear")
resamp_projected <- projectRaster(resamp_ras, crs = target_crs, method = "ngb")
#writeRaster(resamp_projected, filename = here("Rasters_results", paste0("predict_transect_N04.tif")), overwrite = TRUE)
#writeRaster(resampled_raster, here("t.tif"), overwrite = TRUE)


# Extract DSM values and add to the training data frame
#####PREVIOUS METHOD###################################
#Ensure CRS match
if (st_crs(trn.all) != crs(nsm_reprojected)) {
  trn.all <- st_transform(trn.all, crs(nsm_reprojected))
}
# df.trn$nsm <- raster::extract(nsm, trn.all)
extracted_values<- raster::extract(nsm_reprojected, st_coordinates(trn.all))

# Ensure the extracted values are a simple vector
extracted_values <- unlist(extracted_values)

trn.all$nsm <- extracted_values

head(trn.all)

df.trn <- st_drop_geometry(trn.all)
###################################
###buffer 
# library(exactextractr)
# 
# #zonal statistics
# buffer_trn.all <- st_buffer(trn.all, dist = 5)
# buffer_trn.all$nsm_mean <- exact_extract(nsm, buffer_trn.all, 'mean')
# 
# df.trn <- st_drop_geometry(buffer_trn.all) %>%
#   dplyr::select(-Id)
###################################
# Extract image values and add to the training data frame
df.img <- as.data.frame(raster::extract(resamp_projected, trn.all))
colnames(df.img) <- c("blue", "green", "red")
df.trn <- bind_cols(df.trn, df.img)

# Remove rows with NA values
df.trn <- df.trn %>% na.omit()

# Calculate GLI and add to the training data frame
df.trn <- mutate(df.trn, GLI = (2 * green - red - blue) / (2 * green + red + blue), RLI =(2 * red - green - blue) / (2 * red + green + blue) )

# df.trn <- mutate(df.trn, GLI = (2 * green - red - blue) / (2 * green + red + blue))
# df.trn <- mutate(df.trn, RLI =(2 * red - green - blue) / (2 * red + green + blue))


# Plot value boxplots
dfg <- gather(df.trn, measure, value, 2:7)
dfg$measure <- factor(dfg$measure, levels = c("blue", "green", "red","GLI","RLI", "nsm"))
ggplot(dfg, aes(x = class, y = value)) +
  geom_boxplot() + 
  facet_wrap(~measure, scales = "free") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = here("graphs", paste0("site", name_site, "_boxplot.jpg")),height = 6, width = 6)


# Set up the random forest model
mtry <- seq(from = 2, to = ncol(df.trn) - 2, by = 1)
tunegrid <- expand.grid(.mtry = mtry,
                        .splitrule = c("extratrees"),
                        .min.node.size = 1)

ranger_model <- train(
  class ~ .,
  tuneGrid = tunegrid,
  importance = "permutation",
  data = df.trn, 
  method = "ranger",
  trControl = trainControl(method = "cv", 
                           number = 5,
                           verboseIter = FALSE))

# Save the model
saveRDS(ranger_model, file = here("Models", paste0("rf.coverClassName_site", name_site, ".02")))

# Print the confusion matrix
print(confusionMatrix(ranger_model))

# Calculate GLI for the image and stack with other layers
GLI <- ((2 * resamp_projected[[2]] - resamp_projected[[3]] - resamp_projected[[1]]) / (2 * resamp_projected[[2]] + resamp_projected[[3]] + resamp_projected[[1]]))
RLI <- ((2 * resamp_projected[[3]] - resamp_projected[[2]] - resamp_projected[[1]]) / (2 * resamp_projected[[3]] + resamp_projected[[2]] + resamp_projected[[1]]))
img <- stack(resamp_projected, GLI, RLI)

# Resample DSM to match the image resolution and stack with other layers
nsm <- raster(nsm_reprojected)
nsm2 <- resample(nsm, img)
stk <- stack(img, nsm2)
names(stk) <- c("blue", "green", "red", "GLI", "RLI","nsm")

# Predict using the updated RasterStack
rst <- raster::predict(stk, ranger_model)

writeRaster(rst, filename = here("Rasters_results", paste0("predict_transect", name_site, ".tif")), overwrite = TRUE)
