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
nsm <- rast(here("ndsm_transect.tif"))
nsm_reprojected <- project(nsm, target_crs, method = "near")
plot(nsm)
img <- brick(here("ortho_transect.tif"))[[c(1:3)]] 



###graphs about accuracy 

#reading model
resampled_rasters <- readRDS("rasters")

## extracting spectral values and adding into a list
trn.all_list <- list()
extracted_values_list <- list()
df.trn_list <- list()

for (i in 1:length(resampled_rasters)) {
  if (st_crs(trn.all) != st_crs(nsm_reprojected)) {
    trn.all <- st_transform(trn.all, st_crs(nsm_reprojected))
  }
  trn.all_list[[i]] <- trn.all
  extracted_values <- raster::extract(nsm_reprojected, st_coordinates(trn.all_list[[i]]))
  extracted_values_list[[i]] <- extracted_values
  extracted_values <- unlist(extracted_values_list[[i]])
  trn.all$nsm <- extracted_values
  trn.all_list[[i]] <- trn.all
  
  df.trn <- st_drop_geometry(trn.all_list[[i]]) %>%
    dplyr::select(-Id)
  
  df.trn_list[[i]] <- df.trn
  
}


df.img_list <- list()
dfg_list <- list()
models <- list()
for (i in 1:length(resampled_rasters)) {
  df.img <- as.data.frame(raster::extract(resampled_rasters[[i]], trn.all_list[[i]]))
  colnames(df.img) <- c("blue", "green", "red")  # Set column names here
  df.img_list[[i]] <- df.img
  
  df.trn <- bind_cols(df.trn_list[[i]], df.img_list[[i]])
  df.trn <- df.trn %>% na.omit()
  
  df.trn <- df.trn %>% mutate(
    GLI = (2 * green - red - blue) / (2 * green + red + blue),
    RLI = (2 * red - green - blue) / (2 * red + green + blue)
  )  
  df.trn_list[[i]] <- df.trn
  
  
  # Plot value boxplots
  dfg <- df.trn_list[[i]] %>%
    pivot_longer(cols = 2:7, names_to = "measure", values_to = "value")
  dfg$measure <- factor(dfg$measure, levels = c("blue", "green", "red", "GLI", "RLI", "nsm"))
  dfg_list[[i]] <- dfg
  
  # Create and save the boxplot
  p <- ggplot(dfg, aes(x = class, y = value)) +
    geom_boxplot() + 
    facet_wrap(~measure, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  #ggsave(here("graphs", paste0("site027_boxplot_", i, ".jpg")), plot = p, height = 6, width = 6)
  
  mtry <- seq(from = 2, to = ncol(df.trn_list[[i]]) - 2, by = 1)
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
  models[[i]] <- ranger_model
}

#creating data frame with accuracies
ac <- list()
for (i in 1:length(models)){
  resolutions <- seq(0.1, 1, by = 0.05)
  acc <- max(models[[i]]$results$Accuracy)
  ac[[i]] <- acc
}
df_acc <- data.frame(resolutions, accuracy = unlist(ac))
print(df_acc)

p2 <- ggplot(df_acc, aes (x = resolutions, y = accuracy )) +
  geom_line() +
  geom_point() +
  labs(title = "Model Accuracy", x = "Raster Resolution", y = "RF Accuracy") +
  theme_minimal()
ggsave(here("graphs", paste0("accuracy_graph.jpg")), plot = p2, height = 6, width = 6)

print(p2)


##run up to here


############################################################################################################################################################


###resampling loop
resolutions <- seq(0.1, 1, by = 0.05)
print(resolutions)

resampled_rasters <- list()
for (i in 1:length(resolutions)) {
  new_res <- resolutions[i]
  ext <- extent(img)
  template_raster <- raster(nrows=100, ncols=100, xmn=ext[1], xmx=ext[2], ymn=ext[3], ymx=ext[4])
  crs(template_raster) <- crs(img)
  res(template_raster) <- new_res
  resampled_rasters[[i]] <- resample(img, template_raster, method = "bilinear")
  #writeRaster(resampled_rasters[[i]], here::here(paste0("resampled_raster_", new_res, ".tif")), format = "GTiff")
}

#create a temporary file path for the entire list
temp_file <- tempfile(pattern = "resampled_rasters", fileext = ".rds")
#save the entire list to the temporary file
saveRDS(resampled_rasters, here::here("rasters"))
############################################################################################################################################################