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

library(doParallel)
library(foreach)

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


tdir2 <- here("shp")
buffer_list <- list.files(tdir2, pattern = "\\.shp$", full.names = FALSE)

buffer_e_list <- list()
for (i in seq_along(buffer_list)) {
  file_path <- file.path(tdir2, buffer_list[i])
  var_name <- tools::file_path_sans_ext(buffer_list[i])
  
  buffer_zone <- st_read(file_path)%>%
    st_transform(crs = target_crs)
  name_site<- unlist(buffer_zone$Site)
  buffer_e_list[[name_site]] <- buffer_zone
}
#to review
selection_list <- list()
for (i in seq_along(buffer_e_list)){
  selection <- st_intersection(trn.all, buffer_e_list[[i]])  %>%
    dplyr::select(class)
  selection$name_site <- buffer_e_list[[i]]$Site
  selection_list <- selection
}


tdir3 <-here("Transects")
ndsm_l <- list.files(tdir3, pattern = "ndsm_.*\\.tif$", full.names = FALSE)

ndsm_list <- list()
for(i in seq_along(ndsm_l)){
  file_path <- file.path(tdir3, ndsm_l[i])
  ndsm <- rast(file_path)
  ndsm_reprojected <- project(ndsm, target_crs, method = "near")
  ndsm_list[[i]] <- ndsm_reprojected
}

#set number of cores and register parallel backend
num_cores <- 10
cluster_obj <- makeCluster(num_cores)
registerDoParallel(cluster_obj)

transect_dir <- here("Transects")
ortho_files <- list.files(transect_dir, pattern = "ortho_.*\\.tif$", full.names = FALSE)

## https://gis.stackexchange.com/questions/434317/how-to-speed-up-the-intersection-of-two-large-shapefiles-in-r-with-sf-or-any-ot
#https://cran.r-project.org/web/packages/doParallel/vignettes/gettingstartedParallel.pdf
#run parallel processing
ortho_outputs <- foreach(i = seq_along(ortho_files),
                         .packages = c("raster", "terra"),
                         .export = c("target_crs", "transect_dir", "ortho_files")) %dopar% {
                           tryCatch({
                             raster_path <- file.path(transect_dir, ortho_files[i])
                             ortho_raster <- brick(raster_path)[[1:3]]
                             
                             ortho_extent <- extent(ortho_raster)
                             template_grid <- raster(nrows = 100, ncols = 100,
                                                     xmn = xmin(ortho_extent), xmx = xmax(ortho_extent),
                                                     ymn = ymin(ortho_extent), ymx = ymax(ortho_extent))
                             crs(template_grid) <- crs(ortho_raster)
                             res(template_grid) <- 0.10
                             
                             resampled_raster <- resample(ortho_raster, template_grid, method = "bilinear")
                             projected_raster <- projectRaster(resampled_raster, crs = target_crs, method = "ngb")
                             
                             list(original = ortho_raster, template = template_grid, projected = projected_raster)
                           }, error = function(e) {
                             message(sprintf("Error processing file %s: %s", ortho_files[i], e$message))
                             NULL
                           })
                         }

stopCluster(cluster_obj)

ortho_outputs <- Filter(Negate(is.null), ortho_outputs)

# Unpack results into separate lists
ortho_rasters <- lapply(ortho_outputs, `[[`, "original")
template_grids <- lapply(ortho_outputs, `[[`, "template")
projected_rasters <- lapply(ortho_outputs, `[[`, "projected")
###########################################################################################
#running without parallel processing

# ortho_l <- list.files(tdir3, pattern = "ortho_.*\\.tif$", full.names = FALSE)
# 
# ortho_list <- list()
# template_list <- list()
# resamp_list <- list()
# 
# for (i in seq_along(ortho_l)) {
#   file_path <- file.path(tdir3, ortho_l[i])
#   img <- brick(here(file_path))[[1:3]] 
#   ortho_list[[i]] <- img
#   
#   ext <- extent(img)
#   template_raster <- raster(nrows = 100, ncols = 100, 
#                             xmn = xmin(ext), xmx = xmax(ext), 
#                             ymn = ymin(ext), ymx = ymax(ext))
#   crs(template_raster) <- crs(img)
#   res(template_raster) <- 0.10
#   template_list[[i]] <- template_raster
#   
#   resamp_ras <- resample(img, template_raster, method = "bilinear")
#   resamp_projected <- projectRaster(resamp_ras, crs = target_crs, method = "ngb")
#   resamp_list[[i]] <- resamp_projected
# }
###########################################################################################

for (i in seq_along(ndsm_list)){
  if (st_crs(trn.all) != crs(ndsm_list[[1]])) {
    trn.all <- st_transform(trn.all, crs(ndsm_list[[1]]))
  }
}

###
#Creating dataframe extracting the data values
df.trn_list <- list()

for (i in seq_along(ndsm_list)) {
  # Extract nDSM values
  extracted_values <- raster::extract(ndsm_list[[i]], st_coordinates(trn.all))
  extracted_values <- unlist(extracted_values)
  
  trn_copy <- trn.all
  trn_copy$nsm <- extracted_values
  df.trn_base <- st_drop_geometry(trn_copy)
  
  for (j in seq_along(projected_rasters)) {
    df.img <- as.data.frame(raster::extract(projected_rasters[[j]], trn_copy))
    colnames(df.img) <- c("blue", "green", "red")
    
    df.trn <- bind_cols(df.trn_base, df.img)
    df.trn <- df.trn %>% na.omit()
    
    df.trn <- df.trn %>%
      mutate(
        GLI = (2 * green - red - blue) / (2 * green + red + blue),
        RLI = (2 * red - green - blue) / (2 * red + green + blue)
      )
    
    df.trn_list[[paste0("ndsm", i, "_img", j)]] <- df.trn
  }
}

###
#Plotting boxplot
df.trn_list <- df.trn_list[sapply(df.trn_list, function(df) nrow(df) > 0)]
names(df.trn_list) <- names(buffer_e_list)

plot_vars <- c("blue", "green", "red", "GLI", "RLI", "nsm")

for (i in seq_along(df.trn_list)) {
  df <- df.trn_list[[i]]
  df_name <- names(df.trn_list)[i]
  
  # Ensure required columns exist
  if (!all(c("class", plot_vars) %in% colnames(df))) next
  
  # Reshape for plotting
  dfg <- tidyr::pivot_longer(df, cols = all_of(plot_vars), names_to = "measure", values_to = "value")
  dfg$measure <- factor(dfg$measure, levels = plot_vars)
  
  # Plot
  p <- ggplot(dfg, aes(x = class, y = value)) +
    geom_boxplot() +
    facet_wrap(~measure, scales = "free") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle(paste("Boxplot for", df_name))
  
  # Save
  ggsave(
    filename = here("graphs", paste0("site_", df_name, "_boxplot.jpg")),
    plot = p,
    height = 6,
    width = 8
  )
}

stopImplicitCluster()

###
#Training model 

# Create directory for saving models if it doesn't exist
dir.create(here("Models"), showWarnings = FALSE)

model_list <- list()
for (i in seq_along(df.trn_list)) {
  df.trn <- df.trn_list[[i]]
  name_site <- names(df.trn_list)[i]
  
  df.trn$class <- as.factor(df.trn$class)
  
  df.trn <- dplyr::select(df.trn, -Id)
  
  predictor_cols <- setdiff(names(df.trn), "class")
  mtry_vals <- seq(from = 2, to = length(predictor_cols), by = 1)

  tunegrid <- expand.grid(
    .mtry = mtry_vals,
    .splitrule = "extratrees",
    .min.node.size = 1
  )
  
  # Train model
  ranger_model <- train(
    class ~ .,
    data = df.trn,
    method = "ranger",
    importance = "permutation",
    tuneGrid = tunegrid,
    trControl = trainControl(
      method = "cv",
      number = 5,
      allowParallel = FALSE
    )
  )
  
  # Save model
  saveRDS(ranger_model, file = here("Models", paste0("rf.coverClassName_site_", name_site, ".rds")))
  print(confusionMatrix(ranger_model))
  model_list[[i]] <- ranger_model
}
ras <- projected_rasters[[3]] ##here
GLI <- ((2 * ras[[2]] - ras[[3]] - ras[[1]]) / (2 * ras[[2]] + ras[[3]] + ras[[1]]))
RLI <- ((2 * ras[[3]] - ras[[2]] - ras[[1]]) / (2 * ras[[3]] + ras[[2]] + ras[[1]]))
img <- stack(ras, GLI, RLI)

ndsm <- raster::raster(ndsm_list[[3]]) ##here
if (!is.null(intersect(extent(ndsm), extent(img)))) {
  ndsm2 <- resample(ndsm, img)
  stk <- stack(img, ndsm2)
  names(stk) <- c("blue", "green", "red", "GLI", "RLI", "nsm")
  
  ranger_model <- model_list[[3]] ##here
  site_name <- names(df.trn_list)[3] ##here
  
  rst <- raster::predict(stk, ranger_model, progress = "text")
  writeRaster(rst, here("Rasters_results", paste0("prediction_site_", site_name, ".tif")),
              overwrite = TRUE, options = c("COMPRESS=LZW"))
} else {
  message("No spatial overlap between raster and nDSM.")
}

n_items <- length(projected_rasters)

for (i in seq_len(n_items)) {
  ras <- projected_rasters[[i]]
  GLI <- ((2 * ras[[2]] - ras[[3]] - ras[[1]]) / (2 * ras[[2]] + ras[[3]] + ras[[1]]))
  RLI <- ((2 * ras[[3]] - ras[[2]] - ras[[1]]) / (2 * ras[[3]] + ras[[2]] + ras[[1]]))
  img <- stack(ras, GLI, RLI)
  
  ndsm <- raster::raster(ndsm_list[[i]])
  
  if (!is.null(intersect(extent(ndsm), extent(img)))) {
    ndsm2 <- resample(ndsm, img)
    stk <- stack(img, ndsm2)
    names(stk) <- c("blue", "green", "red", "GLI", "RLI", "nsm")
    
    ranger_model <- model_list[[i]]
    site_name <- names(df.trn_list)[i]
    
    rst <- raster::predict(stk, ranger_model, progress = "text")
    writeRaster(rst, here("Rasters_results", paste0("prediction_site_", site_name, ".tif")),
                overwrite = TRUE, options = c("COMPRESS=LZW"))
  } else {
    message(paste("No spatial overlap between raster and nDSM for index", i))
  }
}




