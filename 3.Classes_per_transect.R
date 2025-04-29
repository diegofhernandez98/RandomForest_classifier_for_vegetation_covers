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

#Reading buffer of transect
tdir2 <- here("shp\\")
f.list2 <-list.files(tdir2, pattern = "buffer.*\\.shp$", full.names = TRUE)


buffer_list <- list()

for (i in seq_along(f.list2)) {
  buffer_zone <- st_read(f.list2[i]) %>%
    st_transform(crs = target_crs)
  
  name_site <- buffer_zone$Site
  
  buffer_list[[i]] <- buffer_zone
  names(buffer_list)[i] <- name_site
}

intersection <- list()
class_list <- list()
for (i in seq_along(buffer_list)) {
  selection <- st_intersection(trn.all, buffer_list[[i]]) %>%
    dplyr::select(class)
  
  name_site <- names(buffer_list)[i]
  intersection[[i]] <- selection
  names(intersection)[i] <- name_site

  
  print(unique(intersection[[i]]$class))
  class_list[[i]] <- unique(intersection[[i]]$class)
  names(class_list)[i] <- name_site
  
}

df_t <- data.frame(class = character(), transect = character(), stringsAsFactors = FALSE)

for (i in seq_along(class_list)) {
  temp_df <- data.frame(
    class = class_list[[i]],
    transect = names(class_list[i]),
    stringsAsFactors = FALSE
  )
  df_t <- rbind(df_t, temp_df)
}

print(df_t)
write.csv(df_t,'Dataframe_Classes.csv')


# data frame 1 
# df <- do.call(rbind, lapply(class_list, function(x) {
#   length(x) <- max(sapply(class_list, length))
#   return(x)
# }))
# 
# df <- as.data.frame(df, stringsAsFactors = FALSE)
# t_df <- as.data.frame(t(df))
# colnames(t_df) <-names(intersection)
# write.csv(t_df,'Dataframe_Classes.csv')

