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
#dev.off()
for (i in seq_along(buffer_list)) {
  selection <- st_intersection(trn.all, buffer_list[[i]]) %>%
    dplyr::select(class)
  
  name_site <- names(buffer_list)[i]
  intersection[[i]] <- selection
  names(intersection)[i] <- name_site
  #plottingg each intersection
  #plot(intersection[[i]], key.pos = 1, key.width = lcm(3.34), key.size = lcm(1.5))
}

df_list <- list()
# Count occurrences of each class in each intersection
for (i in seq_along(intersection)) {
  counts <- table(intersection[[i]]$class)
  df <- as.data.frame(counts)
  name_site <- names(buffer_list)[i]
  colnames(df) <- c("class", "count")
  
  df_list[[i]] <- df
  names(df_list)[i] <- name_site
  p <- ggplot(df,aes(x = class, y = count, fill = names(df_list)[i])) +
    geom_bar(stat = "identity", position = position_dodge(df_list)) +
    labs(title = "Classes Counts", x = "Class", y = "Count")
  print(p)
}
combined_df <- bind_rows(df_list, .id = "site")

ggplot(data = combined_df, aes(x = class, y = count, fill = site)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Class Counts Across Sites", x = "Class", y = "Count") +
  theme_minimal()

# Transpose the data frame
transposed_data <- t(combined_df)
