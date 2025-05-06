library(lidR)
library(ggplot2)
library(raster)
library(terra)
library(sf)
library(here)
library(xfun)
library(sp)


#reading las
target_crs <- "EPSG:7850"
tdir <- here::here("Z:\\DEC\\ForestManagementPlan_2014_2023\\DATA\\Working\\SiteResistivity_Balingup\\Outputs\\LAS_FILES\\P02")
las_files <- list.files(tdir, pattern = "\\.las$", full.names = TRUE)

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
study_area <- buffer_list$P02


#LASCatalog processing engine 
#https://r-lidar.github.io/lidRbook/engine.html
ctg <- readLAScatalog(las_files)
ctg
las_check(ctg)

las_clipped <- clip_roi(ctg, study_area)
#writeLAS(las_clipped, here::here("LAS_Processing/las/las_clipped.las"), index = FALSE)
#plot(las_clipped, bg = "white",legend = TRUE)

#############################################################
tdir <- here::here("Z:\\DEC\\ForestManagementPlan_2014_2023\\DATA\\Working\\SiteResistivity_Balingup\\LAS_Processing\\las")
las_files <- list.files(tdir, pattern = "\\.las$", full.names = TRUE)
ctg2 <- readLAScatalog(las_files)
ctg2
las_check(ctg2)

las <- readLAS(file.path(tdir, "las_clipped.las"))
plot(las)
#############################################################
#dalponte2026 method  raster-based
las_n <- normalize_height(las, tin())

chm <- rasterize_canopy(las_n, 0.5, p2r(subcircle = 0.2), pkg = "terra")
ttops <- locate_trees(chm, lmf(ws = 2))
# Apply the algorithm
algo <- dalponte2016(chm, ttops, th_tree = 1.5)


las_seg <- segment_trees(las_n, algo, attribute = "IDdalponte")
crowns_dalponte <- crown_metrics(las_seg, func = NULL, attribute = "IDdalponte", geom = "concave")

las_n_no_ground <- filter_poi(las_n, Classification != 2L)
las_seg <- segment_trees(las_n_no_ground, algo, attribute = "IDdalponte")

# Calculate GLI for the image and stack with other layers
Green <- las_seg$G
Red <-  las_seg$R
Blue <- las_seg$B

GLI <- ((2 * Green - Red - Blue) / (2 * Green + Red +Blue ))
las_seg <- add_attribute(las_seg, GLI, "GLI")
#check attribute added
print(las_seg@data$GLI)





crowns <- crown_metrics(las, func = .stdtreemetrics,output ='sf')
crowns_2d <- st_zm(crowns, drop = TRUE, what = "ZM")


df.las <- as.data.frame(payload(las_seg))
crowns_dalponte <- crown_metrics(las_seg, func = NULL, attribute = "IDdalponte", geom = "concave")
st_write(crowns_2d, here("LAS_Processing/shp_points/Points_P03.shp"), delete_layer = TRUE)


st_write(crowns_dalponte, here("LAS_Processing/shp_points/Polygons_P02.shp"), delete_layer = TRUE)


################################################################
#for LATER
#Segmenting crowns and getting points of each tree
plot(crowns$lyr.1)
plot(las_n_no_ground)

tree1 <- filter_poi(las_seg, IDdalponte == 212)
plot(tree1, size = 8, bg = "white")

df.las <- as.data.frame(payload(tree1)) 
my_sf <- st_as_sf(df.las, coords = c("X", "Y"), crs = target_crs)
st_write(my_sf, here("LAS_Processing/shp_points/Points_tree1.shp"), delete_layer = TRUE)

ggplot(payload(tree1), aes(X,Z, color = Z)) + 
  geom_point(size = 0.5) + 
  coord_equal() + 
  theme_minimal() +
  scale_color_gradientn(colours = height.colors(50))

################################################################

#converting data into a dataframe
#method 1
df.las2 <- as.data.frame(las@data)
#method 2
df.las <- as.data.frame(payload(las))



#dont run it 
#my_sf <- st_as_sf(df.las, coords = c("X", "Y"), crs = target_crs)
#st_write(my_sf, here("LAS_Processing/shp_points/Points_P02.shp"), delete_layer = TRUE)
################################################################
