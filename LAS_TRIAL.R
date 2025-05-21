
#Indivitual tree dectection and segmentation DEFAULT MODEL
#https://r-lidar.github.io/lidRbook/itd.html
chm <- rasterize_canopy(las_clipped, 0.5, pitfree(subcircle = 0.2))

#The number of detected trees is correlated to the ws argument. 
#Small windows sizes usually gives more trees,  
#while large windows size generally miss smaller trees that are “hidden” by big trees
ttops <- locate_trees(las_clipped, lmf(ws = 3))

#plotting the trees points on top of the las file
plot(chm, col = height.colors(50))
plot(sf::st_geometry(ttops), add = TRUE, pch = 3)

x <- plot(las_clipped, bg = "white", size = 4)
add_treetops3d(x, ttops)

#Segmentation of the point-cloud by using lmf and chm default parameters 
algo <- dalponte2016(chm, ttops)
las <- segment_trees(las_clipped, algo) # segment point cloud
plot(las, bg = "white", size = 4, color = "treeID") # visualize trees

tree110 <- filter_poi(las, treeID == 1)
plot(tree110, size = 8, bg = "white")


#Points
crowns <- crown_metrics(las, func = .stdtreemetrics,output ='sf')
plot(crowns["convhull_area"], main = "Crown points (convex hull)")
crowns_2d <- st_zm(crowns, drop = TRUE, what = "ZM")
# Write the 2D points to a shapefile
st_write(crowns_2d, here("LAS_Processing/shp_points/Points_P03.shp"), delete_layer = TRUE)

#polygons
crowns2 <- crown_metrics(las, func = .stdtreemetrics, geom = "convex")
st_write(crowns2, here("LAS_Processing/shp_points/Polygons_P03.shp"), delete_layer = TRUE)

################################################################
chm <- rasterize_canopy(las_clipped, 0.5, pitfree(subcircle = 0.2))

ttops <- locate_trees(las_clipped, lmf(ws = 2))

algo <- dalponte2016(chm, ttops)
las <- segment_trees(las_clipped, algo) # segment point cloud

crowns_dalponte  <- crown_metrics(las, func = .stdtreemetrics, geom = "concave")
st_write(crowns_dalponte, here("LAS_Processing/shp_points/Polygons_P04.shp"), delete_layer = TRUE)



crowns_dalponte <- crown_metrics(las, func = NULL, attribute = "IDdalponte", geom = "concave")
crowns <- crown_metrics(las, func = .stdtreemetrics,output ='sf')

################################################################
# Canopy height model BEST MODEL SO FAR 
#las_n <- normalize_height(las_clipped, tin())

chm_p2r_05 <- rasterize_canopy(las_clipped, 0.5, p2r(subcircle = 0.2), pkg = "terra")

ttops_chm_p2r_05_smoothed <- locate_trees(chm_p2r_05, lmf(ws = 2))

# Apply the algorithm
algo <- dalponte2016(chm_p2r_05, ttops_chm_p2r_05_smoothed, th_tree = 1.5)
#algo <- dalponte2016(chm_p2r_05_smoothed, ttops_chm_p2r_05_smoothed, th_tree = 1.5, th_seed = 0.45, th_cr = 0.55)


las <- segment_trees(las_clipped, algo, attribute = "IDdalponte")
plot(las, size = 8, bg = "white")
plot(las, bg = "white", size = 4, color = "treeID", height = TRUE, treetops = TRUE)

crowns_dalponte <- crown_metrics(las, func = NULL, attribute = "IDdalponte", geom = "concave")
st_write(crowns_dalponte, here("LAS_Processing/shp_points/Polygons3_P01.shp"), delete_layer = TRUE)

tree1 <- filter_poi(las, IDdalponte == 804)
plot(tree1, size = 8, bg = "white")

tree2 <- filter_poi(las, IDdalponte == 796)
plot(tree2, size = 8, bg = "white")

################################################################

#watershed using ForestTools
library(ForestTools)

las_n <- normalize_height(las_clipped, tin())
chm <- rasterize_canopy(las_n, res = 0.5, p2r(subcircle = 0.2))
#clean up na
chm <- app(chm, fun = function(x) ifelse(x < 0, NA, x))
chm <- project(chm, target_crs)
chm <- raster::raster(chm)
#vwf_function <- function(x) { 0.1 * x + 0.5 }
vwf_function <- function(x) { x * 0.06 + 0.5 }

tree_tops <- ForestTools::vwf(chm, winFun = vwf_function, minHeight = 2)
tree_tops_sf <- st_as_sf(tree_tops)
tree_mask <- mcws(chm, tree_tops)
crowns <- mcws(chm, treetops = tree_tops_sf, minHeight = 2)
crowns <- raster::raster(crowns)
crown_polygons <- rasterToPolygons(crowns, dissolve = TRUE)
plot(crown_polygons)

sf_crowns <- st_as_sf(crown_polygons)

st_write(sf_crowns, here("LAS_Processing/shp_points/watershed.shp"), delete_layer = TRUE)

################################################################
#silva2016
las_n <- normalize_height(las, tin())
chm <- rasterize_canopy(las_n, res = 0.5, p2r(subcircle = 0.2))



ttops <- locate_trees(chm, lmf(3))
las   <- segment_trees(las, silva2016(chm, ttops))
crowns <- silva2016(chm, ttops)
las <- segment_trees(las_n,crowns)

plot(las, color = "treeID")
cro <- crown_metrics(las,func = NULL, attribute = "treeID", geom = "concave")
st_write(cro, here("LAS_Processing/shp_points/Polygons3.shp"), delete_layer = TRUE)
plot(cro)



