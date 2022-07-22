################################################################################
########## Julian Wittische - July 2022 - Hoverfly landscape genetics ##########
################################################################################
_______________________________________________________________________________
###### ELEVATION

## Using data downloaded from Copernicus
# EU-DEM 1.1; tiles:E40N30 & E40N20; ESPG: 3035 (ETRS89, LAEA)

# We need to combine two tiles
LUX_DEM_raster_path <- "../ENV_DATA_EUROPE/EU-DEM_v1.1/eu_dem_v11_E40N20/eu_dem_v11_E40N20.TIF"
LUX_DEM <- raster(LUX_DEM_raster_path)
COL_DEM_raster_path <- "../ENV_DATA_EUROPE/EU-DEM_v1.1/eu_dem_v11_E40N30/eu_dem_v11_E40N30.TIF"
COL_DEM <- raster(COL_DEM_raster_path)
#DEMp <- projectRaster(DEM, crs=crs_good) #should already be the right crs so unnecessary here
DEMp1 <- DEM1
DEMp2 <- DEM2
# Crop to raccoon study area
e <- extend(extent(sample.localesMicro_p), c(20000, 20000))
DEMpc1 <- crop(DEMp1, e)
DEMpc2 <- crop(DEMp2, e)
DEM <- merge(DEMpc1, DEMpc2, tolerance=0)
# Checking step #
plot(DEM, col=terrain.colors(100))
points(sample.localesMicro_p)
plot(ger, add=TRUE)
plot(lux, add=TRUE)

#_______________________________________________________________________________
##### IMPERVIOUSNESS
# 
# ## Using data downloaded from Copernicus
# # Imperviousness Density (built-up & sealed soil); 2018; 10m; ; ESPG: 3035 (ETRS89, LAEA); Luxembourg & Germany; NA =255
# # Based on %
# 
# IMP_raster_file_folder_lux <- "C:/Users/Utilisateur/Desktop/Projects/ENV_DATA_EUROPE/ImperviousnessDensity_2018_10m/IMD_2018_010m_lu_03035_v020/DATA"
# IMP_raster_file_folder_ger <- "C:/Users/Utilisateur/Desktop/Projects/ENV_DATA_EUROPE/ImperviousnessDensity_2018_10m/IMD_2018_010m_de_03035_v020/DATA"
# l_lux <- list.files(IMP_raster_file_folder_lux, pattern = "\\.tif$", full.names = TRUE)
# l_ger <- list.files(IMP_raster_file_folder_ger, pattern = "\\.tif$", full.names = TRUE)
# 
# # We need to combine all tiles FOR EACH country, and for that we need a function
# 
# # 1) Let's try with raster::merge()
# # merge_rasters_folder <- function(folder_path){
# #   do.call(merge, lapply(folder_path, raster))
# # }
# #IMP_lux <- merge_rasters_folder(l_lux)
# #IMP_ger <- merge_rasters_folder(l_ger)
# # /!\ raster::merge is too slow for germany here so let's build another function /!\ #
# 
# # 2) Let's try with gdalsUtils::mosaic_rasters()
# template <- raster(extent(DEM_))
# projection(template) <- crs_good
# writeRaster(template, file="temp_IMP_lux.tif", format="GTiff", overwrite=TRUE)
# writeRaster(template, file="temp_IMP_ger.tif", format="GTiff", overwrite=TRUE)
# mosaic_rasters(gdalfile=l_lux, dst_dataset="temp_IMP_lux.tif", of="GTiff")
# mosaic_rasters(gdalfile=l_ger, dst_dataset="temp_IMP_ger.tif", of="GTiff")
# IMP <- merge(raster("temp_IMP_lux.tif"), raster("temp_IMP_ger.tif"), tolerance=0) #slow but works
# NAvalue(IMP) <- 255
# plot(IMP, col=RColorBrewer::brewer.pal(9, "Reds"))
# points(sample.localesMicro_p)
# plot(ger, add=TRUE)
# plot(lux, add=TRUE)

#_______________________________________________________________________________
# ##### TREE COVER
# 
# ## Using data downloaded from Copernicus
# # Tree Cover Density (no type info but available online if we need it); 2018; 10m; ESPG: 3035 (ETRS89, LAEA); NA = 255Luxembourg & Germany
# # Based on %
# 
# TC_raster_file_folder_lux <- "C:/Users/Utilisateur/Desktop/Projects/ENV_DATA_EUROPE/TreeCoverDensity_2018_10m/TCD_2018_010m_lu_03035_v020/DATA"
# TC_raster_file_folder_ger <- "C:/Users/Utilisateur/Desktop/Projects/ENV_DATA_EUROPE/TreeCoverDensity_2018_10m/TCD_2018_010m_de_03035_v020/DATA"
# l_lux <- list.files(TC_raster_file_folder_lux, pattern = "\\.tif$", full.names = TRUE)
# l_ger <- list.files(TC_raster_file_folder_ger, pattern = "\\.tif$", full.names = TRUE)
# 
# writeRaster(template, file="temp_TC_lux.tif", format="GTiff", overwrite=TRUE)
# writeRaster(template, file="temp_TC_ger.tif", format="GTiff", overwrite=TRUE)
# mosaic_rasters(gdalfile=l_lux, dst_dataset="temp_TC_lux.tif", of="GTiff")
# mosaic_rasters(gdalfile=l_ger, dst_dataset="temp_TC_ger.tif", of="GTiff")
# TC <- merge(raster("temp_TC_lux.tif"), raster("temp_TC_ger.tif"), tolerance=0) #slow but works
# NAvalue(TC) <- 255
# plot(TC, col=RColorBrewer::brewer.pal(9, "Greens"))
# points(sample.localesMicro_p)
# plot(ger, add=TRUE)
# plot(lux, add=TRUE)

#_______________________________________________________________________________
# ##### GRASSLAND
# 
# ## Using data downloaded from Copernicus
# # Grassland Density; 2018; 10m; ; ESPG: 3035 (ETRS89, LAEA); Luxembourg & Germany; NA = 255
# # 0: No grassland; 1: Grassland
# 
# GRA_raster_file_folder_lux <- "C:/Users/Utilisateur/Desktop/Projects/ENV_DATA_EUROPE/GrasslandDensity_2018_10m/GRA_2018_010m_lu_03035_v010/DATA"
# GRA_raster_file_folder_ger <- "C:/Users/Utilisateur/Desktop/Projects/ENV_DATA_EUROPE/GrasslandDensity_2018_10m/GRA_2018_010m_de_03035_v010/DATA"
# l_lux <- list.files(GRA_raster_file_folder_lux, pattern = "\\.tif$", full.names = TRUE)
# l_ger <- list.files(GRA_raster_file_folder_ger, pattern = "\\.tif$", full.names = TRUE)
# 
# writeRaster(template, file="temp_GRA_lux.tif", format="GTiff", overwrite=TRUE)
# writeRaster(template, file="temp_GRA_ger.tif", format="GTiff", overwrite=TRUE)
# mosaic_rasters(gdalfile=l_lux, dst_dataset="temp_GRA_lux.tif", of="GTiff")
# mosaic_rasters(gdalfile=l_ger, dst_dataset="temp_GRA_ger.tif", of="GTiff")
# GRA <- merge(raster("temp_GRA_lux.tif"), raster("temp_GRA_ger.tif"), tolerance=0) #slow but works
# NAvalue(GRA) <- 255
# plot(GRA, col=RColorBrewer::brewer.pal(9, "Greys"))
# points(sample.localesMicro_p)
# plot(ger, add=TRUE)
# plot(lux, add=TRUE)

#_______________________________________________________________________________
# ##### WATER AND WETNESS
# 
# ## Using data downloaded from Copernicus
# # Occurrence of water and wet surfaces; 2018; 10m; ; ESPG: 3035 (ETRS89, LAEA); Luxembourg & Germany; NA = 255, 254 (unclassifiable), 253 (seawater)
# # 0: Dry; 1: Permanent water; 2: Temporary water; 3: Permanent wet; 4: Temporary wet
# 
# WAW_raster_file_folder_lux <- "C:/Users/Utilisateur/Desktop/Projects/ENV_DATA_EUROPE/Water&Wetness_2018_10m/WAW_2018_010m_lu_03035_v020/DATA"
# WAW_raster_file_folder_ger <- "C:/Users/Utilisateur/Desktop/Projects/ENV_DATA_EUROPE/Water&Wetness_2018_10m/WAW_2018_010m_de_03035_v020/DATA"
# l_lux <- list.files(WAW_raster_file_folder_lux, pattern = "\\.tif$", full.names = TRUE)
# l_ger <- list.files(WAW_raster_file_folder_ger, pattern = "\\.tif$", full.names = TRUE)
# 
# writeRaster(template, file="temp_WAW_lux.tif", format="GTiff", overwrite=TRUE)
# writeRaster(template, file="temp_WAW_ger.tif", format="GTiff", overwrite=TRUE)
# mosaic_rasters(gdalfile=l_lux, dst_dataset="temp_WAW_lux.tif", of="GTiff")
# mosaic_rasters(gdalfile=l_ger, dst_dataset="temp_WAW_ger.tif", of="GTiff")
# WAW <- merge(raster("temp_WAW_lux.tif"), raster("temp_WAW_ger.tif"), tolerance=0) #slow but works
# WAW <- reclassify(WAW, cbind(4, Inf, NA), right=TRUE)
# plot(WAW, col=RColorBrewer::brewer.pal(9, "Blues"))
# points(sample.localesMicro_p)
# plot(ger, add=TRUE)
# plot(lux, add=TRUE)

## Clean-up and saving the final products
#writeRaster(DEM, file="RACCOON_DEM.grd", overwrite=TRUE)
#writeRaster(IMP, file="RACCOON_IMP.grd", format="raster", overwrite=TRUE)
#writeRaster(TC, file="RACCOON_TC.grd", format="raster", overwrite=TRUE)
#writeRaster(GRA, file="RACCOON_GRA.grd", format="raster", overwrite=TRUE)
#writeRaster(WAW, file="RACCOON_WAW.grd", format="raster", overwrite=TRUE)
#temp_remove<- list.files(pattern = "\\.*temp.*\\.tif$")
#sapply(temp_remove, unlink)

# If WAW is not enough, we can use this (not polygonized and rasterized yet)
# #_______________________________________________________________________________
# ##### RIVER NETWORK
# 
# list_rivers <- list.files("C:/Users/Utilisateur/Desktop/Projects/ENV_DATA_EUROPE/Hydro_gpkg",
#                           pattern = "v013\\.gpkg$", recursive = TRUE, full.names = TRUE)
# list_rivers2 <- lapply(list_rivers, st_read)
# # Selecting perennial and on/above surface rivers
# list_rivers3 <- list_rivers2 %>% lapply(. %>% select(HYP,LOC) %>% filter(HYP==1 & LOC>40) %>% select(HYP) %>% st_zm)
# rivers <- st_union(list_rivers3[[1]], list_rivers3[[2]], list_rivers3[[3]], list_rivers3[[4]],
#                    list_rivers3[[5]], list_rivers3[[6]])
# rivers <- rivers %>% select(HYP)
# st_write(rivers, "RACCOON_Rivers.gpkg", append=FALSE)
# 
# # We need to convert multi-line strings to polygons
# rivers_poly <- st_cast(rivers)
# rivers_poly <- st_polygonize(rivers)
# 
# # We need to rasterize in some way to be able to use it as a resistance surface.
# library(fasterize)
# fasterize(rivers_poly, DEM_)
# # Let's start with 25m just like the DEM