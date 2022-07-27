################################################################################
########## Julian Wittische - July 2022 - Hoverfly landscape genetics ##########
################################################################################
#
# SCRIPT OBJECTIVE:
#
# - creating GIS rasters for relevant environmental variables
#   which are readily usable in ResitanceGA or other softwares
#_______________________________________________________________________________
# TO DO:
#
# - check that no M. florea is out of the area defined by the extent of the SP
#   dataset + 2km
#
# - check that there are no other discrepancies between online metadata file and
#   the files used here
#
# - check if there are no observations in the water
#
# - choose if we keep WAW simplification, and do it
#
# - raster to polygons to extract river from WAW
#
# ISSUES:
#
#  - there was an outlier for LUX, it should be 49.52180 not 49.4218 for sample 
#    SW_109_2_27-07-21_A; I corrected it based on original metada 
#_______________________________________________________________________________
# ###### LIBRARIES AND DATA
# 
# library(raster)
# library(sp)
# library(sf)
# library(gdalUtils)
# library(dplyr)
# 
# # We can derive the extent from the genind files - let's load the SP object
# SP <- readRDS("Data/SP_genind_noSpp141_360_080.RDS")
# SP_LUX <- SP[SP@pop=="SW"|SP@pop=="LU"]
# SP_COL <- SP[SP@pop=="CO"]
# 
# # Study are extent + 2km for cropping/masking
# e_LUX <- extend(extent(SP_LUX@other$xy), c(2000, 2000))
# e_COL <- extend(extent(SP_COL@other$xy), c(2000, 2000))
# 
#_______________________________________________________________________________
# ###### ELEVATION
# 
# ## Using data downloaded from Copernicus
# # EU-DEM 1.1; tiles:E40N30 & E40N20; ESPG: 3035 (ETRS89, LAEA)
# 
# # We need to combine two tiles
# LUX_DEM_raster_path <- "../ENV_DATA_EUROPE/EU-DEM_v1.1/eu_dem_v11_E40N20/eu_dem_v11_E40N20.TIF"
# LUX_DEM <- raster(LUX_DEM_raster_path)
# COL_DEM_raster_path <- "../ENV_DATA_EUROPE/EU-DEM_v1.1/eu_dem_v11_E40N30/eu_dem_v11_E40N30.TIF"
# COL_DEM <- raster(COL_DEM_raster_path)
# 
# # We keep the crs from the various rasters from Copernicus (3035)
# crs_good <- CRS(SRS_string = "EPSG:3035")
# 
# # It should already be the right crs but let's make sure
# LUX_DEMp <- projectRaster(LUX_DEM, crs=crs_good)
# COL_DEMp <- projectRaster(COL_DEM, crs=crs_good)
# 
# # Checkpoint
# st_crs(LUX_DEM)
# st_crs(LUX_DEMp)
# st_crs(COL_DEM)
# st_crs(COL_DEMp)
# 
# # Crop to study area
# LUX_DEMpc <- crop(LUX_DEM, e_LUX)
# COL_DEMpc <- crop(COL_DEM, e_COL)
# 
# # Checkpoint
# plot(LUX_DEMpc, col=terrain.colors(100))
# points(SP_LUX@other$xy)
# plot(COL_DEMpc, col=terrain.colors(100))
# points(SP_COL@other$xy)

# # Clean-up and save the final products
# writeRaster(LUX_DEMpc, file="Data/LUX_DEM.grd", overwrite=TRUE)
# writeRaster(COL_DEMpc, file="Data/COL_DEM.grd", overwrite=TRUE)
# remove(LUX_DEM); remove(COL_DEM); remove(LUX_DEMpc); remove(COL_DEMpc)

#_______________________________________________________________________________
# #### IMPERVIOUSNESS
# 
# ## Using data downloaded from Copernicus
# # Imperviousness Density (built-up & sealed soil); 2018; 10m; ; ESPG: 3035 (ETRS89, LAEA); Luxembourg & Germany; NA =255
# # Based on %
# 
# LUX_IMP_raster_file_folder <- "../ENV_DATA_EUROPE/IMD2018/IMD_2018_010m_lu_03035_v020/DATA"
# COL_IMP_raster_file_folder <- "../ENV_DATA_EUROPE/IMD2018/IMD_2018_010m_de_03035_v020/DATA"
# l_LUX <- list.files(LUX_IMP_raster_file_folder, pattern = "\\.tif$", full.names = TRUE)
# l_COL <- list.files(COL_IMP_raster_file_folder, pattern = "\\.tif$", full.names = TRUE)
# 
# template_LUX <- raster(extent(LUX_DEMpc))
# projection(template_LUX) <- crs_good
# template_COL <- raster(extent(COL_DEMpc))
# projection(template_COL) <- crs_good
# writeRaster(template_LUX, file="temp_IMP_LUX.tif", format="GTiff", overwrite=TRUE)
# writeRaster(template_COL, file="temp_IMP_COL.tif", format="GTiff", overwrite=TRUE)
# mosaic_rasters(gdalfile=l_LUX, dst_dataset="temp_IMP_LUX.tif", of="GTiff")
# mosaic_rasters(gdalfile=l_COL, dst_dataset="temp_IMP_COL.tif", of="GTiff")
# LUX_IMP <- raster("temp_IMP_LUX.tif")
# LUX_IMPpc <- crop(LUX_IMP, e_LUX)
# NAvalue(LUX_IMPpc) <- 255
# COL_IMP <- raster("temp_IMP_COL.tif")
# COL_IMPpc <- crop(COL_IMP, e_COL)
# NAvalue(COL_IMPpc) <- 255
# 
# # Checkpoint
# plot(LUX_IMPpc, col=RColorBrewer::brewer.pal(9, "Reds"))
# points(SP_LUX@other$xy)
# plot(COL_IMPpc, col=RColorBrewer::brewer.pal(9, "Reds"))
# points(SP_COL@other$xy)
# 
# # Clean-up and save the final products
# writeRaster(LUX_IMPpc, file="Data/LUX_IMP.grd", overwrite=TRUE)
# writeRaster(COL_IMPpc, file="Data/COL_IMP.grd", overwrite=TRUE)
# remove(LUX_IMP); remove(COL_IMP); remove(LUX_IMPpc); remove(COL_IMPpc)
# temp_remove<- list.files(pattern = "\\.*temp.*\\.tif$")
# sapply(temp_remove, unlink)


#_______________________________________________________________________________
# ##### TREE COVER
# 
# ## Using data downloaded from Copernicus
# # Tree Cover Density (no type info but available online if we need it); 2018; 10m; ESPG: 3035 (ETRS89, LAEA); NA = 255Luxembourg & Germany
# # Based on %
# 
# LUX_TCD_raster_file_folder <- "../ENV_DATA_EUROPE/TCD2018/TCD_2018_010m_lu_03035_v020/DATA"
# COL_TCD_raster_file_folder <- "../ENV_DATA_EUROPE/TCD2018/TCD_2018_010m_de_03035_v020/DATA"
# l_LUX <- list.files(LUX_TCD_raster_file_folder, pattern = "\\.tif$", full.names = TRUE)
# l_COL <- list.files(COL_TCD_raster_file_folder, pattern = "\\.tif$", full.names = TRUE)
# 
# template_LUX <- raster(extent(LUX_DEMpc))
# projection(template_LUX) <- crs_good
# template_COL <- raster(extent(COL_DEMpc))
# projection(template_COL) <- crs_good
# writeRaster(template_LUX, file="temp_TCD_LUX.tif", format="GTiff", overwrite=TRUE)
# writeRaster(template_COL, file="temp_TCD_COL.tif", format="GTiff", overwrite=TRUE)
# mosaic_rasters(gdalfile=l_LUX, dst_dataset="temp_TCD_LUX.tif", of="GTiff")
# mosaic_rasters(gdalfile=l_COL, dst_dataset="temp_TCD_COL.tif", of="GTiff")
# LUX_TCD <- raster("temp_TCD_LUX.tif")
# LUX_TCDpc <- crop(LUX_TCD, e_LUX)
# NAvalue(LUX_TCDpc) <- 255
# COL_TCD <- raster("temp_TCD_COL.tif")
# COL_TCDpc <- crop(COL_TCD, e_COL)
# NAvalue(COL_TCDpc) <- 255
# 
# # Checkpoint
# plot(LUX_TCDpc, col=RColorBrewer::brewer.pal(9, "Greens"))
# points(SP_LUX@other$xy)
# plot(COL_TCDpc, col=RColorBrewer::brewer.pal(9, "Greens"))
# points(SP_COL@other$xy)
# 
# # Clean-up and save the final products
# writeRaster(LUX_TCDpc, file="Data/LUX_TCD.grd", overwrite=TRUE)
# writeRaster(COL_TCDpc, file="Data/COL_TCD.grd", overwrite=TRUE)
# remove(LUX_TCD); remove(COL_TCD); remove(LUX_TCDpc); remove(COL_TCDpc)
# temp_remove<- list.files(pattern = "\\.*temp.*\\.tif$")
# sapply(temp_remove, unlink)

#_______________________________________________________________________________
# ##### GRASSLAND
# 
# ## Using data downloaded from Copernicus
# # Grassland Density; 2018; 10m; ; ESPG: 3035 (ETRS89, LAEA); Luxembourg & Germany; NA = 255
# # 0: No grassland; 1: Grassland
# 
# LUX_GRA_raster_file_folder <- "../ENV_DATA_EUROPE/GRA2018/GRA_2018_010m_lu_03035_v010/DATA"
# COL_GRA_raster_file_folder <- "../ENV_DATA_EUROPE/GRA2018/GRA_2018_010m_de_03035_v010/DATA"
# l_LUX <- list.files(LUX_GRA_raster_file_folder, pattern = "\\.tif$", full.names = TRUE)
# l_COL <- list.files(COL_GRA_raster_file_folder, pattern = "\\.tif$", full.names = TRUE)
# 
# template_LUX <- raster(extent(LUX_DEMpc))
# projection(template_LUX) <- crs_good
# template_COL <- raster(extent(COL_DEMpc))
# projection(template_COL) <- crs_good
# writeRaster(template_LUX, file="temp_GRA_LUX.tif", format="GTiff", overwrite=TRUE)
# writeRaster(template_COL, file="temp_GRA_COL.tif", format="GTiff", overwrite=TRUE)
# mosaic_rasters(gdalfile=l_LUX, dst_dataset="temp_GRA_LUX.tif", of="GTiff")
# mosaic_rasters(gdalfile=l_COL, dst_dataset="temp_GRA_COL.tif", of="GTiff")
# LUX_GRA <- raster("temp_GRA_LUX.tif")
# LUX_GRApc <- crop(LUX_GRA, e_LUX)
# NAvalue(LUX_GRApc) <- 255
# COL_GRA <- raster("temp_GRA_COL.tif")
# COL_GRApc <- crop(COL_GRA, e_COL)
# NAvalue(COL_GRApc) <- 255
# 
# # Checkpoint
# plot(LUX_GRApc, col=RColorBrewer::brewer.pal(9, "Greys"))
# points(SP_LUX@other$xy)
# plot(COL_GRApc, col=RColorBrewer::brewer.pal(9, "Greys"))
# points(SP_COL@other$xy)
# 
# # Clean-up and save the final products
# writeRaster(LUX_GRApc, file="Data/LUX_GRA.grd", overwrite=TRUE)
# writeRaster(COL_GRApc, file="Data/COL_GRA.grd", overwrite=TRUE)
# remove(LUX_GRA); remove(COL_GRA); remove(LUX_GRApc); remove(COL_GRApc)
# temp_remove<- list.files(pattern = "\\.*temp.*\\.tif$")
# sapply(temp_remove, unlink)

#_______________________________________________________________________________
# ##### WATER AND WETNESS
# 
# ## Using data downloaded from Copernicus
# # Occurrence of water and wet surfaces; 2018; 10m; ; ESPG: 3035 (ETRS89, LAEA); Luxembourg & Germany; NA = 255, 254 (unclassifiable), 253 (seawater)
# # 0: Dry; 1: Permanent water; 2: Temporary water; 3: Permanent wet; 4: Temporary wet
# 
# LUX_WAW_raster_file_folder <- "../ENV_DATA_EUROPE/WAW2018/WAW_2018_010m_lu_03035_v020/DATA"
# COL_WAW_raster_file_folder <- "../ENV_DATA_EUROPE/WAW2018/WAW_2018_010m_de_03035_v020/DATA"
# l_LUX <- list.files(LUX_WAW_raster_file_folder, pattern = "\\.tif$", full.names = TRUE)
# l_COL <- list.files(COL_WAW_raster_file_folder, pattern = "\\.tif$", full.names = TRUE)
# 
# template_LUX <- raster(extent(LUX_DEMpc))
# projection(template_LUX) <- crs_good
# template_COL <- raster(extent(COL_DEMpc))
# projection(template_COL) <- crs_good
# writeRaster(template_LUX, file="temp_WAW_LUX.tif", format="GTiff", overwrite=TRUE)
# writeRaster(template_COL, file="temp_WAW_COL.tif", format="GTiff", overwrite=TRUE)
# mosaic_rasters(gdalfile=l_LUX, dst_dataset="temp_WAW_LUX.tif", of="GTiff")
# mosaic_rasters(gdalfile=l_COL, dst_dataset="temp_WAW_COL.tif", of="GTiff")
# LUX_WAW <- raster("temp_WAW_LUX.tif")
# LUX_WAWpc <- crop(LUX_WAW, e_LUX)
# LUX_WAWpc <- reclassify(LUX_WAWpc, cbind(4, Inf, NA), right=TRUE)
# COL_WAW <- raster("temp_WAW_COL.tif")
# COL_WAWpc <- crop(COL_WAW, e_COL)
# COL_WAWpc <- reclassify(COL_WAWpc, cbind(4, Inf, NA), right=TRUE)
# 
# # Checkpoint
# plot(LUX_WAWpc, col=RColorBrewer::brewer.pal(5, "Blues")[c(1,5,4,3,2)])
# points(SP_LUX@other$xy)
# plot(COL_WAWpc, col=RColorBrewer::brewer.pal(5, "Blues")[c(1,5,4,3,2)])
# points(SP_COL@other$xy)
# 
# # Clean-up and save the final products
# writeRaster(LUX_WAWpc, file="Data/LUX_WAW.grd", overwrite=TRUE)
# writeRaster(COL_WAWpc, file="Data/COL_WAW.grd", overwrite=TRUE)
# remove(LUX_WAW); remove(COL_WAW); remove(LUX_WAWpc); remove(COL_WAWpc)
# temp_remove<- list.files(pattern = "\\.*temp.*\\.tif$")
# sapply(temp_remove, unlink)

# # If WAW is not enough, we can use this (not polygonized and rasterized yet)
# #_______________________________________________________________________________
# ##### RIVER NETWORK
# 
# list_rivers <- list.files("..//ENV_DATA_EUROPE/EU-HYDRO",
#                           pattern = "v013\\.gpkg$", recursive = TRUE, full.names = TRUE)
# list_rivers2 <- lapply(list_rivers, st_read)
# # Selecting perennial and on/above surface rivers
# # We need to use the select from dplyr
# list_rivers3 <- list_rivers2 %>% lapply(. %>% select(HYP,LOC) %>% filter(HYP==1 & LOC>40) %>% select(HYP) %>% st_zm)
# rivers <- st_union(list_rivers3[[1]], list_rivers3[[2]])
# rivers <- rivers %>% select(HYP)
# 
# st_write(rivers, "HOVERFLY_Rivers.gpkg", append=FALSE)
# 
# # We need to convert multi-line strings to polygons
# rivers_poly <- st_cast(rivers_poly, 'MULTIPOLYGON')
# 
# rivers_poly <- 
# 
# p_LUX <- as(e_LUX, 'SpatialPolygons')
# crs(p_LUX) <- crs_good
# LUX_RIV <- st_intersection(rivers, p_LUX)
# 
# rivers <- st_transform(rivers, crs_good)
# 
# # We need to rasterize in some way to be able to use it as a resistance surface.
# library(fasterize)
# fasterize(rivers_poly, LUX_DEMpc)
# # Let's start with 25m just like the DEM
# 
# # Simplification and reduction of overlap with river:
# # Reduce WAW to just permanent non-river water and wetness (permanent and temporary); a binary variable


