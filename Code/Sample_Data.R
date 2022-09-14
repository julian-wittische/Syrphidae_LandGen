################################################################################
######### Julian Wittische - October 2021 - Raccoon landscape genetics #########
################################################################################

#-------------------------------------------------------------------------------
############################### Loading packages ###############################
#-------------------------------------------------------------------------------

##### Needed
library(adegenet)
library(ResistanceGA)
library(sp)
library(raster)
##### Need unchecked
# library(sf)
# library(ggplot2)
# library(rgdal)
# library(dplyr)

# Loading setup script for JuliaCall
JULIA_HOME <- "../../../AppData/Local/Programs/Julia-1.7.3/bin/"
JuliaCall::julia_setup(JULIA_HOME)

#-------------------------------------------------------------------------------
################### Loading genetic and sample location data ###################
#-------------------------------------------------------------------------------

##### Syritta pipiens genind (adegenet)

SP <- readRDS("Data/SP_genind_noSpp141.RDS")
SP_LUX <- SP[SP@pop=="SW"|SP@pop=="LU"]
SP_COL <- SP[SP@pop=="CO"]

SP@pop <- replace(SP@pop, which(SP@pop=="SW"), "LU")
SP@pop  <-droplevels(SP@pop)
table(SP@pop)
genind2structure(SP,  file="Data/SP_STRU.txt", pops=TRUE)

##### Myathropa florea genind (adegenet)

MF <- MF_genind

MF@pop <- replace(MF@pop, which(MF@pop=="SW"), "LU")
MF@pop  <-droplevels(MF@pop)
table(MF@pop)
genind2structure(MF,  file="Data/MF_STRU.txt", pops=TRUE)

#-------------------------------------------------------------------------------
################### Collecting and loading environmental data ##################
#-------------------------------------------------------------------------------

##### Reading our homemade rasters
# see "Env_Rasters.R" to see how those rasters were created

# /!\ RAW GRASSLAND (2 levels), WAW (5 levels) ARE CATEGORICAL /!\

#Luxembourg
LUX_DEM <- raster("Data/LUX_DEM.grd")
LUX_IMP <- raster("Data/LUX_IMP.grd")
LUX_TCD <- raster("Data/LUX_TCD.grd")
LUX_GRA <- raster("Data/LUX_GRA.grd")
LUX_WAW <- raster("Data/LUX_WAW.grd")
LUX_DEM <- resample(LUX_DEM, LUX_IMP, method='bilinear')

#Cologne
COL_DEM <- raster("Data/COL_DEM.grd")
COL_IMP <- raster("Data/COL_IMP.grd")
COL_TCD <- raster("Data/COL_TCD.grd")
COL_GRA <- raster("Data/COL_GRA.grd")
COL_WAW <- raster("Data/COL_WAW.grd")
COL_DEM <- resample(COL_DEM, COL_IMP, method='bilinear')

# If needed we could add river network (should be covered by WAW) and roads (should be covered by imperviousness)
# Small woody features high res raster has not been produced as of 27/07/2022