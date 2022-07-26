################################################################################
######### Julian Wittische - October 2021 - Raccoon landscape genetics #########
################################################################################

#-------------------------------------------------------------------------------
############################### Loading packages ###############################
#-------------------------------------------------------------------------------

# Needed
library(adegenet)
library(ResistanceGA)
library(sp)
library(raster)
# Need unchecked
library(sf)
library(ggplot2)
library(rgdal)

library(dplyr)

# Loading setup script for JuliaCall
JULIA_HOME <- "../../../AppData/Local/Programs/Julia-1.7.3/bin/"
JuliaCall::julia_setup(JULIA_HOME)

#-------------------------------------------------------------------------------
################### Loading genetic and sample location data ###################
#-------------------------------------------------------------------------------

##### Syritta pipiens genind (adegenet)

SP_genind_noSpp141_360_080 <- readRDS("Data/SP_genind_noSpp141_360_080.RDS")
SP <- SP_genind_noSpp141_360_080
SP_LUX <- SP[SP@pop=="SW"|SP@pop=="LU"]
SP_COL <- SP[SP@pop=="CO"]

# # We keep the crs from the various rasters from Copernicus (3035)
# crs_good <- CRS(SRS_string = "EPSG:3035")
# ger <- spTransform(raster::getData("GADM", country="DEU", level=1), CRSobj=crs_good)
# lux <- spTransform(raster::getData("GADM", country="LUX", level=1), CRSobj=crs_good)

#-------------------------------------------------------------------------------
################### Collecting and loading environmental data ##################
#-------------------------------------------------------------------------------
# see "Env_Rasters.R" to see how those rasters were created
# Reading our homemade rasters (code available to recreate below)

# /!\ RAW GRASSLAND (2 levels), WAW (5 levels) ARE CATEGORICAL /!\
DEM_ <- raster("RACCOON_DEM.grd")
IMP_ <- raster("RACCOON_IMP.grd")
TC_ <- raster("RACCOON_TC.grd")
GRA_ <- raster("RACCOON_GRA.grd")
WAW_ <- raster("RACCOON_WAW.grd")
# If needed we could add river network (should be covered by WAW) and roads (should be covered by imperviousness)
# Small woody features high res raster has not been produced as of 21/01/2022
