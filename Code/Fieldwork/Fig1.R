################################################################################
############################ Maps of hoverfly samples ##########################
################################################################################
# Project lead: Julian Wittische
# December 2022

#-------------------------------------------------------------------------------
############################### Loading packages ###############################
#-------------------------------------------------------------------------------

source("C:/Users/jwittische/Desktop/Projects/EnvironmentalData/CLC2018/LoadingCLC2018.R")
source("Code/GenDataPrep.R")

crs(SP_genind_LUX@other$xy) <- "EPSG:3035"
crs(SPproj@other$xy)

SPproj <- SP
SPproj@other$xy <- spTransform(SPproj@other$xy, CRS(SRS_string = "EPSG:3035"))
MFproj <- MF
MFproj@other$xy <- spTransform(MFproj@other$xy, CRS(SRS_string = "EPSG:2169"))

plot(cologne, col=rainbow(30))
points(SPproj@other$xy, col="black", pch=24, bg="white", cex=1)
points(MFproj@other$xy, col="black", pch=25, bg="black", cex=1)

north(xy=NULL, type=2, label="N", angle=0, head=0.1, xpd=TRUE)

cologne <- raster("Data/ascii/cologne_in.asc")
crs(cologne) <- CRS('+init=EPSG:2169')
SPproj <- SP
SPproj@other$xy <- spTransform(SPproj@other$xy, CRS(SRS_string = "EPSG:2169"))
MFproj <- MF
MFproj@other$xy <- spTransform(MFproj@other$xy, CRS(SRS_string = "EPSG:2169"))

plot(cologne, col=rainbow(30))
points(SPproj@other$xy, col="black", pch=24, bg="white", cex=1)
points(MFproj@other$xy, col="black", pch=25, bg="black", cex=1)

north(xy=NULL, type=2, label="N", angle=0, head=0.1, xpd=TRUE)

# # Legacy from older code, I am not sure what is really needed beyond sf
# library(dplyr)
# library(latticeExtra)
# library(mapview)
# library(RColorBrewer)
# library(reshape2)
# library(rgeos)
# library(RgoogleMaps)
# library(sp)
# library(leaflet)
# library(htmlwidgets)
# library(googledrive)
# library(googlesheets4)
# library(ggplot2)
# library(maps)
# library(terra)