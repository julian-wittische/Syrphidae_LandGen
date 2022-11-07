library(dplyr)
library(latticeExtra)
library(mapview)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(reshape2)
library(rgeos)
library(RgoogleMaps)
library(sf)
library(sp)
library(leaflet)
library(htmlwidgets)
library(googledrive)
library(googlesheets4)
library(ggplot2)

source("Code/GenDataPrep.R")

cologne <- raster("Data/ascii/cologne_in.asc")
crs(cologne) <- CRS('+init=EPSG:2169')
SPproj <- SP
SPproj@other$xy <- spTransform(SPproj@other$xy, CRS(SRS_string = "EPSG:2169"))
MFproj <- MF
MFproj@other$xy <- spTransform(MFproj@other$xy, CRS(SRS_string = "EPSG:2169"))

plot(cologne, col=rainbow(30))
points(SPproj@other$xy, col="white", pch=19, cex=1.25)
points(MFproj@other$xy, col="black", pch=15, cex=1.25)

