library(raster)
library(sf)
library(sp)
library(latticeExtra)
library(dplyr)
library(rgeos)
library(rasterVis)
library(mapview)
library(stars)

############################### LUXEMBOURG CITY ###############################

# Load, define CRS, check
lux <- raster("ascii/lux_in.asc")
crs(lux) <- CRS('+init=EPSG:2169')


lux_stars <- read_stars("ascii/lux_in.asc")
