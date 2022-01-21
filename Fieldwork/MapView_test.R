library(raster)
library(sf)
library(sp)
library(latticeExtra)
library(dplyr)
library(rgeos)
library(rasterVis)
library(reshape2)
library(RColorBrewer)
library(RgoogleMaps)
library(mapview)

mapview(cell_number_cologne,
        method = "ngb", 
        na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
        query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
        trim = TRUE,
        legend = FALSE, #no need for legend
        map.types = "Esri.WorldImagery", #get satellite
        alpha.regions = 0.25) #get rid of color

mapview(rtp,
        method = "ngb", 
        na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
        query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
        trim = TRUE,
        legend = FALSE, #no need for legend
        map.types = "Esri.WorldImagery", #get satellite
        alpha.regions = 0) #get rid of color

rtp <- rasterToPolygons(cologne1km, digits=20)

mapview(rtp)

map.types

grid <- SpatialPixels(SpatialPoints(coordinates(cologne1km)[!is.na(values(cologne1km)),]))
crs(grid) <- CRS('+init=EPSG:2169')
plot(grid)


mapview(grid, method = "bilinear",
        na.color = rgb(0, 0, 255, max = 255, alpha = 0),
        query.type = "click",
        trim = TRUE,
        grid=TRUE,
        legend = FALSE)



g <- as(cologne1km, 'SpatialGridDataFrame')
crs(g) <- CRS('+init=EPSG:2169')

mapview(g)

g <- spTransform(g, CRS("+proj=longlat +datum=WGS84"))

# teams_latlong <- projectRaster(teams, crs = CRS("+proj=longlat +datum=WGS84"))
# grid_latlong <- as(teams_latlong, 'SpatialGridDataFrame')
# 
# mapview(grid_latlong, method = "ngb",
#         na.color = rgb(0, 0, 255, max = 255, alpha = 0),
#         query.type = "click",
#         trim = TRUE, 
#         legend = FALSE, grid = TRUE)
