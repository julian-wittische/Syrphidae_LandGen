##### Julian Wittische 
##### April 2021

# /!\ I included checking steps for your own curiosity, no need to run them /!\

###### Installing and loading packages

### Please update your packages regularly, I am just auto-installing missing packages here
list.of.packages <- c( "RColorBrewer", "raster", "sf", "dplyr",
                       "mapview", "reshape2", "RgoogleMaps",
                      "rasterVis", "latticeExtra", "rgeos", "sp")

### This will install the packages you have not installed yet
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
install.packages(new.packages)

### /!\ you MIGHT need to install gdal or Rtools if your computer does not have them /!\

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

##### Load the ascii file from Alain, define CRS, resample
# to simplify it for you I put this interactive line, pick the .asc file
cologne <- raster(file.choose())

crs(cologne) <- CRS('+init=EPSG:2169') # this is accorfing to Alain
### Checking steps
# cologne
# nlayers(cologne)
# plot(cologne)
# res(cologne)
# sum(is.na(values(cologne)))
# sum(!is.na(values(cologne)))

### Change resolution to 1 kilometer, check
cologne1km <- aggregate(cologne, fact=1000/res(cologne))
#### Checking steps
# res(cologne1km) #checking step
# dim(cologne1km)
# ncell(cologne1km)
# plot(cologne1km)

##### Create polygon from raster to be able to plot grid on interactive map
rtp <- rasterToPolygons(cologne1km, digits=20)


##### Plotting procedure
mapview(rtp,
        method = "ngb", 
        na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
        query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
        trim = TRUE,
        legend = FALSE, #no need for legend
        map.types = "Esri.WorldImagery", # CHANGE TO "OpenStreetMap" IF YOU WANT
        alpha.regions = 0) #get rid of color

### Cheers