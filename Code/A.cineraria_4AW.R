################################################################################
########################### Julian Wittische - MNHNL ###########################
################################################################################

################################################################################
### Loading the INCOMPLETE list of A. cineraria females from various projects ##
################################################################################

# Libraries 
library(googledrive)
library(googlesheets4)
library(sf)
library(mapview)

# Important as fgb is not recognised by pandoc
mapviewOptions(fgb = FALSE)

# Loading from google drive document (Hinatea)
table_AC <- read_sheet("https://docs.google.com/spreadsheets/d/1_FW7Efzs22RwNe_9M9FvZcrQtUCAsw_G23qZQH5pl1U/edit?usp=sharing",
                       sheet = 4,
                       na="")
1

# Sorting, checking, and cleaning
table_AC <- table_AC[which(table_AC$`A. cineraria` > 0), ]
table_AC <- table_AC[, c("Longitude", "Latitude")]
str(table_AC) # columns should not be numeric
table_AC$Longitude # 107 is not latlong
table_AC <- table_AC[c(-107),] # get rid of it
table_AC <- as.data.frame(apply(table_AC, 2, as.numeric))
dim(table_AC)

# Get rid of samples far from Luxembourg
table_AC <- table_AC[which(table_AC$Latitude < 50.2),] 
dim(table_AC)

# Format change
AC_sf <- st_as_sf(table_AC[,1:2], coords = c("Longitude", "Latitude"), crs = 4326)

# Use a projection appropriate for Luxembourg
AC_sf <- st_transform(AC_sf, crs = 2169)

# Plotting
ACmap <- mapview(AC_sf, na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
                 query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
                 trim = TRUE,
                 legend = FALSE, #no need for legend
                 map.types = "OpenStreetMap",
                 alpha.regions = 0.5,
                 col.regions = "blue",
                 lwd = 1,
                 color="blue",
                 cex = 3.75)
ACmap

mapshot(ACmap, url="ACmap_INCOMPLETE_osm.html")


################################################################################
# NOTE: Not all of those are required, I will sort it out once
library(dplyr)
library(latticeExtra)

library(raster)
library(rasterVis)
library(RColorBrewer)
library(reshape2)
library(rgeos)
library(RgoogleMaps)

library(sp)
library(leaflet)
library(htmlwidgets)
library(ggplot2)
