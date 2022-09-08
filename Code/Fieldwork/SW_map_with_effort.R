##### Julian Wittische 
##### April 2021


# /!\ I included checking steps for your own curiosity, no need to run them /!\


### Please update your packages regularly, I am just auto-installing missing packages here
list.of.packages <- c( "RColorBrewer", "raster", "sf", "dplyr",
                       "mapview", "reshape2", "RgoogleMaps",
                       "rasterVis", "latticeExtra", "rgeos", "sp", "ggplot2"
                       , "googlesheets4", "googledrive")

### This will install the packages you have not installed yet
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
install.packages(new.packages)

### /!\ you MIGHT need to install gdal or Rtools if your computer does not have them /!\

###### Installing and loading packages
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

###### Important as fgb is not recognised by pandoc
mapviewOptions(fgb = FALSE)

################################# Authentication #might be weird the first time, just run everyhting several times

sheets_auth(scope = "https://www.googleapis.com/auth/drive")
1
drive_auth(email="jwittische@gmail.com")
drive_auth(token = sheets_token())
1

##### Load the ascii file from Alain, define CRS, resample
# to simplify it for you I put this interactive line, pick the .asc file
SW <- raster("Data/ascii/esch_in.asc")

crs(SW) <- CRS('+init=EPSG:2169') # this is according to Alain
### Checking steps
# SW
# nlayers(SW)
# plot(SW)
# res(SW)
# sum(is.na(values(SW)))
# sum(!is.na(values(SW)))

### Change resolution to 1 kilometer, check
SW1km <- aggregate(SW, fact=1000/res(SW))
#### Checking steps
# res(SW1km) #checking step
# dim(SW1km)
# ncell(SW1km)
# plot(SW1km)

cell_number_SW <- SW1km
cell_number_SW[!is.na(cell_number_SW)] <- 1:length(cell_number_SW[!is.na(cell_number_SW)])

##### Create polygon from raster to be able to plot grid on interactive map
rtp <- rasterToPolygons(cell_number_SW, digits=20)

##### Read from our online document
table_SW <- read_sheet("https://docs.google.com/spreadsheets/d/16uQnmVqKACVurcTPiauHRVpJj-C85lDpAcUeNXWjjxo/edit#gid=0", na="")
table_SW <- as.data.frame(table_SW)

cells_done_SW <- which(as.numeric(unlist(table_SW$`S. pipiens` ))>=1)

effort <- cell_number_SW
effort[] <- NA
effort[which(values(cell_number_SW)%in%cells_done_SW)] <- cells_done_SW

rtp_effort <- rasterToPolygons(effort, digits=20)

##### Plotting procedure
m_SW <- mapview(rtp,
             method = "ngb", 
             na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
             query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
             trim = TRUE,
             legend = FALSE, #no need for legend
             map.types =  "Esri.WorldImagery",#,#"OpenStreetMap",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
             alpha.regions = 0,
             lwd=2,
             color="red") #get rid of color

eff_SW <- mapview(rtp_effort,
               method = "ngb", 
               na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
               query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
               trim = TRUE,
               legend = FALSE, #no need for legend
               map.types = "Esri.WorldImagery",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
               alpha.regions = 0.25,
               col.regions = "blue",
               lwd=2,
               color="blue") #get rid of color

comb_SW <- m_SW+eff_SW

m_osm_SW <- mapview(rtp,
                 method = "ngb", 
                 na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
                 query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
                 trim = TRUE,
                 legend = FALSE, #no need for legend
                 map.types =  "OpenStreetMap",#,#"OpenStreetMap",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
                 alpha.regions = 0,
                 lwd=2,
                 color="red") #get rid of color

eff_osm_SW <- mapview(rtp_effort,
                   method = "ngb", 
                   na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
                   query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
                   trim = TRUE,
                   legend = FALSE, #no need for legend
                   map.types = "OpenStreetMap",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
                   alpha.regions = 0.25,
                   col.regions = "blue",
                   lwd=2,
                   color="blue") #get rid of color

comb_osm_SW <- m_osm_SW + eff_osm_SW

#mapshot(comb, url="SWmap_effort_Syritta_satellite.html")

### Cheers

sum(table_SW$`S. pipiens`,na.rm=TRUE)
length(which(as.numeric(unlist(table_SW$`S. pipiens` ))>=1))
length(which(as.numeric(unlist(table_SW$`M. florea` ))>=1))
################################################################################