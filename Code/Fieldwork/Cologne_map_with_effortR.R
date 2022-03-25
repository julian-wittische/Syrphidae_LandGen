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
cologne <- raster("ascii/cologne_in.asc")

crs(cologne) <- CRS('+init=EPSG:2169') # this is according to Alain
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

cell_number_cologne <- cologne1km
cell_number_cologne[!is.na(cell_number_cologne)] <- 1:length(cell_number_cologne[!is.na(cell_number_cologne)])

##### Create polygon from raster to be able to plot grid on interactive map
rtp <- rasterToPolygons(cell_number_cologne, digits=20)

##### Read from our online document
table_Col <- read_sheet("https://docs.google.com/spreadsheets/d/1eey0sTASSw02a6Yvze8-UPziqRF1VGnqPX5RLcsWTmk/edit#gid=0", na="")
table_Col <- as.data.frame(table_Col)
table_Col <- table_Col[-nrow(table_Col),]

cells_done_Col <- which(as.numeric(as.character(table_Col$`S. pipiens`))>=1)
#cells_done_Col <- cells_done_Col[-length(cells_done_Col)]
effort <- cell_number_cologne
effort[] <- NA
effort[which(values(cell_number_cologne)%in%cells_done_Col)] <- cells_done_Col

rtp_effort <- rasterToPolygons(effort, digits=20)
#
cells_done_Col_MF <- which(as.numeric(as.character(table_Col$`M. florea` ))>=1)
#
effort_MF <- cell_number_cologne
effort_MF[] <- NA
effort_MF[which(values(cell_number_cologne)%in%cells_done_Col_MF)] <- cells_done_Col_MF

rtp_effort_MF <- rasterToPolygons(effort_MF, digits=20)
#
cells_done_Col_both <- which(as.numeric(as.character(table_Col$`M. florea`))>=1 & 
                               as.numeric(as.character(table_Col$`S. pipiens`))>=1)
#
effort_both <- cell_number_cologne
effort_both[] <- NA
effort_both[which(values(cell_number_cologne)%in%cells_done_Col_both)] <- cells_done_Col_both

rtp_effort_both <- rasterToPolygons(effort_both, digits=20)

##### Plotting procedure
m <- mapview(rtp,
             method = "ngb", 
             na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
             query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
             trim = TRUE,
             legend = FALSE, #no need for legend
             map.types =  "Esri.WorldImagery",#,#"OpenStreetMap",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
             alpha.regions = 0,
             lwd=2,
             color="red") #get rid of color

eff <- mapview(rtp_effort,
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

eff_MF <- mapview(rtp_effort_MF,
               method = "ngb", 
               na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
               query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
               trim = TRUE,
               legend = FALSE, #no need for legend
               map.types = "Esri.WorldImagery",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
               alpha.regions = 0.25,
               col.regions = "yellow",
               lwd=2,
               color="yellow") #get rid of color

eff_both <- mapview(rtp_effort_both,
                  method = "ngb", 
                  na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
                  query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
                  trim = TRUE,
                  legend = FALSE, #no need for legend
                  map.types = "Esri.WorldImagery",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
                  alpha.regions = 0.25,
                  col.regions = "green",
                  lwd=2,
                  color="green") #get rid of color

comb_Col <- m+eff

m_osm <- mapview(rtp,
             method = "ngb", 
             na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
             query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
             trim = TRUE,
             legend = FALSE, #no need for legend
             map.types =  "OpenStreetMap",#,#"OpenStreetMap",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
             alpha.regions = 0,
             lwd=2,
             color="red") #get rid of color

eff_osm <- mapview(rtp_effort,
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

eff_osm_MF <- mapview(rtp_effort_MF,
                   method = "ngb", 
                   na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
                   query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
                   trim = TRUE,
                   legend = FALSE, #no need for legend
                   map.types = "OpenStreetMap",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
                   alpha.regions = 0.25,
                   col.regions = "yellow",
                   lwd=2,
                   color="yellow") #get rid of color

eff_osm_both <- mapview(rtp_effort_both,
                      method = "ngb", 
                      na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
                      query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
                      trim = TRUE,
                      legend = FALSE, #no need for legend
                      map.types = "OpenStreetMap",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
                      alpha.regions = 0.25,
                      col.regions = "green",
                      lwd=2,
                      color="green") #get rid of color

comb_osm_Col <- m_osm + eff_osm
cmap <- comb_osm_Col+eff_MF+eff_both
mapshot(cmap, url="colognemap_effort_osm.html")
cmap
### Cheers
sum(table_Col$`S. pipiens`,na.rm=TRUE)
sum(as.numeric(unlist(table_Col$`M. florea` )), na.rm=TRUE)
length(which(as.numeric(unlist(table_Col$`S. pipiens` ))>=1))
length(which(as.numeric(as.character(table_Col$`M. florea`))>=1))
#length(which(as.numeric(unlist(table_Col$`DONE?`))>=1))
################################################################################