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
lux <- raster("Data/ascii/lux_in.asc")

crs(lux) <- CRS('+init=EPSG:2169') # this is according to Alain
### Checking steps
# lux
# nlayers(lux)
# plot(lux)
# res(lux)
# sum(is.na(values(lux)))
# sum(!is.na(values(lux)))

### Change resolution to 1 kilometer, check
lux1km <- aggregate(lux, fact=1000/res(lux))
#### Checking steps
# res(lux1km) #checking step
# dim(lux1km)
# ncell(lux1km)
# plot(lux1km)

cell_number_lux <- lux1km
cell_number_lux[!is.na(cell_number_lux)] <- 1:length(cell_number_lux[!is.na(cell_number_lux)])

##### Create polygon from raster to be able to plot grid on interactive map
rtp <- rasterToPolygons(cell_number_lux, digits=20)

##### Read from our online document
table_Lux <- read_sheet("https://docs.google.com/spreadsheets/d/1spY_XlLYqT1DtNS_-Eqcz7Ezh2qW5cBISHSpY8tw3aI/edit?usp=sharing", na="")
table_Lux <- as.data.frame(table_Lux)

cells_done_Lux <- which(as.numeric(unlist(table_Lux$`S. pipiens`))>=1)

effort <- cell_number_lux
effort[] <- NA
effort[which(values(cell_number_lux)%in%cells_done_Lux)] <- cells_done_Lux

rtp_effort <- rasterToPolygons(effort, digits=20)

##### Plotting procedure
m_Lux <- mapview(rtp,
             method = "ngb", 
             na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
             query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
             trim = TRUE,
             legend = FALSE, #no need for legend
             map.types =  "Esri.WorldImagery",#,#"OpenStreetMap",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
             alpha.regions = 0,
             lwd=2,
             color="red") #get rid of color

eff_Lux <- mapview(rtp_effort,
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

comb_Lux <- m_Lux+eff_Lux

m_osm_Lux <- mapview(rtp,
                 method = "ngb", 
                 na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
                 query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
                 trim = TRUE,
                 legend = FALSE, #no need for legend
                 map.types =  "OpenStreetMap",#,#"OpenStreetMap",#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
                 alpha.regions = 0,
                 lwd=2,
                 color="red") #get rid of color

eff_osm_Lux <- mapview(rtp_effort,
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

comb_osm_Lux <- m_osm_Lux + eff_osm_Lux

#mapshot(comb, url="luxmap_effort_Volucella.html")

### Cheers

sum(table_Lux$`S. pipiens`,na.rm=TRUE)
sum(as.numeric(unlist(table_Lux$`M. florea` )),na.rm=TRUE)
length(which(as.numeric(unlist(table_Lux$`S. pipiens` ))>=1))
length(which(as.numeric(unlist(table_Lux$`M. florea` ))>=1))
################################################################################