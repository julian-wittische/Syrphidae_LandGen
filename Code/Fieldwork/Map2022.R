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

table_swlux <- read_sheet("https://docs.google.com/spreadsheets/d/1hifgIF_Yw765mo0xqFZNHAfZ0U6C9fA46Bmeriz1Xe4/edit?usp=sharing", na="")
1

###### Important as fgb is not recognised by pandoc
mapviewOptions(fgb = FALSE)

sw <- raster("Data/ascii/esch_in.asc")
lux <- raster("Data/ascii/lux_in.asc")
crs(sw) <- CRS('+init=EPSG:2169')
crs(lux) <- CRS('+init=EPSG:2169')
swlux <- raster::merge(sw, lux)

swlux1km <- aggregate(swlux, fact=1000/res(swlux))
plot(swlux1km)

cell_number_swlux <- swlux1km
cell_number_swlux[!is.na(cell_number_swlux)] <- 1:length(cell_number_swlux[!is.na(cell_number_swlux)])
length(swlux1km[!is.na(swlux1km)])

rtp <- rasterToPolygons(cell_number_swlux, digits=20)

table_swlux <- as.data.frame(table_swlux)
#table_swlux <- table_swlux[-nrow(table_swlux),]

cells_done_swlux <- which(as.numeric(as.character(table_swlux$`A. cineraria`))>=1)

effort <- cell_number_swlux
effort[] <- NA
effort[which(values(cell_number_swlux)%in%cells_done_swlux)] <- cells_done_swlux

rtp_effort <- rasterToPolygons(effort, digits=20)

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

comb <- m+eff
comb_osm <- m_osm + eff_osm
cmap <- comb

cmap_osm <- comb_osm

mapshot(cmap, url="swluxmap2022_effort_sat.html")
mapshot(cmap_osm, url="swluxmap2022_effort_osm.html")

# ### Cheers
# sum(table_Col$`S. pipiens`,na.rm=TRUE)
# sum(as.numeric(unlist(table_Col$`M. florea` )), na.rm=TRUE)
# length(which(as.numeric(unlist(table_Col$`S. pipiens` ))>=1))
# length(which(as.numeric(as.character(table_Col$`M. florea`))>=1))
# #length(which(as.numeric(unlist(table_Col$`DONE?`))>=1))