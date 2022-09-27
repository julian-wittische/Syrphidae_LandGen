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

tab <- apply(as.matrix(table_swlux[,2:5]), 2, unlist)
sum(as.numeric(tab$`A. cineraria`), na.rm=T)


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

cells_done_swlux_BM <- which(as.numeric(as.character(table_swlux$`B. major`))>=1)
effort_BM <- cell_number_swlux
effort_BM[] <- NA
effort_BM[which(values(cell_number_swlux)%in%cells_done_swlux_BM)] <- cells_done_swlux_BM
rtp_effort_BM <- rasterToPolygons(effort_BM, digits=20)

cells_done_swlux_OB <- which(as.numeric(as.character(table_swlux$`O. bicornis`))>=1)
effort_OB <- cell_number_swlux
effort_OB[] <- NA
effort_OB[which(values(cell_number_swlux)%in%cells_done_swlux_OB)] <- cells_done_swlux_OB
rtp_effort_OB <- rasterToPolygons(effort_OB, digits=20)

cells_done_swlux_OC <- which(as.numeric(as.character(table_swlux$`O. cornuta`))>=1)
effort_OC <- cell_number_swlux
effort_OC[] <- NA
effort_OC[which(values(cell_number_swlux)%in%cells_done_swlux_OC)] <- cells_done_swlux_OC
rtp_effort_OC <- rasterToPolygons(effort_OC, digits=20)

cells_done_swlux_AC <- which(as.numeric(as.character(table_swlux$`A. cineraria`))>=1)
effort_AC <- cell_number_swlux
effort_AC[] <- NA
effort_AC[which(values(cell_number_swlux)%in%cells_done_swlux_AC)] <- cells_done_swlux_AC
rtp_effort_AC <- rasterToPolygons(effort_AC, digits=20)

satorosm <- "OpenStreetMap" #"Esri.WorldImagery" #  "OpenStreetMap" #

##### Plotting procedure
m <- mapview(rtp,
             method = "ngb", 
             na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
             query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
             trim = TRUE,
             legend = FALSE, #no need for legend
             map.types = satorosm ,#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
             alpha.regions = 0,
             lwd=2,
             color="red") #get rid of color

eff_BM <- mapview(rtp_effort_BM,
                  method = "ngb", 
                  na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
                  query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
                  trim = TRUE,
                  legend = FALSE, #no need for legend
                  map.types = satorosm,#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
                  alpha.regions = 0.25,
                  col.regions = "blue",
                  lwd=2,
                  color="blue") #get rid of color

eff_OB <- mapview(rtp_effort_OB,
               method = "ngb", 
               na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
               query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
               trim = TRUE,
               legend = FALSE, #no need for legend
               map.types = satorosm,#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
               alpha.regions = 0.25,
               col.regions = "black",
               lwd=2,
               color="black") #get rid of color

eff_OC <- mapview(rtp_effort_OC,
                  method = "ngb", 
                  na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
                  query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
                  trim = TRUE,
                  legend = FALSE, #no need for legend
                  map.types = satorosm,#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
                  alpha.regions = 0.25,
                  col.regions = "yellow",
                  lwd=2,
                  color="yellow") #get rid of color

eff_AC <- mapview(rtp_effort_AC,
                  method = "ngb", 
                  na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
                  query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
                  trim = TRUE,
                  legend = FALSE, #no need for legend
                  map.types = satorosm,#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
                  alpha.regions = 0.5,
                  col.regions = "red",
                  lwd=2,
                  color="red") #get rid of color

comb_BM <- m + eff_BM
comb_BM

comb_OB <- m + eff_OB
comb_OB

comb_OC <- m + eff_OC
comb_OC

comb_AC <- m + eff_AC
comb_AC

comb_all <- m + eff_BM + eff_OB + eff_OC + eff_AC
comb_all

mapshot(comb_BM, url="swluxmap2022BM_effort_osm.html")


# B. major
sum(unlist(as.numeric(as.character(table_swlux$`B. major`))), na.rm=TRUE)
length(which(unlist(as.numeric(as.character(table_swlux$`B. major`)))>=1))
length(which(unlist(as.numeric(as.character(table_swlux$`B. major`)))>=1))/457

# O. bicornis
sum(unlist(as.numeric(as.character(table_swlux$`O. bicornis`))), na.rm=TRUE)
length(which(unlist(as.numeric(as.character(table_swlux$`O. bicornis`)))>=1))
length(which(unlist(as.numeric(as.character(table_swlux$`O. bicornis`)))>=1))/457

# O. cornuta
sum(unlist(as.numeric(as.character(table_swlux$`O. cornuta`))), na.rm=TRUE)
length(which(unlist(as.numeric(as.character(table_swlux$`O. cornuta`)))>=1))
length(which(unlist(as.numeric(as.character(table_swlux$`O. cornuta`)))>=1))/457

# A. cineraria
sum(unlist(as.numeric(as.character(table_swlux$`A. cineraria`))), na.rm=TRUE)
length(which(unlist(as.numeric(as.character(table_swlux$`A. cineraria`)))>=1))
length(which(unlist(as.numeric(as.character(table_swlux$`A. cineraria`)))>=1))/457

# People
species <- table_swlux[, 2:5]
species  <- apply(species, 2, FUN=function(x){as.numeric(as.character(x))})
table(table_swlux[rowSums(species, na.rm=TRUE) != 0,"BOOKED BY"])

### Cheers

AC21 <- read.csv("Data/Acineraria2021.csv")
AC21 <- st_as_sf(AC21[,2:3], coords = c("Longitude", "Latitude"), crs = 4326)
AC21 <- st_transform(AC21, crs=2169)

AC21map <- mapview(AC21, na.color = rgb(0, 0, 255, max = 255, alpha = 0), #get rid of color
                      query.type = "click", #CLICK ON A PLACE TO KNOW WHICH CELL YOU ARE IN
                      trim = TRUE,
                      legend = FALSE, #no need for legend
                      map.types = satorosm,#"Esri.WorldImagery",#, # CHANGE TO "Esri.WorldImagery" IF YOU WANT
                      alpha.regions = 0.5,
                      col.regions = "blue",
                      lwd=2,
                      color="blue")

comb_AC + AC21map



                    