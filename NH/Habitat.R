######################### Julian Wittische - MNHNL
######################### January 2022

# This subproject aims to investigate what explained netting success and by
# extension on abundance for S. pipiens and M. florea

# Libraries
library("raster")
library("sf")

##### Load the ascii file from Alain, define CRS, resample
# to simplify it for you I put this interactive line, pick the .asc file
cologne <- raster("ascii/cologne_in.asc")
#Issue somewhere so I am doing it myself

# Fieldwork metadata
mdcol <- read.csv("NH/Fieldwork_metadata - SW.csv")

# From Julian or Hinatea
mdcol <- mdcol[mdcol$Collector.name%in%c("Julian Wittische", "Hinatea Ariey",
                                         "HInatea Ariey"),]

# Get only relevant columns
mdcolrel <- mdcol[,c("Date", "Latitude", "Longitude", "Sampling.station",
                     "M..florea", "S..pipiens")]
SP <-mdcolrel[mdcol$S..pipiens>0,c(1:4,6)]
MF <-mdcolrel[mdcol$M..florea>0,1:5]
SP[,2:4] <- apply(SP[,2:4],2,as.numeric)
MF[,2:4] <- apply(MF[,2:4],2,as.numeric)

# Project sampling coordinates to the right CRS
crs(cologne) <- CRS('+init=EPSG:2169') # this is according to Alain
SP[,2:3] <- as.data.frame(spTransform(SpatialPoints(SP[,3:2], CRS(SRS_string = "EPSG:4326")),
                     CRS(SRS_string = "EPSG:3035")))
MF[,2:3] <- as.data.frame(spTransform(SpatialPoints(MF[,3:2], CRS(SRS_string = "EPSG:4326")),
                        CRS(SRS_string = "EPSG:3035")))

plot(cologne)
points(SP[,2:3])
