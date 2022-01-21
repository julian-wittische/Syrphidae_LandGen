######################### Julian Wittische - MNHNL
######################### January 2022

# This subproject aims to investigate what explained netting success and by
# extension on abundance for S. pipiens and M. florea

# Libraries
library("raster")

##### Load the ascii file from Alain, define CRS, resample
# to simplify it for you I put this interactive line, pick the .asc file
cologne <- raster("ascii/cologne_in.asc")


# Fieldwork metadata
mdcol <- read.csv("Fieldwork_metadata - SW.csv")

# From Julian or Hinatea
mdcol <- mdcol[mdcol$Collector.name%in%c("Julian Wittische", "Hinatea Ariey", "HInatea Ariey"),]

# Get only relevant columns
mdcolrel <- mdcol[,c("Date", "Latitude", "Longitude")]
SP <-mdcolrel[mdcol$S..pipiens>0,]
MF <-mdcolrel[mdcol$M..florea>0,]

