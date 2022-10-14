################################################################################
######### Julian Wittische - October 2021 - Raccoon landscape genetics #########
################################################################################

#-------------------------------------------------------------------------------
############################### Loading packages ###############################
#-------------------------------------------------------------------------------

##### Needed
library(adegenet)
library(sp)
library(raster)
library(graph4lg)
source("Code/genind2structure.R")

# # Loading setup script for JuliaCall
# JULIA_HOME <- "../../../AppData/Local/Programs/Julia-1.7.3/bin/"
# JuliaCall::julia_setup(JULIA_HOME)

#-------------------------------------------------------------------------------
################### Loading genetic and sample location data ###################
#-------------------------------------------------------------------------------

##### Myathropa florea genind (adegenet)

# Read
MF_df_raw <- as.data.frame(readxl::read_excel("Data/MFlorea_raw.xlsx",
                                              sheet=1, .name_repair="minimal"))

# Set first column as row names and remove
MF_df <- MF_df_raw[,c(-1, -2, -3, -4, -5)]
MF_geo <-  MF_df_raw[,3:4]

# Fill fragment lengths with 0 to have three characters
MF_df <- apply(MF_df, 2, FUN=function(x){sprintf("%03d", x)})

# Create a column with the required format (combine columns by locus)
for (i in seq(1,ncol(MF_df),2)){
  MF_df[,i] <- as.character(paste(MF_df[,i], MF_df[,i+1], sep= "/"));
}

# Remove the now useless second column of each locus
MF_df <- MF_df[,-seq(2, ncol(MF_df), 2)]

# Restore locus names
colnames(MF_df) <- gsub( "\\..*$", "", colnames(MF_df))
row.names(MF_df) <- as.vector(MF_df_raw[,2])

# Transform into genind object
MF_genind <- df2genind(MF_df, ncode=3, ploidy=2, sep="/",
                       type="codom", NA.char=" NA/ NA")

# Add pop information (study area)
MF_genind@pop <- as.factor(substr(row.names(MF_df), 1, 2))

# Add coordinates and reproject them
MF_geo_MF <- SpatialPoints(MF_geo, CRS(SRS_string = "EPSG:4326"))
MF_geo_MF <-spTransform(MF_geo_MF, CRS(SRS_string = "EPSG:3035"))
MF_genind@other$xy <- MF_geo_MF

MF_genind@loc.n.all
MF_genind@all.names
stepcheck <- data.frame(apply(MF_genind@tab, 2, as.numeric))
(lol <- stepcheck[,grep("MF269", colnames(stepcheck))])
colSums(lol, na.rm=TRUE)
colSums(lol, na.rm=TRUE)[sort(names(colSums(lol, na.rm=TRUE)))]


MF_df_raw[which(stepcheck$MF70.169>0),2]


