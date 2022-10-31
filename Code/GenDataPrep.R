################################################################################
######### Julian Wittische - Hoverfly landscape genetics #########
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

##### Syritta pipiens genind (adegenet)
source("Code/Syritta_fix.R")

# Read
SP_df_raw <- SP_df_raw_fix

# Set first column as row names and remove
SP_df <- SP_df_raw[,c(-1, -2, -3, -4, -5)]
SP_geo <-  SP_df_raw[,3:4]

# Fill fragment lengths with 0 to have three characters
SP_df <- apply(SP_df, 2, FUN=function(x){sprintf("%03d", x)})

# Create a column with the required format (combine columns by locus)
for (i in seq(1,ncol(SP_df),2)){
  SP_df[,i] <- as.character(paste(SP_df[,i], SP_df[,i+1], sep= "/"));
}

# Remove the now useless second column of each locus
SP_df <- SP_df[,-seq(2, ncol(SP_df), 2)]

# Restore locus names
colnames(SP_df) <- gsub( "\\..*$", "", colnames(SP_df))
row.names(SP_df) <- as.vector(SP_df_raw[,2])

# Transform into genind object
SP_genind <- df2genind(SP_df, ncode=3, ploidy=2, sep="/",
                       type="codom", NA.char=" NA/ NA")

# Add pop information (study area)
SP_genind@pop <- as.factor(substr(row.names(SP_df), 1, 2))

# Add coordinates and reproject them
SP_geo_sp <- SpatialPoints(SP_geo, CRS(SRS_string = "EPSG:4326"))
SP_geo_sp <-spTransform(SP_geo_sp, CRS(SRS_string = "EPSG:3035"))
SP_genind@other$xy <- SP_geo_sp

################################################################################
################################################################################
################################################################################

SP_genind@pop <- replace(SP_genind@pop, which(SP_genind@pop=="SW"), "LU")
SP_genind@pop  <-droplevels(SP_genind@pop)
table(SP_genind@pop)

SP_genind_LUX <- SP_genind[SP_genind@pop=="LU"]
SP_genind_COL <- SP_genind[SP_genind@pop=="CO"]

# Converst to run STRUCTURE
# Get rid of Spp141 (LD, HWE, high Fis)
SP <- SP_genind[loc=c("Spp010", "Spp053", "Spp080", "Spp142",
                                      "Spp231", "Spp273", "Spp476", "Spp051",
                                      "Spp108", "Spp313", "Spp360", "Spp391",
                                      "Spp416"), drop=TRUE]

genind2structure(SP,  file="Data/SP_STRU.txt", pops=TRUE)


##### Myathropa florea genind (adegenet)

# MF <- MF_genind
# 
# MF@pop <- replace(MF@pop, which(MF@pop=="SW"), "LU")
# MF@pop  <-droplevels(MF@pop)
# table(MF@pop)
# genind2structure(MF,  file="Data/MF_STRU.txt", pops=TRUE)

