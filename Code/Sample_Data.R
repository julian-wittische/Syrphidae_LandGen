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

##### Syritta pipiens genind (adegenet)

# Read
SP_df_raw <- as.data.frame(readxl::read_excel("Data/SPipiens_raw.xlsx",
                                              sheet=1, .name_repair="minimal"))

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

SP_genind@loc.n.all
SP_genind@all.names
stepcheck <- data.frame(apply(SP_genind@tab, 2, as.numeric))
(lol <- stepcheck[,grep("Spp051", colnames(stepcheck))])
colSums(lol, na.rm=TRUE)
colSums(lol, na.rm=TRUE)[sort(names(colSums(lol, na.rm=TRUE)))]


SP_df_raw[which(stepcheck$Spp273.279>0),2]

stepcheck 


SP_genind@loc.n.all
colSums(stepcheck, na.rm=TRUE)[11:27]

################################################################################
################################################################################
################################################################################

# Get rid of Spp141 (LD, HWE, high Fis)
SP_genind_noSpp141 <- SP_genind[loc=c("Spp010", "Spp053", "Spp080", "Spp142",
                                      "Spp231", "Spp273", "Spp476", "Spp051",
                                      "Spp108", "Spp313", "Spp360", "Spp391",
                                      "Spp416"), drop=TRUE]

SP <- SP_genind_noSpp141
SP@pop <- replace(SP@pop, which(SP@pop=="SW"), "LU")
SP@pop  <-droplevels(SP@pop)
table(SP@pop)
SP




SP_LUX <- SP[SP@pop=="LU"]
SP_COL <- SP[SP@pop=="CO"]

# Converst to run STRUCTURE
genind2structure(SP,  file="Data/SP_STRU.txt", pops=TRUE)
# Convert to genepop to run Migraine
genind_to_genepop(SP_LUX, output = "Data/SP_LUX_no141.txt")

##### Myathropa florea genind (adegenet)

MF <- MF_genind

MF@pop <- replace(MF@pop, which(MF@pop=="SW"), "LU")
MF@pop  <-droplevels(MF@pop)
table(MF@pop)
genind2structure(MF,  file="Data/MF_STRU.txt", pops=TRUE)

#-------------------------------------------------------------------------------
################### Collecting and loading environmental data ##################
#-------------------------------------------------------------------------------

##### Reading our homemade rasters
# see "Env_Rasters.R" to see how those rasters were created

# /!\ RAW GRASSLAND (2 levels), WAW (5 levels) ARE CATEGORICAL /!\

#Luxembourg
LUX_DEM <- raster("Data/LUX_DEM.grd")
LUX_IMP <- raster("Data/LUX_IMP.grd")
LUX_TCD <- raster("Data/LUX_TCD.grd")
LUX_GRA <- raster("Data/LUX_GRA.grd")
LUX_WAW <- raster("Data/LUX_WAW.grd")
LUX_DEM <- resample(LUX_DEM, LUX_IMP, method='bilinear')

#Cologne
COL_DEM <- raster("Data/COL_DEM.grd")
COL_IMP <- raster("Data/COL_IMP.grd")
COL_TCD <- raster("Data/COL_TCD.grd")
COL_GRA <- raster("Data/COL_GRA.grd")
COL_WAW <- raster("Data/COL_WAW.grd")
COL_DEM <- resample(COL_DEM, COL_IMP, method='bilinear')

# If needed we could add river network (should be covered by WAW) and roads (should be covered by imperviousness)
# Small woody features high res raster has not been produced as of 27/07/2022