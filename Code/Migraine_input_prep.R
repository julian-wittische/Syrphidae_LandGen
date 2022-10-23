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

# Get rid of Spp141 (LD, HWE, high Fis)
SP_genind_noSpp141 <- SP_genind[loc=c("Spp010", "Spp053", "Spp080", "Spp142",
                                      "Spp231", "Spp273", "Spp476", "Spp051",
                                      "Spp108", "Spp313", "Spp360", "Spp391",
                                      "Spp416"), drop=TRUE]

SP <- SP_genind_noSpp141
SP@pop <- replace(SP@pop, which(SP@pop=="SW"), "LU")
SP@pop  <-droplevels(SP@pop)

# Convert to run STRUCTURE
genind2structure(SP,  file="Data/SP_STRU.txt", pops=TRUE)

# SP_LUX <- SP[SP@pop=="LU"]
# SP_COL <- SP[SP@pop=="CO"]
# 
# 
# # Convert to genepop to run Migraine
# genind_to_genepop(SP_LUX, output = "Data/SP_LUX_no141.txt")

##### Myathropa florea genind (adegenet)
source("Code/Myathropa_fix.R")

# Read
MF_df_raw <- MF_df_raw_fix

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
MF_geo_sp <- SpatialPoints(MF_geo, CRS(SRS_string = "EPSG:4326"))
MF_geo_sp <-spTransform(MF_geo_sp, CRS(SRS_string = "EPSG:3035"))
MF_genind@other$xy <- MF_geo_sp

MF <- MF_genind
MF@pop <- replace(MF@pop, which(MF@pop=="SW"), "LU")
MF@pop  <-droplevels(MF@pop)
table(MF@pop)

MF_LUX <- MF[MF@pop=="LU"]
MF_COL <- MF[MF@pop=="CO"]

# Converst to run STRUCTURE
source("Code/genind2structure.R")
genind2structure(MF,  file="Data/MF_STRU.txt", pops=TRUE)
# Convert to genepop to run Migraine
genind_to_genepop(MF_LUX, output = "Data/MF_LUX_no141.txt")

