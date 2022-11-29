################################################################################
######### Julian Wittische - Hoverfly landscape genetics #########
################################################################################

#-------------------------------------------------------------------------------
############################### Loading packages ###############################
#-------------------------------------------------------------------------------

# Source dta loading script (SP_genind is the object with all 14 loci)

source("Code/GenDataPrep.R")

#-------------------------------------------------------------------------------
############### Hardy-Weinberg equilibrium within study areas ##################
#-------------------------------------------------------------------------------

# HW_SP_LU <- hw.test(SP_genind[SP_genind@pop=="LU"], B=1000)
# HW_SP_CO <- hw.test(SP_genind[SP_genind@pop=="CO"], B=1000)
# 
# # Which loci are in HW disequilibrium after fdr correction?
# names(which(p.adjust(HW_SP_LU[,4], "BH") < 0.05 & p.adjust(HW_SP_CO[,4], "BH") < 0.05))

# RESULT: Spp142, Spp051, Spp108, Spp141, Spp360  are in disequilibrium
# but let us keep them because Alain has subsampled 30 individuals 10 times
# and they pass (spring 2022)
# 
# HW_LU <- hw.test(MF_genind[MF_genind@pop=="LU"], B=1000)
# HW_CO <- hw.test(MF_genind[MF_genind@pop=="CO"], B=1000)
# 
# names(which(HW_LU[,4] < 0.05 & HW_CO[,4] < 0.05))
# RESULT: MFp142, MFp051, MFp108, MFp141, MFp360  are somewhat in disequilibrium
# but let us keep them because Alain has subsampled 30 individuals 10 times and they pass


#-------------------------------------------------------------------------------
########################### Stricter HWE testing ###############################
#-------------------------------------------------------------------------------

########## Removing duplicates

SP_genind_nodup <- SP_genind[!duplicated(as.data.frame(SP_genind@other$xy)),]

# ### Testing
# HW_SP_LU <- hw.test(SP_genind_nodup[SP_genind_nodup@pop=="LU"], B=1000)
# HW_SP_CO <- hw.test(SP_genind_nodup[SP_genind_nodup@pop=="CO"], B=1000)
# 
# ### Which loci are in HW disequilibrium after fdr correction?
# names(which(p.adjust(HW_SP_LU[,4], "BH") < 0.05 & p.adjust(HW_SP_CO[,4], "BH") < 0.05))

# RESULT: Spp142, Spp051, Spp108, Spp141, Spp360 - SAME AS WITH DUPLICATES

MF_genind_nodup <- MF_genind[!duplicated(as.data.frame(MF_genind@other$xy)),]

# ### Testing
# HW_MF_LU <- hw.test(MF_genind_nodup[MF_genind_nodup@pop=="LU"], B=100000)
# HW_MF_CO <- hw.test(MF_genind_nodup[MF_genind_nodup@pop=="CO"], B=100000)
# 
# ### Which loci are in HW disequilibrium after fdr correction?
# names(which(p.adjust(HW_MF_LU[,4], "BH") < 0.05 & p.adjust(HW_MF_CO[,4], "BH") < 0.05))

# RESULT:  - SAME AS WITH DUPLICATES

#-------------------------------------------------------------------------------
#### SP Hardy-Weinberg equilibrium within Luxembourg communes
#-------------------------------------------------------------------------------

lux_communes <- geoboundaries("Luxembourg", "adm2")
# adm2 is what we want (adm1 cantons)
crs(lux_communes)

# Transform into the same CRS used by the points
lux_communes <- st_transform(lux_communes, crs="EPSG:3035")
SP_genind_nodup@other$xy <- st_transform(st_as_sf(SP_genind_nodup@other$xy), crs="EPSG:3035")

# Plot to check samples distribution
ggplot(data = lux_communes) +
  geom_sf() +
  geom_sf(data=SP_genind_nodup[SP_genind_nodup@pop=="LU"]@other$xy,
          aes(col="red")) + theme(legend.position = "none")

# Looks like a decent subdivision for the subsampling technique suggested by AF
SP_nodup_communes <-st_join(SP_genind_nodup@other$xy, lux_communes, join = st_within)

# Count the number of individuals per commune
enough <- as.data.frame(count(as_tibble(SP_nodup_communes), shapeName))

# Select only communes with more than twenty individuals 
enough <- enough[enough$n>20,]

# Select only genotypes from individuals from those 20+ communes
SP_nodup_compop <- SP_genind_nodup
SP_nodup_compop@pop <- as.factor(SP_nodup_communes$shapeName)
#  Last one is NA
SP_nodup_compop <- SP_nodup_compop[SP_nodup_compop@pop %in% enough$shapeName[-11]]

# Testing
HWE.test <- data.frame(sapply(seppop(SP_nodup_compop),
                              function(ls) hw.test(ls, B=1000000)[,4]))
HWE.test.exact <- t(data.matrix(HWE.test))
{cat("Exact p-values (Monte Carlo procedure):", "\n")
  round(HWE.test.exact,5)}

write.csv(HWE.test.exact, "HWE_SP.csv")
t(apply(HWE.test.exact, 2, FUN=function(x) p.adjust(x, "BH")))


# ### Subsample 20 in each 20+ commune
# 
# HWE.test <- data.frame(sapply(seppop(SP_nodup_compop), 
#                               function(ls) hw.test(ls[sample(1:nrow(ls@tab),20),], B=100000)[,4]))
# HWE.test.exact <- t(data.matrix(HWE.test))
# HWE.test.exact <- apply(HWE.test.exact, 2, FUN=function(x) p.adjust(x, "BH"))
# {cat("Exact adjusted p-values (Monte Carlo procedure):", "\n")
#   round(HWE.test.exact,3)}

#-------------------------------------------------------------------------------
### SP Hardy-Weinberg equilibrium within Cologne districts
#-------------------------------------------------------------------------------

cologne_districts <- sf::st_read(dsn = "Data/Stadtteile_-_K%C3%B6ln/Stadtteile_Koeln.shp")
cologne_districts_agg <- aggregate(cologne_districts, FUN=min, by=list(cologne_districts$NR_STADTBE))

ggplot(data = cologne_districts_agg) +
  geom_sf()

cologne_districts <- st_transform(cologne_districts, crs="EPSG:3035")
cologne_districts_agg <- st_transform(cologne_districts_agg, crs="EPSG:3035")

SP_nodup_cologne <-st_join(SP_genind_nodup@other$xy, cologne_districts, join = st_within)

enough_COL <- as.data.frame(count(as_tibble(SP_nodup_cologne), STADTBEZIR))
enough_COL <- enough_COL[enough_COL$n>20,]

# Select only genotypes from individuals from those 20+ districts
SP_nodup_compop_COL <- SP_genind_nodup
SP_nodup_compop_COL@pop <- as.factor(SP_nodup_cologne$STADTBEZIR)
#  Last one is NA
SP_nodup_compop_COL <- SP_nodup_compop_COL[SP_nodup_compop_COL@pop %in% enough_COL$STADTBEZIR[-10]]

# Testing
HWE.test_SP_COL <- data.frame(sapply(seppop(SP_nodup_compop_COL),
                              function(ls) hw.test(ls, B=1000000)[,4]))
HWE.test.exact_SP_COL <- t(data.matrix(HWE.test_SP_COL))
{cat("Exact p-values (Monte Carlo procedure):", "\n")
  round(HWE.test.exact_SP_COL,5)}


write.csv(HWE.test.exact_SP_COL, "HWE_SP_COL.csv")
t(apply(HWE.test.exact_SP_COL, 2, FUN=function(x) p.adjust(x, "BH")))# ### Subsample 20 in each 20+ district
# 
# HWE.test_COL <- data.frame(sapply(seppop(SP_nodup_compop_COL), 
#                               function(ls) hw.test(ls[sample(1:nrow(ls@tab),20),], B=100000)[,4]))
# HWE.test.exact_COL <- t(data.matrix(HWE.test_COL))
# HWE.test.exact_COL <- apply(HWE.test.exact_COL, 2, FUN=function(x) p.adjust(x, "BH"))
# {cat("Exact adjusted p-values (Monte Carlo procedure):", "\n")
#   round(HWE.test.exact_COL,3)}

#-------------------------------------------------------------------------------
#### MF Hardy-Weinberg equilibrium within Luxembourg communes
#-------------------------------------------------------------------------------

MF_genind_nodup@other$xy <- st_transform(st_as_sf(MF_genind_nodup@other$xy), crs="EPSG:3035")

# Plot to check samples distribution
ggplot(data = lux_communes) +
  geom_sf() +
  geom_sf(data=MF_genind_nodup[MF_genind_nodup@pop=="LU"]@other$xy,
          aes(col="red")) + theme(legend.position = "none")

# Looks like a decent subdivision for the subsampling technique suggested by AF
MF_nodup_communes <-st_join(MF_genind_nodup@other$xy, lux_communes, join = st_within)

# Count the number of individuals per commune
enough <- as.data.frame(count(as_tibble(MF_nodup_communes), shapeName))

# Select only communes with more than twenty individuals 
enough <- enough[enough$n>11,]

# Select only genotypes from individuals from those 20+ communes
MF_nodup_compop <- MF_genind_nodup
MF_nodup_compop@pop <- as.factor(MF_nodup_communes$shapeName)
#  Last one is NA
MF_nodup_compop <- MF_nodup_compop[MF_nodup_compop@pop %in% enough$shapeName[-8]]

# Testing
HWE.test <- data.frame(sapply(seppop(MF_nodup_compop),
                              function(ls) hw.test(ls, B=1000000)[,4]))
HWE.test.exact <- t(data.matrix(HWE.test))
{cat("Exact p-values (Monte Carlo procedure):", "\n")
  round(HWE.test.exact,5)}

write.csv(HWE.test.exact, "HWE_MF.csv")
t(apply(HWE.test.exact, 2, FUN=function(x) p.adjust(x, "BH")))


# ### Subsample 20 in each 20+ commune
# 
# HWE.test <- data.frame(sapply(seppop(MF_nodup_compop), 
#                               function(ls) hw.test(ls[sample(1:nrow(ls@tab),20),], B=100000)[,4]))
# HWE.test.exact <- t(data.matrix(HWE.test))
# HWE.test.exact <- apply(HWE.test.exact, 2, FUN=function(x) p.adjust(x, "BH"))
# {cat("Exact adjusted p-values (Monte Carlo procedure):", "\n")
#   round(HWE.test.exact,3)}

#-------------------------------------------------------------------------------
### MF Hardy-Weinberg equilibrium within Cologne districts
#-------------------------------------------------------------------------------

MF_nodup_cologne <-st_join(MF_genind_nodup@other$xy, cologne_districts, join = st_within)

enough_COL <- as.data.frame(count(as_tibble(MF_nodup_cologne), STADTBEZIR))
enough_COL <- enough_COL[enough_COL$n>11,]

# Select only genotypes from individuals from those 20+ districts
MF_nodup_compop_COL <- MF_genind_nodup
MF_nodup_compop_COL@pop <- as.factor(MF_nodup_cologne$STADTBEZIR)
#  Last one is NA
MF_nodup_compop_COL <- MF_nodup_compop_COL[MF_nodup_compop_COL@pop %in% enough_COL$STADTBEZIR[-10]]

# Testing
HWE.test_MF_COL <- data.frame(sapply(seppop(MF_nodup_compop_COL),
                                     function(ls) hw.test(ls, B=1000000)[,4]))
HWE.test.exact_MF_COL <- t(data.matrix(HWE.test_MF_COL))
{cat("Exact p-values (Monte Carlo procedure):", "\n")
  round(HWE.test.exact_MF_COL,5)}


write.csv(HWE.test.exact_MF_COL, "HWE_MF_COL.csv")
t(apply(HWE.test.exact_MF_COL, 2, FUN=function(x) p.adjust(x, "BH")))# ### Subsample 20 in each 20+ district
# 
# HWE.test_COL <- data.frame(sapply(seppop(MF_nodup_compop_COL), 
#                               function(ls) hw.test(ls[sample(1:nrow(ls@tab),20),], B=100000)[,4]))
# HWE.test.exact_COL <- t(data.matrix(HWE.test_COL))
# HWE.test.exact_COL <- apply(HWE.test.exact_COL, 2, FUN=function(x) p.adjust(x, "BH"))
# {cat("Exact adjusted p-values (Monte Carlo procedure):", "\n")
#   round(HWE.test.exact_COL,3)}
                            