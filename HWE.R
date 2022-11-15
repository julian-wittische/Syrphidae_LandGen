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

HW_SP_LU <- hw.test(SP_genind[SP_genind@pop=="LU"], B=10000)
HW_SP_CO <- hw.test(SP_genind[SP_genind@pop=="CO"], B=10000)

# Which loci are in HW disequilibrium after fdr correction?
names(which(p.adjust(HW_SP_LU[,4], "BH") < 0.05 & p.adjust(HW_SP_CO[,4], "BH") < 0.05))

# RESULT: Spp142, Spp051, Spp108, Spp141, Spp360  are in disequilibrium
# but let us keep them because Alain has subsampled 30 individuals 10 times
# and they pass (spring 2022)

#-------------------------------------------------------------------------------
########################### Stricter HWE testing ###############################
#-------------------------------------------------------------------------------

########## Removing duplicates

SP_genind_nodup <- SP_genind[!duplicated(as.data.frame(SP_genind@other$xy)),]

### Testing
HW_SP_LU <- hw.test(SP_genind_nodup[SP_genind_nodup@pop=="LU"], B=10000)
HW_SP_CO <- hw.test(SP_genind_nodup[SP_genind_nodup@pop=="CO"], B=10000)

### Which loci are in HW disequilibrium after fdr correction?
names(which(p.adjust(HW_SP_LU[,4], "BH") < 0.05 & p.adjust(HW_SP_CO[,4], "BH") < 0.05))

# RESULT: Spp142, Spp051, Spp108, Spp141, Spp360 - SAME AS WITH DUPLICATES

########## Hardy-Weinberg equilibrium within communes within study areas
lux_communes <- geoboundaries("Luxembourg", "adm2")
# adm2 is what we want (adm1 cantons)
crs(lux_communes)

# Transform into the same CRS used by the points
lux_communes <- st_transform(lux_communes, crs="EPSG:3035")
SP_genind_nodup@other$xy <- st_transform(SP_genind_nodup@other$xy, crs="EPSG:3035")

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
                              function(ls) hw.test(ls, B=1000)[,4]))
HWE.test.exact <- t(data.matrix(HWE.test))
{cat("Exact p-values (Monte Carlo procedure):", "\n")
  round(HWE.test.exact,3)}

### Subsample 20 in each 20+ commune



