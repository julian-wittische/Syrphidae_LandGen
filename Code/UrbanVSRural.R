LUX_IMP <- raster("Data/LUX_IMP.grd")
LUX_GRA <- raster("Data/LUX_GRA.grd")
COL_IMP <- raster("Data/COL_IMP.grd")
COL_GRA <- raster("Data/COL_GRA.grd")

plot(LUX_IMP)
plot(LUX_GRA)
plot(COL_IMP)
plot(COL_GRA)

par(mar=c(3,3,3,3))
hist(values(LUX_IMP))
hist(values(LUX_GRA))
hist(values(COL_IMP))
hist(values(COL_GRA))

plot(LUX_IMP>50)
plot(COL_IMP>50)
SP_genind_LUX@other$xy
which(extract(LUX_IMP, SP_genind_LUX@other$xy)>90)%in%which(extract(LUX_GRA, SP_genind_LUX@other$xy)>0)
which(extract(LUX_GRA, SP_genind_LUX@other$xy)>0)%in%which(extract(LUX_IMP, SP_genind_LUX@other$xy)>90)

SP_LUX_urban <- SP_genind_LUX[which(extract(LUX_IMP, SP_genind_LUX@other$xy)>90)]

SP_LUX_grass <- SP_genind_LUX[which(extract(LUX_GRA, SP_genind_LUX@other$xy)>0)]

wc(SP_LUX_urban)
wc(SP_LUX_grass)

SP_LUX_urban_summary <- summary(SP_LUX_urban)

mean(SP_LUX_urban_summary$Hexp)
sd(SP_LUX_urban_summary$Hexp)
mean(SP_LUX_urban_summary$Hobs)
sd(SP_LUX_urban_summary$Hobs)

SP_LUX_grass_summary <- summary(SP_LUX_grass)

mean(SP_LUX_grass_summary$Hexp)
sd(SP_LUX_grass_summary$Hexp)
mean(SP_LUX_grass_summary$Hobs)
sd(SP_LUX_grass_summary$Hobs)

