source("Code/GenDataPrep.R")
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

################################################################################

SP_LUX_urban <- SP_genind_LUX[which(extract(LUX_IMP, SP_genind_LUX@other$xy)>90)]

SP_LUX_grass <- SP_genind_LUX[which(extract(LUX_GRA, SP_genind_LUX@other$xy)>0)]

SP_COL_urban <- SP_genind_COL[which(extract(COL_IMP, SP_genind_COL@other$xy)>90)]

SP_COL_grass <- SP_genind_COL[which(extract(COL_GRA, SP_genind_COL@other$xy)>0)]

wc(SP_LUX_urban)
wc(SP_LUX_grass)

SP_LUX_urban_summary <- summary(SP_LUX_urban)

mean(SP_LUX_urban_summary$Hexp)
sd(SP_LUX_urban_summary$Hexp)
mean(SP_LUX_urban_summary$Hobs)
sd(SP_LUX_urban_summary$Hobs)

mean(allelicrichness(as.loci(SP_LUX_urban)))
sd(allelicrichness(as.loci(SP_LUX_urban)))

SP_LUX_grass_summary <- summary(SP_LUX_grass)

mean(SP_LUX_grass_summary$Hexp)
sd(SP_LUX_grass_summary$Hexp)
mean(SP_LUX_grass_summary$Hobs)
sd(SP_LUX_grass_summary$Hobs)

SP_COL_urban_summary <- summary(SP_COL_urban)

mean(SP_COL_urban_summary$Hexp)
sd(SP_COL_urban_summary$Hexp)
mean(SP_COL_urban_summary$Hobs)
sd(SP_COL_urban_summary$Hobs)

mean(allelicrichness(as.loci(SP_COL_urban)))
sd(allelicrichness(as.loci(SP_COL_urban)))

SP_COL_grass_summary <- summary(SP_COL_grass)

mean(SP_COL_grass_summary$Hexp)
sd(SP_COL_grass_summary$Hexp)
mean(SP_COL_grass_summary$Hobs)
sd(SP_COL_grass_summary$Hobs)

mean(allelicrichness(as.loci(SP_COL_grass)))
sd(allelicrichness(as.loci(SP_COL_grass)))

################################################################################

MF_LUX_urban <- MF_genind_LUX[which(extract(LUX_IMP, MF_genind_LUX@other$xy)>90)]

MF_LUX_grass <- MF_genind_LUX[which(extract(LUX_GRA, MF_genind_LUX@other$xy)>0)]

MF_COL_urban <- MF_genind_COL[which(extract(COL_IMP, MF_genind_COL@other$xy)>90)]

MF_COL_grass <- MF_genind_COL[which(extract(COL_GRA, MF_genind_COL@other$xy)>0)]

MF_LUX_urban_summary <- summary(MF_LUX_urban)

mean(MF_LUX_urban_summary$Hexp)
sd(MF_LUX_urban_summary$Hexp)
mean(MF_LUX_urban_summary$Hobs)
sd(MF_LUX_urban_summary$Hobs)

mean(allelicrichness(as.loci(MF_LUX_urban)))
sd(allelicrichness(as.loci(MF_LUX_urban)))

MF_LUX_grass_summary <- summary(MF_LUX_grass)

mean(MF_LUX_grass_summary$Hexp)
sd(MF_LUX_grass_summary$Hexp)
mean(MF_LUX_grass_summary$Hobs)
sd(MF_LUX_grass_summary$Hobs)

mean(allelicrichness(as.loci(MF_LUX_grass)))
sd(allelicrichness(as.loci(MF_LUX_grass)))

MF_COL_urban_summary <- summary(MF_COL_urban)

mean(MF_COL_urban_summary$Hexp)
sd(MF_COL_urban_summary$Hexp)
mean(MF_COL_urban_summary$Hobs)
sd(MF_COL_urban_summary$Hobs)

mean(allelicrichness(as.loci(MF_COL_urban)))
sd(allelicrichness(as.loci(MF_COL_urban)))

MF_COL_grass_summary <- summary(MF_COL_grass)

mean(MF_COL_grass_summary$Hexp)
sd(MF_COL_grass_summary$Hexp)
mean(MF_COL_grass_summary$Hobs)
sd(MF_COL_grass_summary$Hobs)

mean(allelicrichness(as.loci(MF_COL_grass)))
sd(allelicrichness(as.loci(MF_COL_grass)))

mean(allelicrichness(as.loci(SP_genind)))
sd(allelicrichness(as.loci(SP_genind)))

################################################################################
SP_UR_LUX <- SP_genind_LUX[c(which(extract(LUX_IMP, SP_genind_LUX@other$xy)>90),
                           which(extract(LUX_GRA, SP_genind_LUX@other$xy)>0))]

URpops <- rep("lol", 314)
URpops[which(extract(LUX_IMP, SP_UR_LUX@other$xy)>90)] <- "urban"
URpops[which(extract(LUX_GRA, SP_UR_LUX@other$xy)>0)] <- "rural"
URpops <- as.factor(URpops)
  
pop(SP_UR_LUX)<-URpops
SP_LUX.urban.hierf<-genind2hierfstat(SP_UR_LUX)
allrich<-allelic.richness(SP_UR_LUX)
t.test(allrich$Ar$rural, allrich$Ar$urban, pair=T, var.equal = T, alter="greater")
mean(allrich$Ar$rural)
mean(allrich$Ar$urban)

poppr(SP_UR_LUX)
Hobs <- t(sapply(seppop(SP_UR_LUX), function(ls) summary(ls)$Hobs))
Hexp <- t(sapply(seppop(SP_UR_LUX), function(ls) summary(ls)$Hexp))
cat("Expected heterozygosity (Hexp):", "\n")
round(Hexp, 2)
cat("\n", "Observed heterozygosity (Hobs):", "\n")
round(Hobs, 2)

rowMeans(Hexp)
rowMeans(Hobs)
t.test(Hobs[1,], Hobs[2,], pair=T, var.equal = T, alter="greater")
t.test(Hexp[1,], Hexp[2,], pair=T, var.equal = T, alter="greater")

#using GENEPOP to test for differentiation between urban and rural
writeGenPop(SP_UR_LUX, "SP_LUX_filtered.gen", comment="infile")
test_diff("SP_LUX_filtered.gen")

################################################################################
SP_UR_COL <- SP_genind_COL[c(which(extract(COL_IMP, SP_genind_COL@other$xy)>90),
                             which(extract(COL_GRA, SP_genind_COL@other$xy)>0))]

URpops <- rep("lol", 258)
URpops[which(extract(COL_IMP, SP_UR_COL@other$xy)>90)] <- "urban"
URpops[which(extract(COL_GRA, SP_UR_COL@other$xy)>0)] <- "rural"
URpops <- as.factor(URpops)

pop(SP_UR_COL)<-URpops
SP_COL.urban.hierf<-genind2hierfstat(SP_UR_COL)
allrich<-allelic.richness(SP_UR_COL)
t.test(allrich$Ar$rural, allrich$Ar$urban, pair=T, var.equal = T, alter="greater")
mean(allrich$Ar$rural)
mean(allrich$Ar$urban)

poppr(SP_UR_COL)
Hobs <- t(sapply(seppop(SP_UR_COL), function(ls) summary(ls)$Hobs))
Hexp <- t(sapply(seppop(SP_UR_COL), function(ls) summary(ls)$Hexp))
cat("Expected heterozygosity (Hexp):", "\n")
round(Hexp, 2)
cat("\n", "Observed heterozygosity (Hobs):", "\n")
round(Hobs, 2)

rowMeans(Hexp)
rowMeans(Hobs)
t.test(Hobs[1,], Hobs[2,], pair=T, var.equal = T, alter="greater") # here
t.test(Hexp[1,], Hexp[2,], pair=T, var.equal = T, alter="greater")

#using GENEPOP to test for differentiation between urban and rural
writeGenPop(SP_UR_COL, "SP_COL_filtered.gen", comment="infile")
test_diff("SP_COL_filtered.gen")

################################################################################
MF_UR_LUX <- MF_genind_LUX[c(which(extract(LUX_IMP, MF_genind_LUX@other$xy)>90),
                             which(extract(LUX_GRA, MF_genind_LUX@other$xy)>0))]

URpops <- rep("lol", 75)
URpops[which(extract(LUX_IMP, MF_UR_LUX@other$xy)>90)] <- "urban"
URpops[which(extract(LUX_GRA, MF_UR_LUX@other$xy)>0)] <- "rural"
URpops <- as.factor(URpops)

pop(MF_UR_LUX)<-URpops
MF_LUX.urban.hierf<-genind2hierfstat(MF_UR_LUX)
allrich<-allelic.richness(MF_UR_LUX)
t.test(allrich$Ar$rural, allrich$Ar$urban, pair=T, var.equal = T, alter="greater")
mean(allrich$Ar$rural)
mean(allrich$Ar$urban)

poppr(MF_UR_LUX)
Hobs <- t(sapply(seppop(MF_UR_LUX), function(ls) summary(ls)$Hobs))
Hexp <- t(sapply(seppop(MF_UR_LUX), function(ls) summary(ls)$Hexp))
cat("Expected heterozygosity (Hexp):", "\n")
round(Hexp, 2)
cat("\n", "Observed heterozygosity (Hobs):", "\n")
round(Hobs, 2)

rowMeans(Hexp)
rowMeans(Hobs)
t.test(Hobs[1,], Hobs[2,], pair=T, var.equal = T, alter="greater")
t.test(Hexp[1,], Hexp[2,], pair=T, var.equal = T, alter="greater")

#using GENEPOP to test for differentiation between urban and rural
writeGenPop(MF_UR_LUX, "MF_LUX_filtered.gen", comment="infile")
test_diff("MF_LUX_filtered.gen")

################################################################################
MF_UR_COL <- MF_genind_COL[c(which(extract(COL_IMP, MF_genind_COL@other$xy)>90),
                             which(extract(COL_GRA, MF_genind_COL@other$xy)>0))]

URpops <- rep("lol", 138)
URpops[which(extract(COL_IMP, MF_UR_COL@other$xy)>90)] <- "urban"
URpops[which(extract(COL_GRA, MF_UR_COL@other$xy)>0)] <- "rural"
URpops <- as.factor(URpops)

pop(MF_UR_COL)<-URpops
MF_COL.urban.hierf<-genind2hierfstat(MF_UR_COL)
allrich<-allelic.richness(MF_UR_COL)
t.test(allrich$Ar$rural, allrich$Ar$urban, pair=T, var.equal = T, alter="greater")
mean(allrich$Ar$rural)
mean(allrich$Ar$urban)

poppr(MF_UR_COL)
Hobs <- t(sapply(seppop(MF_UR_COL), function(ls) summary(ls)$Hobs))
Hexp <- t(sapply(seppop(MF_UR_COL), function(ls) summary(ls)$Hexp))
cat("Expected heterozygosity (Hexp):", "\n")
round(Hexp, 2)
cat("\n", "Observed heterozygosity (Hobs):", "\n")
round(Hobs, 2)

rowMeans(Hexp)
rowMeans(Hobs)
t.test(Hobs[1,], Hobs[2,], pair=T, var.equal = T, alter="greater") # here
t.test(Hexp[1,], Hexp[2,], pair=T, var.equal = T, alter="greater")

#using GENEPOP to test for differentiation between urban and rural
writeGenPop(MF_UR_COL, "MF_COL_filtered.gen", comment="infile")
test_diff("MF_COL_filtered.gen")

