################################################################################
########## Julian Wittische - July 2022 - Hoverfly landscape genetics ##########
################################################################################
#
# SCRIPT OBJECTIVE:
#
# - load genetic data and change it to a readable format
#
# - perform basic popgen analyses, IBD, PCA, and DAPC
#_______________________________________________________________________________

################################################################################
############### Loading the INCOMPLETE S. pipiens microsatellite dataset

library(adegenet)
library(pegas)
library(hierfstat)
library(sp)
library(EcoGenetics)
library(MASS)
library(vegan)

################################################################################
############### Basic exploration
source("Code/GenDataPrep.R")

MF_genind_summary <- summary(MF_genind)

### Alleles and sample sizes
par(mfrow=c(2,2))
plot(MF_genind_summary$n.by.pop, MF_genind_summary$pop.n.all, xlab="Study area sample size",ylab="Number of alleles",main="Alleles numbers and sample sizes",type="n")
text(MF_genind_summary$n.by.pop,MF_genind_summary$pop.n.all,lab=names(MF_genind_summary$n.by.pop))
barplot(MF_genind_summary$loc.n.all, ylab="Number of alleles",main="Number of alleles per locus")
barplot(MF_genind_summary$Hexp-MF_genind_summary$Hobs, main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")
barplot(MF_genind_summary$n.by.pop, main="Sample sizes per population",ylab="Number of genotypes",las=3)
barplot(MF_genind_summary$Hexp, main="Heterozygosity: expected",ylab="Hexp")
barplot(MF_genind_summary$Hobs, main="Heterozygosity: observed",ylab="Hobs")

mean(MF_genind_summary$Hexp)
sd(MF_genind_summary$Hexp)
mean(MF_genind_summary$Hobs)
sd(MF_genind_summary$Hobs)

### Allelic richness
barplot(allelicrichness(as.loci(MF_genind)), beside = TRUE)

### Is mean observed H significantly lower than mean expected H ?
par(mfrow = c(1,1))
bartlett.test(list(MF_genind_summary$Hexp,MF_genind_summary$Hobs))
t.test(MF_genind_summary$Hexp,MF_genind_summary$Hobs,pair=T,var.equal=TRUE,alter="greater")
# RESULT: Yes! (mean difference is reasonable - similar ton other data sets I have seen)

### Hardy-Weinberg equilibrium
hw.test(MF_genind, B=1000)
HW_LU <- hw.test(MF_genind[MF_genind@pop=="LU"], B=1000)
HW_CO <- hw.test(MF_genind[MF_genind@pop=="CO"], B=1000)

names(which(HW_LU[,4] < 0.05 & HW_CO[,4] < 0.05))
# RESULT: MFp142, MFp051, MFp108, MFp141, MFp360  are somewhat in disequilibrium
# but let us keep them because Alain has subsampled 30 individuals 10 times and they pass

## Weir and Cockerham F statistics

wc(MF_genind)
# RESULT: global Fst is very low; Fis is rather high

# Let's build confidence intervals for that
MF_g2h <- genind2hierfstat(MF_genind)
boot.vc(MF_g2h[1], MF_g2h[-1])$ci
# RESULT: Fst is very low and 0 is within the CI

# Let's check per locus
Fperlocus <- Fst(as.loci(MF_genind))
Fperlocus # 
colMeans(Fperlocus) # Check: similar to calculations above

# Pairwise Fst
Fst <- genet.dist(MF_genind, method = "Nei87")
Fst
is.euclid(Fst)

###### Check linkage disequilibrium

MF_LD_CHECK <- poppr::ia(MF_genind, sample=199)
MF_LD_CHECK
# There is statistically significant association among the markers but the overall correlation is very low.
MF_LD_pair <- poppr::pair.ia(MF_genind)
MF_LD_pair

###### Check null alleles

library(PopGenReport)
Null.alleles <- PopGenReport::null.all(MF_genind)

{cat(" summary1 (Chakraborty et al. 1994):", "\n")
  round(Null.alleles$null.allele.freq$summary1, 2)}

{cat("summary2 (Brookfield et al. 1996):", "\n")
  round(Null.alleles$null.allele.freq$summary2, 2)}

# We cannot exclude that MFp142, MFp051, MFp108, MFp141, MFp313, MFp360 and MFp416 bear null alleles
# MFp141 is the highest here, again
# https://doi.org/10.1093/molbev/msl191 suggests we could keep MFp108 and MFp360
# https://dx.doi.org/10.7717%2Fpeerj.3188 warns about the relationship between
# null alleles and genetic structure
# NOTE: MFp142, MFp051, MFp108, MFp141, MFp360 were somewhat in HWE
# False positive possible explanations: incomplete genotypes, genetic structure,
# breeding patterns, differences in allele sizes, and unique alleles that are
# fixed or nearly fixed locally

################################################################################
############### NonMFatial population structure

######## DAPC

##### Prior populations: Luxembourg and Cologne
MFxval1_priorpop <- xvalDapc(tab(MF, NA.method = "mean"), pop(MF),
                             n.pca = 1:100, n.rep = 100,
                             parallel = "snow", ncpus = 4L)

# highest success: 62; lowest MSE: 53; consider +-8
MFxval2_priorpop <- xvalDapc(tab(MF, NA.method = "mean"), pop(MF),
                             n.pca = 80:85, n.rep = 1000,
                             parallel = "snow", ncpus = 4L)
# result: 83 is the best

dapc_priorpop_MF <- dapc(MF, pop(MF), n.pca=83, n.da=1)

scatter(dapc_priorpop_MF, col = c("#FF7F00", "#8F00FF"), cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

##### De novo populations
grpres <- numeric(1000)

for (i in 1:1000){
  grp <- find.clusters(MF, max.n.clust=40, n.iter=100, n.pca=1000,
                       choose.n.clust=FALSE, criterion="diffNgroup")
  grpres[i] <- names(grp$stat)
  print(paste(i, "%"))
}
table(grpres)
#6 is best but barely

grp <- find.clusters(MF, max.n.clust=40, n.iter=1000000, n.pca=1000)
6

MFxval1_denovopop <- xvalDapc(tab(MF, NA.method = "mean"), grp$grp,
                              n.pca = 1:100, n.rep = 100,
                              parallel = "snow", ncpus = 4L)

MFxval1_denovopop2 <- xvalDapc(tab(MF, NA.method = "mean"), grp$grp,
                               n.pca = 1:30, n.rep = 1000,
                               parallel = "snow", ncpus = 4L)

dapc_denovopop <- dapc(MF, grp$grp, n.pca=21, n.da=5)

scatter(dapc_denovopop, col = c("#000000", "#FF0000", "#66FF66", "#F17925",
                                "#CCAA14", "#DEA6CF"), cex = 6, legend = TRUE,
        clabel = TRUE, posi.leg = "bottomleft", scree.pca = FALSE,
        posi.pca = "bottomright", cleg = 2, xax = 1, yax = 2, scree.da=FALSE,
        solid=0.5, mstree=FALSE, cstar=1)
lvls <- levels(dapc_denovopop$grp) <- c("#000000", "#FF0000", "#66FF66", "#F17925",
                                        "#CCAA14", "#DEA6CF")

plot(MF@other$xy, col=lvls)
plot(MF@other$xy, col=dapc_denovopop$grp)

################################################################################
############### MFatial analyses

### Classic IBD

empir_geo_dist_MF <- as.matrix(dist(as.data.frame(MF@other$xy)))
empir_geo_dist_MF_LU <- as.matrix(dist(as.data.frame(MF[MF@pop=="LU"]@other$xy)))
empir_geo_dist_MF_CO <- as.matrix(dist(as.data.frame(MF[MF@pop=="CO"]@other$xy)))

empir_geo_dist_MF2 <-empir_geo_dist_MF
empir_geo_dist_MF_LU2 <-empir_geo_dist_MF_LU
empir_geo_dist_MF_CO2 <-empir_geo_dist_MF_CO

empir_geo_dist_MF2[empir_geo_dist_MF2==0] <- NA
empir_geo_dist_MF_LU2[empir_geo_dist_MF_LU2==0] <- NA
empir_geo_dist_MF_CO2[empir_geo_dist_MF_CO2==0] <- NA

MF_nogeo <- MF
MF_nogeo@other <- NULL
MF_nogeo@pop <- NULL

empirLoiselle_EcoGenetics_MF <- eco.kin.loiselle(genind2ecogen(MF_nogeo))
empirLoiselle_EcoGenetics_MF_LU <- eco.kin.loiselle(genind2ecogen(MF_nogeo[MF@pop=="LU"]))
empirLoiselle_EcoGenetics_MF_CO <- eco.kin.loiselle(genind2ecogen(MF_nogeo[MF@pop=="CO"]))

IBD_MF <- lm(c(as.dist(empirLoiselle_EcoGenetics_MF))~log(c(as.dist(empir_geo_dist_MF2))))
IBD_MF_LU <- lm(c(as.dist(empirLoiselle_EcoGenetics_MF_LU))~log(c(as.dist(empir_geo_dist_MF_LU2))))
IBD_MF_CO <- lm(c(as.dist(empirLoiselle_EcoGenetics_MF_CO))~log(c(as.dist(empir_geo_dist_MF_CO2))))

summary(IBD_MF)
summary(IBD_MF_LU)
summary(IBD_MF_CO)

mantel.randtest(as.dist(empir_geo_dist_MF), as.dist(1-empirLoiselle_EcoGenetics_MF), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_MF_LU), as.dist(1-empirLoiselle_EcoGenetics_MF_LU), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_MF_CO), as.dist(1-empirLoiselle_EcoGenetics_MF_CO), nrepet = 9999)

plot(mantel.randtest(as.dist(empir_geo_dist_MF), as.dist(1-empirLoiselle_EcoGenetics_MF), nrepet = 9999))
plot(mantel.randtest(as.dist(empir_geo_dist_MF_CO), as.dist(1-empirLoiselle_EcoGenetics_MF_CO), nrepet = 9999))
plot(mantel.randtest(as.dist(empir_geo_dist_MF_LU), as.dist(1-empirLoiselle_EcoGenetics_MF_LU), nrepet = 9999))

plot(log(empir_geo_dist_MF2), empirLoiselle_EcoGenetics_MF)
abline(IBD_MF, col="red")
plot(log(empir_geo_dist_MF_LU2), empirLoiselle_EcoGenetics_MF_LU)
abline(IBD_MF_LU, col="red")
plot(log(empir_geo_dist_MF_CO2), empirLoiselle_EcoGenetics_MF_CO)
abline(IBD_MF_CO, col="red")

mantel.correlog(empirLoiselle_EcoGenetics_MF,empir_geo_dist_MF)
mantel.correlog(empirLoiselle_EcoGenetics_MF_LU,empir_geo_dist_MF_LU)
mantel.correlog(empirLoiselle_EcoGenetics_MF_CO,empir_geo_dist_MF_CO)
