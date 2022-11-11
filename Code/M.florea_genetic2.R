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
mean(MF_genind_summary$Hobs)

### Allelic richness
barplot(allelicrichness(as.loci(MF_genind)), beside = TRUE)

### Is mean observed H significantly lower than mean expected H ?
par(mfrow = c(1,1))
bartlett.test(list(MF_genind_summary$Hexp,MF_genind_summary$Hobs))
t.test(MF_genind_summary$Hexp,MF_genind_summary$Hobs,pair=T,var.equal=TRUE,alter="greater")
# RESULT: Yes! (mean difference is reasonable - similar ton other data sets I have seen)

### Hardy-Weinberg equilibrium

HW_MF_LU <- hw.test(MF_genind[MF_genind@pop=="LU"], B=1000000)
HW_MF_CO <- hw.test(MF_genind[MF_genind@pop=="CO"], B=1000000)

p.adjust(HW_MF_LU[,4], method="BH")

names(which(p.adjust(HW_MF_LU[,4], method="BH") < 0.05 & p.adjust(HW_MF_CO[,4], method="BH") < 0.05))
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

E_MF_LU <- as.vector(empirLoiselle_EcoGenetics_MF_LU)
E_MF_LU[is.na(E_MF_LU)] <- 0

dens <- kde2d(as.vector(empir_geo_dist_MF_LU), E_MF_LU, n=300)
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(empir_geo_dist_MF_LU, empirLoiselle_EcoGenetics_MF_LU, pch=20,cex=.5)
image(dens, col=tranMF(myPal(300),.7), add=TRUE)
dist_lm <- lm(as.vector(empirLoiselle_EcoGenetics_MF_LU) ~ as.vector(empir_geo_dist_MF_LU))
abline(dist_lm)
title("Isolation by distance plot")

MF_genind_CO <- MF_genind[MF_genind@pop=="CO"]
geo_jittered <- apply(as.data.frame(MF_genind[MF_genind@pop=="CO"]@other$xy), 2, FUN= function(x) jitter(x,amount = 1000))
MF_genind_CO@other$xy <- geo_jittered
MF_genind_CO@tab <- tab(MF_genind_CO, freq = TRUE, NA.method = "mean")

myMFca <- MFca(MF_genind_CO, type=2,ask=FALSE,scannf=FALSE)
barplot(myMFca$eig,main="Eigenvalues of MFCA", col=rep(c("red","grey"),c(1,100)))
barplot(myMFca$eig, main="A variant of the plot\n of MFCA eigenvalues",col=MFectral(length(myMFca$eig)))
legend("topright", fill=MFectral(2),leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")
screeplot(myMFca)
MF_genind_CO@tab <- tab(MF_genind_CO, freq = TRUE, NA.method = "mean")
myGtest <- global.rtest(MF_genind_CO$tab, myMFca$lw, nperm=999)
myGtest
myLtest <- local.rtest(MF_genind_CO$tab,myMFca$lw, nperm=999)
myLtest

plot(myMFca)
colorplot(myMFca,cex=3,main="colorplot of myMFca, first global score")

library(akima)
x <- other(MF_genind_CO)$xy[,1]
y <- other(MF_genind_CO)$xy[,2]
temp <- interp(x, y, myMFca$li[,1])
image(temp, col=azur(100))
points(x,y)
interpX <- seq(min(x),max(x),le=200)
interpY <- seq(min(y),max(y),le=200)
temp <- interp(x, y, myMFca$ls[,1], xo=interpX, yo=interpY)
image(temp, col=azur(100))
points(x,y)


myPal <- colorRampPalette(c("firebrick2", "white", "lightslateblue"))
annot <- function(){
  title("MFCA - interpolated map of individual scores")
  points(x,y)
}

filled.contour(temp, color.pal=myPal, nlev=50,key.title=title("lagged\nscore 1"), plot.title=annot())


### LU
MF_genind_LU <- MF_genind[MF_genind@pop=="SW"|MF_genind@pop=="LU"]
geo_jittered <- apply(as.data.frame(MF_genind_LU@other$xy), 2, FUN= function(x) jitter(x,amount = 500))
MF_genind_LU@other$xy <- geo_jittered
MF_genind_LU@tab <- tab(MF_genind_LU, freq = TRUE, NA.method = "mean")

myMFca <- MFca(MF_genind_LU, type=2,ask=FALSE,scannf=FALSE)
barplot(myMFca$eig,main="Eigenvalues of MFCA", col=rep(c("red","grey"),c(1,100)))
barplot(myMFca$eig, main="A variant of the plot\n of MFCA eigenvalues",col=MFectral(length(myMFca$eig)))
legend("topright", fill=MFectral(2),leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")
screeplot(myMFca)
MF_genind_CO@tab <- tab(MF_genind_LU, freq = TRUE, NA.method = "mean")
myGtest <- global.rtest(MF_genind_LU$tab, myMFca$lw, nperm=999)
myGtest
myLtest <- local.rtest(MF_genind_LU$tab, myMFca$lw, nperm=999)
myLtest

plot(myMFca)
colorplot(myMFca,cex=3,main="colorplot of myMFca, first global score")

library(akima)
x <- other(MF_genind_LU)$xy[,1]
y <- other(MF_genind_LU)$xy[,2]
temp <- interp(x, y, myMFca$li[,1])
image(temp, col=azur(100))
points(x,y)
interpX <- seq(min(x),max(x),le=200)
interpY <- seq(min(y),max(y),le=200)
temp <- interp(x, y, myMFca$ls[,1], xo=interpX, yo=interpY)
image(temp, col=azur(100))
points(x,y)


myPal <- colorRampPalette(c("firebrick2", "white", "lightslateblue"))
annot <- function(){
  title("MFCA - interpolated map of individual scores")
  points(x,y)
}

filled.contour(temp, color.pal=myPal, nlev=50,key.title=title("lagged\nscore 1"), plot.title=annot())