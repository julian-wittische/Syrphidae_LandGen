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

################################################################################
############### Basic exploration
source("Code/GenDataPrep.R")

SP_genind_summary <- summary(SP_genind)

### Alleles and sample sizes
par(mfrow=c(2,2))
plot(SP_genind_summary$n.by.pop, SP_genind_summary$pop.n.all, xlab="Study area sample size",ylab="Number of alleles",main="Alleles numbers and sample sizes",type="n")
text(SP_genind_summary$n.by.pop,SP_genind_summary$pop.n.all,lab=names(SP_genind_summary$n.by.pop))
barplot(SP_genind_summary$loc.n.all, ylab="Number of alleles",main="Number of alleles per locus")
barplot(SP_genind_summary$Hexp-SP_genind_summary$Hobs, main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")
barplot(SP_genind_summary$n.by.pop, main="Sample sizes per population",ylab="Number of genotypes",las=3)
barplot(SP_genind_summary$Hexp, main="Heterozygosity: expected",ylab="Hexp")
barplot(SP_genind_summary$Hobs, main="Heterozygosity: observed",ylab="Hobs")

mean(SP_genind_summary$Hexp)
mean(SP_genind_summary$Hobs)

### Allelic richness
barplot(allelicrichness(as.loci(SP_genind)), beside = TRUE)

### Is mean observed H significantly lower than mean expected H ?
par(mfrow = c(1,1))
bartlett.test(list(SP_genind_summary$Hexp,SP_genind_summary$Hobs))
t.test(SP_genind_summary$Hexp,SP_genind_summary$Hobs,pair=T,var.equal=TRUE,alter="greater")
# RESULT: Yes! (mean difference is reasonable - similar ton other data sets I have seen)

### Hardy-Weinberg equilibrium
hw.test(SP_genind, B=1000)
HW_LU <- hw.test(SP_genind[SP_genind@pop=="LU"], B=1000)
HW_CO <- hw.test(SP_genind[SP_genind@pop=="CO"], B=1000)

names(which(HW_LU[,4] < 0.05 & HW_CO[,4] < 0.05))
# RESULT: Spp142, Spp051, Spp108, Spp141, Spp360  are somewhat in disequilibrium
# but let us keep them because Alain has subsampled 30 individuals 10 times and they pass

## Weir and Cockerham F statistics

wc(SP_genind)
# RESULT: global Fst is very low; Fis is rather high

# Let's build confidence intervals for that
SP_g2h <- genind2hierfstat(SP_genind)
boot.vc(SP_g2h[1], SP_g2h[-1])$ci
# RESULT: Fst is very low and 0 is within the CI

# Let's check per locus
Fperlocus <- Fst(as.loci(SP_genind))
Fperlocus # Spp141 has a high Fis value 0.24
colMeans(Fperlocus) # Check: similar to calculations above

# Pairwise Fst
Fst <- genet.dist(SP_genind, method = "Nei87")
Fst
# RESULT: SW and LU are closer to each other than they are from CO
is.euclid(Fst) #FALSE because of missing values

###### Check linkage disequilibrium

SP_LD_CHECK <- poppr::ia(SP_genind, sample=199)
SP_LD_CHECK
# There is statistically significant association among the markers but the overall correlation is very low.
SP_LD_pair <- poppr::pair.ia(SP_genind)
SP_LD_pair
# RESULT: /!\ HIGH /!\ between Spp141 and Spp051 ! Those two are also somewhat in HWE.

###### Check null alleles

library(PopGenReport)
Null.alleles <- PopGenReport::null.all(SP_genind)

{cat(" summary1 (Chakraborty et al. 1994):", "\n")
  round(Null.alleles$null.allele.freq$summary1, 2)}

{cat("summary2 (Brookfield et al. 1996):", "\n")
  round(Null.alleles$null.allele.freq$summary2, 2)}

# We cannot exclude that Spp142, Spp051, Spp108, Spp141, Spp313, Spp360 and Spp416 bear null alleles
# Spp141 is the highest here, again
# https://doi.org/10.1093/molbev/msl191 suggests we could keep Spp108 and Spp360
# https://dx.doi.org/10.7717%2Fpeerj.3188 warns about the relationship between
# null alleles and genetic structure
# NOTE: Spp142, Spp051, Spp108, Spp141, Spp360 were somewhat in HWE
# False positive possible explanations: incomplete genotypes, genetic structure,
# breeding patterns, differences in allele sizes, and unique alleles that are
# fixed or nearly fixed locally

################################################################################
############### Nonspatial population structure

# Please note that SP is an object without Spp141 (see GenDataPrep.R)
colours(1)

######## DAPC

##### Prior populations: Luxembourg and Cologne
SPxval1_priorpop <- xvalDapc(tab(SP, NA.method = "mean"), pop(SP),
                             n.pca = 1:100, n.rep = 100,
                             parallel = "snow", ncpus = 4L)

# highest success: 62; lowest MSE: 53; consider +-8
SPxval2_priorpop <- xvalDapc(tab(SP, NA.method = "mean"), pop(SP),
                             n.pca = 45:70, n.rep = 1000,
                             parallel = "snow", ncpus = 4L)
# result: 50 is the best

dapc_priorpop <- dapc(SP, pop(SP), n.pca=50, n.da=1)

scatter(dapc_priorpop, col = rainbow(as.numeric(pop(SP))), cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

scatter(dapc_priorpop, ratio.pca=0.3, bg="white", pch=20, cell=0 ,cstar=0, col=rainbow(2),
        solid=.4, cex=3, clab=0,mstree=FALSE, scree.da=FALSE,
        posi.pca="bottomright", leg=TRUE, txt.leg=paste("Cluster",1:2))


##### De novo populations
grp <- find.clusters(SP, max.n.clust=40, n.iter=1000000, n.pca=1000)
16

grp <- find.clusters(SP, max.n.clust=40, n.iter=1000000, n.pca=1000,
                     choose.n.clust=FALSE, criterion="diffNgroup")

grpres <- numeric(100)

for (i in 1:100){
  grp <- find.clusters(SP, max.n.clust=40, n.iter=100000, n.pca=1000,
                       choose.n.clust=FALSE, criterion="diffNgroup")
  grpres[i] <- names(grp$stat) 
}


SPxval1_denovopop <- xvalDapc(tab(SP, NA.method = "mean"), grp$grp,
                             n.pca = 1:100, n.rep = 100,
                             parallel = "snow", ncpus = 4L)

SPxval1_denovopop <- xvalDapc(tab(SP, NA.method = "mean"), grp$grp,
                              n.pca = 1:100, n.rep = 1000,
                              parallel = "snow", ncpus = 4L)




dapc_denovopop <- dapc(SP, grp$grp)
73
4



scatter(pramx$DAPC, col = as.numeric(pop(SP)), cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)

plot(SP@other$xy, col=dapc1$grp)

scatter(dapc1)

myCol <- c("darkblue", "purple", "green", "orange", "red", "blue")
scatter(dapc1, ratio.pca=0.3, bg="white", pch=20, cell=0 ,cstar=0, col=myCol,
        solid=.4, cex=3, clab=0,mstree=FALSE, scree.da=FALSE,
        posi.pca="bottomright", leg=TRUE, txt.leg=paste("Cluster",1:5))
par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,cex=3, lwd=8, col="black")
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,cex=3, lwd=2, col=myCol)
myInset <- function(){
  temp <- dapc1$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"), c(dapc1$n.pca,1000)),
       ylim=c(0,100), xlab="PCA axis", ylab="Cumulated variance (%)",cex=1,
       pch=20, type="h", lwd=2)
}
add.scatter(myInset(), posi="bottomright",inset=c(-0.03,-0.01), ratio=.28,
            bg=transp("white"))

scatter(dapc1,1,1, col=myCol, bg="white",scree.da=FALSE, legend=TRUE, solid=.4)

set.seed(4)
contrib <- loadingplot(dapc1$var.contr, axis=2, thres=.07, lab.jitter=1)
# SP323 was not among the loci deviating a bit from HWE

round(head(dapc1$posterior),3)
summary(dapc1)
assignplot(dapc1, subset=1:50)
compoplot(dapc1, subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:3), lab="",ncol=2, xlab="individuals")

# a-score analysis
dapc2 <- dapc(SP_genind, n.da=100, n.pca=10)

## DAPC by study area
# SWLU
grp_SWLU <- find.clusters(SP_genind_noSpp141_360_080[SP_genind@pop=="SW"|SP_genind@pop=="LU"], max.n.clust=40)
1000
6

dapc_SWLU <- dapc(SP_genind_noSpp141_360_080[SP_genind@pop=="SW"|SP_genind@pop=="LU"], grp_SWLU$grp)
50
5

plot(as.data.frame(SP_genind_noSpp141_360_080[SP_genind@pop=="SW"|SP_genind@pop=="LU"]@other$xy), col=dapc_SWLU$grp)

scatter(dapc_SWLU)

myCol <- c("darkblue", "purple", "green", "orange", "red", "blue")
scatter(dapc_SWLU, ratio.pca=0.3, bg="white", pch=20, cell=0 ,cstar=0, col=myCol,
        solid=.4, cex=3, clab=0,mstree=TRUE, scree.da=FALSE,
        posi.pca="bottomright", leg=TRUE, txt.leg=paste("Cluster",1:6))
par(xpd=TRUE)
points(dapc_SWLU$grp_SWLU.coord[,1], dapc_SWLU$grp_SWLU.coord[,2], pch=4,cex=3, lwd=8, col="black")
points(dapc_SWLU$grp_SWLU.coord[,1], dapc_SWLU$grp_SWLU.coord[,2], pch=4,cex=3, lwd=2, col=myCol)
myInset <- function(){
  temp <- dapc_SWLU$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"), c(dapc_SWLU$n.pca,1000)),
       ylim=c(0,100), xlab="PCA axis", ylab="Cumulated variance (%)",cex=1,
       pch=20, type="h", lwd=2)
}
add.scatter(myInset(), posi="bottomright",inset=c(-0.03,-0.01), ratio=.28,
            bg=transp("white"))

scatter(dapc_SWLU,1,1, col=myCol, bg="white",scree.da=FALSE, legend=TRUE, solid=.4)

set.seed(4)
contrib <- loadingplot(dapc_SWLU$var.contr, axis=2, thres=.07, lab.jitter=1)

round(head(dapc_SWLU$posterior),3)
summary(dapc_SWLU)
assignplot(dapc_SWLU, subset=1:50)
compoplot(dapc_SWLU, subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:6), lab="",ncol=2, xlab="individuals")

# CO
grp_CO <- find.clusters(SP_genind_noSpp141_360_080[SP_genind_noSpp141_360_080@pop=="CO"], max.n.clust=40)
1000
6

dapc_CO <- dapc(SP_genind_noSpp141_360_080[SP_genind_noSpp141_360_0801@pop=="CO"], grp_CO$grp)
40
5

plot(as.data.frame(SP_genind_noSpp141_360_080[SP_genind_noSpp141_360_080@pop=="CO"]@other$xy), col=dapc_CO$grp)

scatter(dapc_CO)

myCol <- c("darkblue", "purple", "green", "orange", "red", "blue")
scatter(dapc_CO, ratio.pca=0.3, bg="white", pch=20, cell=0 ,cstar=0, col=myCol,
        solid=.4, cex=3, clab=0,mstree=TRUE, scree.da=FALSE,
        posi.pca="bottomright", leg=TRUE, txt.leg=paste("Cluster",1:6))
par(xpd=TRUE)
points(dapc_CO$grp_CO.coord[,1], dapc_CO$grp_CO.coord[,2], pch=4,cex=3, lwd=8, col="black")
points(dapc_CO$grp_CO.coord[,1], dapc_CO$grp_CO.coord[,2], pch=4,cex=3, lwd=2, col=myCol)
myInset <- function(){
  temp <- dapc_CO$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"), c(dapc_CO$n.pca,1000)),
       ylim=c(0,100), xlab="PCA axis", ylab="Cumulated variance (%)",cex=1,
       pch=20, type="h", lwd=2)
}
add.scatter(myInset(), posi="bottomright",inset=c(-0.03,-0.01), ratio=.28,
            bg=transp("white"))

scatter(dapc_CO,1,1, col=myCol, bg="white",scree.da=FALSE, legend=TRUE, solid=.4)

set.seed(4)
contrib <- loadingplot(dapc_CO$var.contr, axis=2, thres=.07, lab.jitter=1)
# SP323 was not among the loci deviating a bit from HWE

round(head(dapc_CO$posterior),3)
summary(dapc_CO)
assignplot(dapc_CO, subset=1:50)
compoplot(dapc_CO, subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:6), lab="",ncol=2, xlab="individuals")

plot(as.data.frame(SP_genind_noSpp141_360_080[SP_genind_noSpp141_360_080@pop=="CO"]@other$xy),
          col=dapc_CO$grp, cex=1, asp=1, pch=19)
geo <- as.data.frame(SP_genind_noSpp141_360_080[SP_genind_noSpp141_360_080@pop=="CO"]@other$xy)

#poppr::poppr.amova(SP_genind_noSpp141_360_080, method="pegas")

################################################################################
############### Spatial analyses

### Classic IBD

empir_geo_dist_SP <- as.matrix(dist(as.data.frame(SP@other$xy)))
empir_geo_dist_SP_LU <- as.matrix(dist(as.data.frame(SP[SP@pop=="LU"]@other$xy)))
empir_geo_dist_SP_CO <- as.matrix(dist(as.data.frame(SP[SP@pop=="CO"]@other$xy)))

empir_geo_dist_SP2 <-empir_geo_dist_SP
empir_geo_dist_SP_LU2 <-empir_geo_dist_SP_LU
empir_geo_dist_SP_CO2 <-empir_geo_dist_SP_CO

empir_geo_dist_SP2[empir_geo_dist_SP2==0] <- NA
empir_geo_dist_SP_LU2[empir_geo_dist_SP_LU2==0] <- NA
empir_geo_dist_SP_CO2[empir_geo_dist_SP_CO2==0] <- NA

SP_nogeo <- SP
SP_nogeo@other <- NULL
SP_nogeo@pop <- NULL

empirLoiselle_EcoGenetics_SP <- eco.kin.loiselle(genind2ecogen(SP_nogeo))
empirLoiselle_EcoGenetics_SP_LU <- eco.kin.loiselle(genind2ecogen(SP_nogeo[SP@pop=="LU"]))
empirLoiselle_EcoGenetics_SP_CO <- eco.kin.loiselle(genind2ecogen(SP_nogeo[SP@pop=="CO"]))

IBD_SP <- lm(c(as.dist(empirLoiselle_EcoGenetics_SP))~log(c(as.dist(empir_geo_dist_SP2))))
IBD_SP_LU <- lm(c(as.dist(empirLoiselle_EcoGenetics_SP_LU))~log(c(as.dist(empir_geo_dist_SP_LU2))))
IBD_SP_CO <- lm(c(as.dist(empirLoiselle_EcoGenetics_SP_CO))~log(c(as.dist(empir_geo_dist_SP_CO2))))

summary(IBD_SP)
summary(IBD_SP_LU)
summary(IBD_SP_CO)

mantel.randtest(as.dist(empir_geo_dist_SP), as.dist(1-empirLoiselle_EcoGenetics_SP), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_SP_LU), as.dist(1-empirLoiselle_EcoGenetics_SP_LU), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_SP_CO), as.dist(1-empirLoiselle_EcoGenetics_SP_CO), nrepet = 9999)

plot(mantel.randtest(as.dist(empir_geo_dist_SP), as.dist(1-empirLoiselle_EcoGenetics_SP), nrepet = 9999))
plot(mantel.randtest(as.dist(empir_geo_dist_SP_CO), as.dist(1-empirLoiselle_EcoGenetics_SP_CO), nrepet = 9999))
plot(mantel.randtest(as.dist(empir_geo_dist_SP_LU), as.dist(1-empirLoiselle_EcoGenetics_SP_LU), nrepet = 9999))

plot(log(empir_geo_dist_SP2), empirLoiselle_EcoGenetics_SP)
abline(IBD_SP, col="red")
plot(log(empir_geo_dist_SP_LU2), empirLoiselle_EcoGenetics_SP_LU)
abline(IBD_SP_LU, col="red")
plot(log(empir_geo_dist_SP_CO2), empirLoiselle_EcoGenetics_SP_CO)
abline(IBD_SP_CO, col="red")

E_SP_LU <- as.vector(empirLoiselle_EcoGenetics_SP_LU)
E_SP_LU[is.na(E_SP_LU)] <- 0

dens <- kde2d(as.vector(empir_geo_dist_SP_LU), E_SP_LU, n=300)
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(empir_geo_dist_SP_LU, empirLoiselle_EcoGenetics_SP_LU, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
dist_lm <- lm(as.vector(empirLoiselle_EcoGenetics_SP_LU) ~ as.vector(empir_geo_dist_SP_LU))
abline(dist_lm)
title("Isolation by distance plot")

SP_genind_CO <- SP_genind[SP_genind@pop=="CO"]
geo_jittered <- apply(as.data.frame(SP_genind[SP_genind@pop=="CO"]@other$xy), 2, FUN= function(x) jitter(x,amount = 1000))
SP_genind_CO@other$xy <- geo_jittered
SP_genind_CO@tab <- tab(SP_genind_CO, freq = TRUE, NA.method = "mean")

mySpca <- spca(SP_genind_CO, type=2,ask=FALSE,scannf=FALSE)
barplot(mySpca$eig,main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(1,100)))
barplot(mySpca$eig, main="A variant of the plot\n of sPCA eigenvalues",col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")
screeplot(mySpca)
SP_genind_CO@tab <- tab(SP_genind_CO, freq = TRUE, NA.method = "mean")
myGtest <- global.rtest(SP_genind_CO$tab, mySpca$lw, nperm=999)
myGtest
myLtest <- local.rtest(SP_genind_CO$tab,mySpca$lw, nperm=999)
myLtest

plot(mySpca)
colorplot(mySpca,cex=3,main="colorplot of mySpca, first global score")

library(akima)
x <- other(SP_genind_CO)$xy[,1]
y <- other(SP_genind_CO)$xy[,2]
temp <- interp(x, y, mySpca$li[,1])
image(temp, col=azur(100))
points(x,y)
interpX <- seq(min(x),max(x),le=200)
interpY <- seq(min(y),max(y),le=200)
temp <- interp(x, y, mySpca$ls[,1], xo=interpX, yo=interpY)
image(temp, col=azur(100))
points(x,y)


myPal <- colorRampPalette(c("firebrick2", "white", "lightslateblue"))
annot <- function(){
  title("sPCA - interpolated map of individual scores")
  points(x,y)
}

filled.contour(temp, color.pal=myPal, nlev=50,key.title=title("lagged\nscore 1"), plot.title=annot())


### LU
SP_genind_LU <- SP_genind[SP_genind@pop=="SW"|SP_genind@pop=="LU"]
geo_jittered <- apply(as.data.frame(SP_genind_LU@other$xy), 2, FUN= function(x) jitter(x,amount = 500))
SP_genind_LU@other$xy <- geo_jittered
SP_genind_LU@tab <- tab(SP_genind_LU, freq = TRUE, NA.method = "mean")

mySpca <- spca(SP_genind_LU, type=2,ask=FALSE,scannf=FALSE)
barplot(mySpca$eig,main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(1,100)))
barplot(mySpca$eig, main="A variant of the plot\n of sPCA eigenvalues",col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")
screeplot(mySpca)
SP_genind_CO@tab <- tab(SP_genind_LU, freq = TRUE, NA.method = "mean")
myGtest <- global.rtest(SP_genind_LU$tab, mySpca$lw, nperm=999)
myGtest
myLtest <- local.rtest(SP_genind_LU$tab, mySpca$lw, nperm=999)
myLtest

plot(mySpca)
colorplot(mySpca,cex=3,main="colorplot of mySpca, first global score")

library(akima)
x <- other(SP_genind_LU)$xy[,1]
y <- other(SP_genind_LU)$xy[,2]
temp <- interp(x, y, mySpca$li[,1])
image(temp, col=azur(100))
points(x,y)
interpX <- seq(min(x),max(x),le=200)
interpY <- seq(min(y),max(y),le=200)
temp <- interp(x, y, mySpca$ls[,1], xo=interpX, yo=interpY)
image(temp, col=azur(100))
points(x,y)


myPal <- colorRampPalette(c("firebrick2", "white", "lightslateblue"))
annot <- function(){
  title("sPCA - interpolated map of individual scores")
  points(x,y)
}

filled.contour(temp, color.pal=myPal, nlev=50,key.title=title("lagged\nscore 1"), plot.title=annot())