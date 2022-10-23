################################################################################
########################### Julian Wittische - MNHNL ###########################
################################################################################

################################################################################
############### Loading the INCOMPLETE M. florea microsatellite dataset

library(adegenet)
library(pegas)
library(hierfstat)
library(sp)
library(EcoGenetics)
library(MASS)

# Read
MF_df_raw <- as.data.frame(readxl::read_excel("Data/Mflorea_raw.xlsx",
                                              sheet=1, .name_repair="minimal"))

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
# row.names(MF_df)[which(duplicated(row.names(MF_df)))] <- "LUX_98_2_08-09-21_A_bis" # CHANGED IN RAW DATA

# Transform into genind object
MF_genind <- df2genind(MF_df, ncode=3, ploidy=2, sep="/",
                       type="codom", NA.char=" NA/ NA")

# Add pop information (study area)
MF_genind@pop <- as.factor(substr(row.names(MF_df), 1, 2))

# Add coordinates and reproject them
MF_geo_sp <- SpatialPoints(MF_geo, CRS(SRS_string = "EPSG:4326"))
MF_geo_sp <-spTransform(MF_geo_sp, CRS(SRS_string = "EPSG:3035"))
MF_genind@other$xy <- MF_geo_sp

################################################################################
############### Basic exploration
#adegenetTutorial( which = c("basics"))
  
MF_genind_summary <- summary(MF_genind)

### Alleles and sample sizes
par(mfrow=c(2,2))
plot(MF_genind_summary$n.by.pop, MF_genind_summary$pop.n.all, xlab="Study area sample size",ylab="Number of alleles",main="Alleles numbers and sample sizes",type="n")
text(MF_genind_summary$n.by.pop,MF_genind_summary$pop.n.all,lab=names(MF_genind_summary$n.by.pop))
barplot(MF_genind_summary$loc.n.all, ylab="Number of alleles",main="Number of alleles per locus")
barplot(MF_genind_summary$Hexp-MF_genind_summary$Hobs, main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")
barplot(MF_genind_summary$n.by.pop, main="Sample sizes per population",ylab="Number of genotypes",las=3)

### Allelic richness
barplot(allelicrichness(as.loci(MF_genind)), beside = TRUE)

### Is mean observed H significantly lower than mean expected H ?
par(mfrow = c(1,1))
bartlett.test(list(MF_genind_summary$Hexp,MF_genind_summary$Hobs))
t.test(MF_genind_summary$Hexp,MF_genind_summary$Hobs,pair=T,var.equal=TRUE,alter="greater")
# Yes! (same mean difference as in other data sets I have seen)

### Hardy-Weinberg equilibrium
hw.test(MF_genind, B=1000)
HW_LU <- hw.test(MF_genind[MF_genind@pop=="LU"], B=1000)
HW_CO <- hw.test(MF_genind[MF_genind@pop=="CO"], B=1000)
HW_SW <- hw.test(MF_genind[MF_genind@pop=="SW"], B=1000)

which(HW_LU[,4] < 0.05 & HW_CO[,4] < 0.05 & HW_SW[,4] < 0.05)
# MF239, MF303, MF103, MF28, MF70 are somewhat in disequilibrium but let us keep them
# Alain has subsampled 30 individuals 10 times and they pass

################################################################################
############### Nonspatial population structure

## Weir and Cockerham F statistics

wc(MF_genind) # global Fst is pretty low high but Fis is high (weird because no HW des)

# Let's build confidence intervals for that
MF_g2h <- genind2hierfstat(MF_genind)
boot.vc(MF_g2h[1], MF_g2h[-1])$ci

# Let's check per locus
Fperlocus <- Fst(as.loci(MF_genind))
Fperlocus # super high Fis values (0.96!) and very low Fst values
colMeans(Fperlocus)

# Pairwise Fst
Fst <- genet.dist(MF_genind, method = "Nei87")
is.euclid(Fst) #FALSE because of missing values

###### Check linkage disequilibrium

MF_LD_CHECK <- poppr::ia(MF_genind, sample=199)
MF_LD_CHECK
# There is statistically significant association among the markers but the overall correlation is very low.
MF_LD_pair <- poppr::pair.ia(MF_genind)
MF_LD_pair
# RESULT: /!\ HIGH /!\ between Spp141 and Spp051 ! Those two are also somewhat in HWE.



###### Check null alleles

library(PopGenReport)
Null.alleles <- PopGenReport::null.all(MF_genind)

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

### PCA
# Let us replace those NA values
sum(is.na(MF_genind$tab))
PCAdf <- tab(MF_genind, freq = TRUE, NA.method = "mean")

# Check
class(PCAdf)
dim(PCAdf)
PCAdf[1:5,1:5]


PCA <- dudi.pca(PCAdf, scale = FALSE, scannf = FALSE, nf = 3)
# Eigenvalues
barplot(PCA$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

# Viz 1 - General PCA plot
s.label(PCA$li)
title("PCA - M. florea \naxes 1-2")
add.scatter.eig(PCA$eig[1:20], 3, 1, 2)
# Three groups as originally shown by Alain

# Viz 2 - Study area focused
s.class(PCA$li, pop(MF_genind))
title("PCA - M. florea \naxes 1-2")
add.scatter.eig(PCA$eig[1:20], 3, 1, 2)
# The three groups are obviously not geographical

# Viz 3 - Third axis (it does not seem very important)
s.class(PCA$li,pop(MF_genind),xax=1,yax=3,sub="PCA 1-3",csub=2)
title("PCA - M. florea \naxes 1-3")
add.scatter.eig(PCA$eig[1:20], nf=3, xax=1, yax=3)

col <- funky(5)
s.class(PCA$li, pop(MF_genind), xax=1, yax=3, col=transp(col,.6), axesell=FALSE, 
        cstar=0, cpoint=3, grid=FALSE)

colorplot(PCA$li, PCA$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA - M. florea\naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

colorplot(PCA$li[c(1,3)], PCA$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
title("PCA - M. florea\naxes 1-3")
abline(v=0,h=0,col="grey", lty=2)

### CA
obj <- genind2genpop(MF_genind)
ca1 <- dudi.coa(tab(obj),scannf=FALSE,nf=3)
barplot(ca1$eig,main="Correspondance Analysis eigenvalues",col=heat.colors(length(ca1$eig)))
s.label(ca1$li, sub="CA 1-2",csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=1,yax=2,posi="bottomright")
# s.label(ca1$li,xax=2,yax=3,lab=popNames(obj),sub="CA 1-3",csub=2)
# add.scatter.eig(ca1$eig,nf=3,xax=2,yax=3,posi="topleft"

### DAPC
# NOTE: we should do this temporally too
grp <- find.clusters(MF_genind, max.n.clust=40)
1000
6

dapc1 <- dapc(MF_genind, grp$grp)
60
5

plot(MF_genind@other$xy, col=dapc1$grp)

scatter(dapc1)

myCol <- c("darkblue", "purple", "green", "orange", "red", "blue")
scatter(dapc1, ratio.pca=0.3, bg="white", pch=20, cell=0 ,cstar=0, col=myCol,
        solid=.4, cex=3, clab=0,mstree=TRUE, scree.da=FALSE,
        posi.pca="bottomright", leg=TRUE, txt.leg=paste("Cluster",1:6))
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
# MF323 was not among the loci deviating a bit from HWE

round(head(dapc1$posterior),3)
summary(dapc1)
assignplot(dapc1, subset=1:50)
compoplot(dapc1, subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:6), lab="",ncol=2, xlab="individuals")

# a-score analysis
dapc2 <- dapc(MF_genind, n.da=100, n.pca=10)

## DAPC by study area
# SWLU
grp_SWLU <- find.clusters(MF_genind[MF_genind@pop=="SW"|MF_genind@pop=="LU"], max.n.clust=40)
1000
2

dapc_SWLU <- dapc(MF_genind[MF_genind@pop=="SW"|MF_genind@pop=="LU"], grp_SWLU$grp)
60
2

plot(as.data.frame(MF_genind[MF_genind@pop=="SW"|MF_genind@pop=="LU"]@other$xy), col=dapc_SWLU$grp)

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

round(head(dapc1$posterior),3)
summary(dapc1)
assignplot(dapc1, subset=1:50)
compoplot(dapc1, subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:6), lab="",ncol=2, xlab="individuals")

# CO
grp_CO <- find.clusters(MF_genind[MF_genind@pop=="CO"], max.n.clust=40)
1000
2

dapc_CO <- dapc(MF_genind[MF_genind@pop=="CO"], grp_CO$grp)
60
1

plot(as.data.frame(MF_genind[MF_genind@pop=="CO"]@other$xy), col=dapc_CO$grp)

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
contrib <- loadingplot(dapc1$var.contr, axis=2, thres=.07, lab.jitter=1)
# MF323 was not among the loci deviating a bit from HWE

round(head(dapc1$posterior),3)
summary(dapc1)
assignplot(dapc1, subset=1:50)
compoplot(dapc1, subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:6), lab="",ncol=2, xlab="individuals")

################################################################################
############### Spatial analyses

### Classic IBD

empir_geo_dist_MF <- as.matrix(dist(as.data.frame(MF_genind@other$xy)))
empir_geo_dist_MF_SW <- as.matrix(dist(as.data.frame(MF_genind[MF_genind@pop=="SW"]@other$xy)))
empir_geo_dist_MF_LU <- as.matrix(dist(as.data.frame(MF_genind[MF_genind@pop=="LU"]@other$xy)))
empir_geo_dist_MF_CO <- as.matrix(dist(as.data.frame(MF_genind[MF_genind@pop=="CO"]@other$xy)))
empir_geo_dist_MF_SWLU <- as.matrix(dist(as.data.frame(MF_genind[MF_genind@pop=="SW"|MF_genind@pop=="LU"]@other$xy)))

empir_geo_dist_MF2 <-empir_geo_dist_MF
empir_geo_dist_MF_SW2 <-empir_geo_dist_MF_SW
empir_geo_dist_MF_LU2 <-empir_geo_dist_MF_LU
empir_geo_dist_MF_CO2 <-empir_geo_dist_MF_CO
empir_geo_dist_MF_SWLU2 <-empir_geo_dist_MF_SWLU

empir_geo_dist_MF2[empir_geo_dist_MF2==0] <- NA
empir_geo_dist_MF_SW2[empir_geo_dist_MF_SW2==0] <- NA
empir_geo_dist_MF_LU2[empir_geo_dist_MF_LU2==0] <- NA
empir_geo_dist_MF_CO2[empir_geo_dist_MF_CO2==0] <- NA
empir_geo_dist_MF_SWLU2[empir_geo_dist_MF_SWLU2==0] <- NA

MF_genind_nogeo <- MF_genind
MF_genind_nogeo@other <- NULL
MF_genind_nogeo@pop <- NULL

empirLoiselle_EcoGenetics_MF <- eco.kin.loiselle(genind2ecogen(MF_genind_nogeo))
empirLoiselle_EcoGenetics_MF_SW <- eco.kin.loiselle(genind2ecogen(MF_genind_nogeo[MF_genind@pop=="SW"]))
empirLoiselle_EcoGenetics_MF_LU <- eco.kin.loiselle(genind2ecogen(MF_genind_nogeo[MF_genind@pop=="LU"]))
empirLoiselle_EcoGenetics_MF_CO <- eco.kin.loiselle(genind2ecogen(MF_genind_nogeo[MF_genind@pop=="CO"]))
empirLoiselle_EcoGenetics_MF_SWLU <- eco.kin.loiselle(genind2ecogen(MF_genind_nogeo[MF_genind@pop=="SW"|MF_genind@pop=="LU"]))

IBD_MF <- lm(c(as.dist(empirLoiselle_EcoGenetics_MF))~log(c(as.dist(empir_geo_dist_MF2))))
IBD_MF_SW <- lm(c(as.dist(empirLoiselle_EcoGenetics_MF_SW))~log(c(as.dist(empir_geo_dist_MF_SW2))))
IBD_MF_LU <- lm(c(as.dist(empirLoiselle_EcoGenetics_MF_LU))~log(c(as.dist(empir_geo_dist_MF_LU2))))
IBD_MF_CO <- lm(c(as.dist(empirLoiselle_EcoGenetics_MF_CO))~log(c(as.dist(empir_geo_dist_MF_CO2))))
IBD_MF_SWLU <- lm(c(as.dist(empirLoiselle_EcoGenetics_MF_SWLU))~log(c(as.dist(empir_geo_dist_MF_SWLU2))))

summary(IBD_MF)
summary(IBD_MF_SW)
summary(IBD_MF_LU)
summary(IBD_MF_CO)
summary(IBD_MF_SWLU)

mantel.randtest(as.dist(empir_geo_dist_MF), as.dist(1-empirLoiselle_EcoGenetics_MF), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_MF_SW), as.dist(1-empirLoiselle_EcoGenetics_MF_SW), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_MF_LU), as.dist(1-empirLoiselle_EcoGenetics_MF_LU), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_MF_CO), as.dist(1-empirLoiselle_EcoGenetics_MF_CO), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_MF_SWLU), as.dist(1-empirLoiselle_EcoGenetics_MF_SWLU), nrepet = 9999)

plot(mantel.randtest(as.dist(empir_geo_dist_MF_CO), as.dist(1-empirLoiselle_EcoGenetics_MF_CO), nrepet = 9999))
plot(mantel.randtest(as.dist(empir_geo_dist_MF_SWLU), as.dist(1-empirLoiselle_EcoGenetics_MF_SWLU), nrepet = 9999))

plot(log(empir_geo_dist_MF2), empirLoiselle_EcoGenetics_MF)
abline(IBD_MF, col="red")
plot(log(empir_geo_dist_MF_SW2), empirLoiselle_EcoGenetics_MF_SW)
abline(IBD_MF_SW, col="red")
plot(log(empir_geo_dist_MF_LU2), empirLoiselle_EcoGenetics_MF_LU)
abline(IBD_MF_LU, col="red")
plot(log(empir_geo_dist_MF_CO2), empirLoiselle_EcoGenetics_MF_CO)
abline(IBD_MF_CO, col="red")
plot(log(empir_geo_dist_MF_SWLU2), empirLoiselle_EcoGenetics_MF_SWLU)
abline(IBD_MF_SWLU, col="red")

E_MF_SWLU <- as.vector(empirLoiselle_EcoGenetics_MF_SWLU)
E_MF_SWLU[is.na(E_MF_SWLU)] <- 0

dens <- kde2d(as.vector(empir_geo_dist_MF_SWLU), E_MF_SWLU, n=300)
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(empir_geo_dist_MF_SWLU, empirLoiselle_EcoGenetics_MF_SWLU, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
dist_lm <- lm(as.vector(empirLoiselle_EcoGenetics_MF_SWLU) ~ as.vector(empir_geo_dist_MF_SWLU))
abline(dist_lm)
title("Isolation by distance plot")

# ### Monmonnier # DOES NOT WORK FOR M. FLOREA (EITHER STRUCTURE TOO WEAK OR WORKS ONLY WELL WITH POPULATIONS NOT INDIVIDUALS)
# # We need to jitter coords a bit
# geo_jittered <- apply(as.data.frame(MF_genind[MF_genind@pop=="CO"]@other$xy), 2, FUN= function(x) jitter(x,amount = 1000))
# 
# gab_MF_CO <- adegenet::chooseCN(geo_jittered, ask=FALSE, type=2, check.duplicates = FALSE)
# 
# mon1 <- monmonier(geo_jittered, as.dist(1-empirLoiselle_EcoGenetics_MF_CO), gab_MF_CO)
# 
# D_MF_CO <- tab(MF_genind[MF_genind@pop=="CO"], freq = TRUE, NA.method = "mean")
# D_MF_CO <- dist(D_MF_CO)
# pco1 <- dudi.pco(D_MF_CO, scannf=FALSE, nf=2)
# barplot(pco1$eig, main="Eigenvalues")
# D_MF_CO <- dist(pco1$li)
# mon1 <- monmonier(geo_jittered, D_MF_CO, gab_MF_CO)

MF_genind_CO <- MF_genind[MF_genind@pop=="CO"]
geo_jittered <- apply(as.data.frame(MF_genind[MF_genind@pop=="CO"]@other$xy), 2, FUN= function(x) jitter(x,amount = 1000))
MF_genind_CO@other$xy <- geo_jittered
MF_genind_CO@tab <- tab(MF_genind_CO, freq = TRUE, NA.method = "mean")

mySpca <- spca(MF_genind_CO, type=2,ask=FALSE,scannf=FALSE)
barplot(mySpca$eig,main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(1,100)))
barplot(mySpca$eig, main="A variant of the plot\n of sPCA eigenvalues",col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")
screeplot(mySpca)
MF_genind_CO@tab <- tab(MF_genind_CO, freq = TRUE, NA.method = "mean")
myGtest <- global.rtest(MF_genind_CO$tab, mySpca$lw, nperm=999)
myGtest
myLtest <- local.rtest(MF_genind_CO$tab,mySpca$lw, nperm=999)
myLtest

plot(mySpca)
colorplot(mySpca,cex=3,main="colorplot of mySpca, first global score")

library(akima)
x <- other(MF_genind_CO)$xy[,1]
y <- other(MF_genind_CO)$xy[,2]
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


### SWLU
MF_genind_SWLU <- MF_genind[MF_genind@pop=="SW"|MF_genind@pop=="LU"]
geo_jittered <- apply(as.data.frame(MF_genind_SWLU@other$xy), 2, FUN= function(x) jitter(x,amount = 500))
MF_genind_SWLU@other$xy <- geo_jittered
MF_genind_SWLU@tab <- tab(MF_genind_SWLU, freq = TRUE, NA.method = "mean")

mySpca <- spca(MF_genind_SWLU, type=2,ask=FALSE,scannf=FALSE)
barplot(mySpca$eig,main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(1,100)))
barplot(mySpca$eig, main="A variant of the plot\n of sPCA eigenvalues",col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")
screeplot(mySpca)
MF_genind_CO@tab <- tab(MF_genind_SWLU, freq = TRUE, NA.method = "mean")
myGtest <- global.rtest(MF_genind_SWLU$tab, mySpca$lw, nperm=999)
myGtest
myLtest <- local.rtest(MF_genind_SWLU$tab, mySpca$lw, nperm=999)
myLtest

plot(mySpca)
colorplot(mySpca,cex=3,main="colorplot of mySpca, first global score")

library(akima)
x <- other(MF_genind_SWLU)$xy[,1]
y <- other(MF_genind_SWLU)$xy[,2]
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
