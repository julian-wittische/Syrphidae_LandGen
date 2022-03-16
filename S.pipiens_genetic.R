################################################################################
########################### Julian Wittische - MNHNL ###########################
################################################################################

################################################################################
############### Loading the INCOMPLETE S. pipiens microsatellite dataset

library(adegenet)
library(pegas)
library(hierfstat)
library(sp)
library(EcoGenetics)
library(MASS)

# Read
SP_df_raw <- as.data.frame(readxl::read_excel("SPipiens_raw.xlsx",
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

################################################################################
############### Basic exploration

SP_genind_summary <- summary(SP_genind)

### Alleles and sample sizes
par(mfrow=c(2,2))
plot(SP_genind_summary$n.by.pop, SP_genind_summary$pop.n.all, xlab="Study area sample size",ylab="Number of alleles",main="Alleles numbers and sample sizes",type="n")
text(SP_genind_summary$n.by.pop,SP_genind_summary$pop.n.all,lab=names(SP_genind_summary$n.by.pop))
barplot(SP_genind_summary$loc.n.all, ylab="Number of alleles",main="Number of alleles per locus")
barplot(SP_genind_summary$Hexp-SP_genind_summary$Hobs, main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")
barplot(SP_genind_summary$n.by.pop, main="Sample sizes per population",ylab="Number of genotypes",las=3)

### Allelic richness
barplot(allelicrichness(as.loci(SP_genind)), beside = TRUE)

### Is mean observed H significantly lower than mean expected H ?
par(mfrow = c(1,1))
bartlett.test(list(SP_genind_summary$Hexp,SP_genind_summary$Hobs))
t.test(SP_genind_summary$Hexp,SP_genind_summary$Hobs,pair=T,var.equal=TRUE,alter="greater")
# Yes! (same mean difference as in other data sets I have seen)

### Hardy-Weinberg equilibrium
hw.test(SP_genind, B=1000)
HW_LU <- hw.test(SP_genind[SP_genind@pop=="LU"], B=1000)
HW_CO <- hw.test(SP_genind[SP_genind@pop=="CO"], B=1000)
HW_SW <- hw.test(SP_genind[SP_genind@pop=="SW"], B=1000)

which(HW_LU[,4] < 0.05 & HW_CO[,4] < 0.05 & HW_SW[,4] < 0.05)
# Spp142, Spp051, Spp108, Spp141, Spp360  are somewhat in disequilibrium but let us keep them
# Alain has subsampled 30 individuals 10 times and they pass

################################################################################
############### Nonspatial population structure

## Weir and Cockerham F statistics

wc(SP_genind) # global Fst is pretty low high but Fis is high

# Let's build confidence intervals for that
SP_g2h <- genind2hierfstat(SP_genind)
boot.vc(SP_g2h[1], SP_g2h[-1])$ci

# Let's check per locus
Fperlocus <- Fst(as.loci(SP_genind))
Fperlocus # super high Fis values (0.96!) and very low Fst values
colMeans(Fperlocus)

# Pairwise Fst
Fst <- genet.dist(SP_genind, method = "Nei87")
is.euclid(Fst) #FALSE because of missing values

### PCA
# Let us replace those NA values
sum(is.na(SP_genind$tab))
PCAdf <- tab(SP_genind, freq = TRUE, NA.method = "mean")

# Check
class(PCAdf)
dim(PCAdf)
PCAdf[1:5,1:5]


PCA <- dudi.pca(PCAdf, scale = FALSE, scannf = FALSE, nf = 3)
# Eigenvalues
barplot(PCA$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

# Viz 1 - General PCA plot
s.label(PCA$li)
title("PCA - S. pipiens \naxes 1-2")
add.scatter.eig(PCA$eig[1:20], 3, 1, 2)
# Three groups as originally shown by Alain

# Viz 2 - Study area focused
s.class(PCA$li, pop(SP_genind))
title("PCA - S. pipiens \naxes 1-2")
add.scatter.eig(PCA$eig[1:20], 3, 1, 2)
# The three groups are obviously not geographical

# Viz 3 - Third axis (it does not seem very important)
s.class(PCA$li,pop(SP_genind),xax=1,yax=3,sub="PCA 1-3",csub=2)
title("PCA - S. pipiens \naxes 1-3")
add.scatter.eig(PCA$eig[1:20], nf=3, xax=1, yax=3)

col <- funky(5)
s.class(PCA$li, pop(SP_genind), xax=1, yax=3, col=transp(col,.6), axesell=FALSE, 
        cstar=0, cpoint=3, grid=FALSE)

# Weird vertical pattern let's test temporal

dates <- rownames(SP_genind@tab)
dates <- gsub("_[^_]*$|^[^_]*_", "", dates, perl=T)
library(stringr)
dates <- str_sub(dates, start= -8)
dates <- as.Date(dates, tryFormats = c("%d-%m-%Y"))
class(dates)
dates <- as.numeric(format(dates, "%j"))
hist(dates)

colorplot(PCA$li, PCA$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA - S. pipiens \naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

colorplot(PCA$li[c(1,3)], PCA$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
title("PCA - S. pipiens \naxes 1-3")
abline(v=0,h=0,col="grey", lty=2)

### CA
obj <- genind2genpop(SP_genind)
ca1 <- dudi.coa(tab(obj),scannf=FALSE,nf=3)
barplot(ca1$eig,main="Correspondance Analysis eigenvalues",col=heat.colors(length(ca1$eig)))
s.label(ca1$li, sub="CA 1-2",csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=1,yax=2,posi="bottomright")
# s.label(ca1$li,xax=2,yax=3,lab=popNames(obj),sub="CA 1-3",csub=2)
# add.scatter.eig(ca1$eig,nf=3,xax=2,yax=3,posi="topleft"

### DAPC
# NOTE: we should do this temporally too
grp <- find.clusters(SP_genind, max.n.clust=40)
100
3

dapc1 <- dapc(SP_genind, grp$grp)
50
2

plot(SP_genind@other$xy, col=dapc1$grp)

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
# SP323 was not among the loci deviating a bit from HWE

round(head(dapc1$posterior),3)
summary(dapc1)
assignplot(dapc1, subset=1:50)
compoplot(dapc1, subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:6), lab="",ncol=2, xlab="individuals")

# a-score analysis
dapc2 <- dapc(SP_genind, n.da=100, n.pca=10)

## DAPC by study area
# SWLU
grp_SWLU <- find.clusters(SP_genind[SP_genind@pop=="SW"|SP_genind@pop=="LU"], max.n.clust=40)
1000
2

dapc_SWLU <- dapc(SP_genind[SP_genind@pop=="SW"|SP_genind@pop=="LU"], grp_SWLU$grp)
60
2

plot(as.data.frame(SP_genind[SP_genind@pop=="SW"|SP_genind@pop=="LU"]@other$xy), col=dapc_SWLU$grp)

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
grp_CO <- find.clusters(SP_genind[SP_genind@pop=="CO"], max.n.clust=40)
1000
7

dapc_CO <- dapc(SP_genind[SP_genind@pop=="CO"], grp_CO$grp)
40
6

plot(as.data.frame(SP_genind[SP_genind@pop=="CO"]@other$xy), col=dapc_CO$grp)

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
# SP323 was not among the loci deviating a bit from HWE

round(head(dapc1$posterior),3)
summary(dapc1)
assignplot(dapc1, subset=1:50)
compoplot(dapc1, subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:6), lab="",ncol=2, xlab="individuals")

################################################################################
############### Spatial analyses

### Classic IBD

empir_geo_dist_SP <- as.matrix(dist(as.data.frame(SP_genind@other$xy)))
empir_geo_dist_SP_SW <- as.matrix(dist(as.data.frame(SP_genind[SP_genind@pop=="SW"]@other$xy)))
empir_geo_dist_SP_LU <- as.matrix(dist(as.data.frame(SP_genind[SP_genind@pop=="LU"]@other$xy)))
empir_geo_dist_SP_CO <- as.matrix(dist(as.data.frame(SP_genind[SP_genind@pop=="CO"]@other$xy)))
empir_geo_dist_SP_SWLU <- as.matrix(dist(as.data.frame(SP_genind[SP_genind@pop=="SW"|SP_genind@pop=="LU"]@other$xy)))

empir_geo_dist_SP2 <-empir_geo_dist_SP
empir_geo_dist_SP_SW2 <-empir_geo_dist_SP_SW
empir_geo_dist_SP_LU2 <-empir_geo_dist_SP_LU
empir_geo_dist_SP_CO2 <-empir_geo_dist_SP_CO
empir_geo_dist_SP_SWLU2 <-empir_geo_dist_SP_SWLU

empir_geo_dist_SP2[empir_geo_dist_SP2==0] <- NA
empir_geo_dist_SP_SW2[empir_geo_dist_SP_SW2==0] <- NA
empir_geo_dist_SP_LU2[empir_geo_dist_SP_LU2==0] <- NA
empir_geo_dist_SP_CO2[empir_geo_dist_SP_CO2==0] <- NA
empir_geo_dist_SP_SWLU2[empir_geo_dist_SP_SWLU2==0] <- NA

SP_genind_nogeo <- SP_genind
SP_genind_nogeo@other <- NULL
SP_genind_nogeo@pop <- NULL

empirLoiselle_EcoGenetics_SP <- eco.kin.loiselle(genind2ecogen(SP_genind_nogeo))
empirLoiselle_EcoGenetics_SP_SW <- eco.kin.loiselle(genind2ecogen(SP_genind_nogeo[SP_genind@pop=="SW"]))
empirLoiselle_EcoGenetics_SP_LU <- eco.kin.loiselle(genind2ecogen(SP_genind_nogeo[SP_genind@pop=="LU"]))
empirLoiselle_EcoGenetics_SP_CO <- eco.kin.loiselle(genind2ecogen(SP_genind_nogeo[SP_genind@pop=="CO"]))
empirLoiselle_EcoGenetics_SP_SWLU <- eco.kin.loiselle(genind2ecogen(SP_genind_nogeo[SP_genind@pop=="SW"|SP_genind@pop=="LU"]))

IBD_SP <- lm(c(as.dist(empirLoiselle_EcoGenetics_SP))~log(c(as.dist(empir_geo_dist_SP2))))
IBD_SP_SW <- lm(c(as.dist(empirLoiselle_EcoGenetics_SP_SW))~log(c(as.dist(empir_geo_dist_SP_SW2))))
IBD_SP_LU <- lm(c(as.dist(empirLoiselle_EcoGenetics_SP_LU))~log(c(as.dist(empir_geo_dist_SP_LU2))))
IBD_SP_CO <- lm(c(as.dist(empirLoiselle_EcoGenetics_SP_CO))~log(c(as.dist(empir_geo_dist_SP_CO2))))
IBD_SP_SWLU <- lm(c(as.dist(empirLoiselle_EcoGenetics_SP_SWLU))~log(c(as.dist(empir_geo_dist_SP_SWLU2))))

summary(IBD_SP)
summary(IBD_SP_SW)
summary(IBD_SP_LU)
summary(IBD_SP_CO)
summary(IBD_SP_SWLU)

mantel.randtest(as.dist(empir_geo_dist_SP), as.dist(1-empirLoiselle_EcoGenetics_SP), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_SP_SW), as.dist(1-empirLoiselle_EcoGenetics_SP_SW), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_SP_LU), as.dist(1-empirLoiselle_EcoGenetics_SP_LU), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_SP_CO), as.dist(1-empirLoiselle_EcoGenetics_SP_CO), nrepet = 9999)
mantel.randtest(as.dist(empir_geo_dist_SP_SWLU), as.dist(1-empirLoiselle_EcoGenetics_SP_SWLU), nrepet = 9999)

plot(mantel.randtest(as.dist(empir_geo_dist_SP_CO), as.dist(1-empirLoiselle_EcoGenetics_SP_CO), nrepet = 9999))
plot(mantel.randtest(as.dist(empir_geo_dist_SP_SWLU), as.dist(1-empirLoiselle_EcoGenetics_SP_SWLU), nrepet = 9999))

plot(log(empir_geo_dist_SP2), empirLoiselle_EcoGenetics_SP)
abline(IBD_SP, col="red")
plot(log(empir_geo_dist_SP_SW2), empirLoiselle_EcoGenetics_SP_SW)
abline(IBD_SP_SW, col="red")
plot(log(empir_geo_dist_SP_LU2), empirLoiselle_EcoGenetics_SP_LU)
abline(IBD_SP_LU, col="red")
plot(log(empir_geo_dist_SP_CO2), empirLoiselle_EcoGenetics_SP_CO)
abline(IBD_SP_CO, col="red")
plot(log(empir_geo_dist_SP_SWLU2), empirLoiselle_EcoGenetics_SP_SWLU)
abline(IBD_SP_SWLU, col="red")

E_SP_SWLU <- as.vector(empirLoiselle_EcoGenetics_SP_SWLU)
E_SP_SWLU[is.na(E_SP_SWLU)] <- 0

dens <- kde2d(as.vector(empir_geo_dist_SP_SWLU), E_SP_SWLU, n=300)
myPal <- colorRampPalette(c("white", "blue", "gold", "orange", "red"))
plot(empir_geo_dist_SP_SWLU, empirLoiselle_EcoGenetics_SP_SWLU, pch=20,cex=.5)
image(dens, col=transp(myPal(300),.7), add=TRUE)
dist_lm <- lm(as.vector(empirLoiselle_EcoGenetics_SP_SWLU) ~ as.vector(empir_geo_dist_SP_SWLU))
abline(dist_lm)
title("Isolation by distance plot")

# ### Monmonnier # DOES NOT WORK FOR M. FLOREA (EITHER STRUCTURE TOO WEAK OR WORKS ONLY WELL WITH POPULATIONS NOT INDIVIDUALS)
# # We need to jitter coords a bit
# geo_jittered <- apply(as.data.frame(SP_genind[SP_genind@pop=="CO"]@other$xy), 2, FUN= function(x) jitter(x,amount = 1000))
# 
# gab_SP_CO <- adegenet::chooseCN(geo_jittered, ask=FALSE, type=2, check.duplicates = FALSE)
# 
# mon1 <- monmonier(geo_jittered, as.dist(1-empirLoiselle_EcoGenetics_SP_CO), gab_SP_CO)
# 
# D_SP_CO <- tab(SP_genind[SP_genind@pop=="CO"], freq = TRUE, NA.method = "mean")
# D_SP_CO <- dist(D_SP_CO)
# pco1 <- dudi.pco(D_SP_CO, scannf=FALSE, nf=2)
# barplot(pco1$eig, main="Eigenvalues")
# D_SP_CO <- dist(pco1$li)
# mon1 <- monmonier(geo_jittered, D_SP_CO, gab_SP_CO)

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


### SWLU
SP_genind_SWLU <- SP_genind[SP_genind@pop=="SW"|SP_genind@pop=="LU"]
geo_jittered <- apply(as.data.frame(SP_genind_SWLU@other$xy), 2, FUN= function(x) jitter(x,amount = 500))
SP_genind_SWLU@other$xy <- geo_jittered
SP_genind_SWLU@tab <- tab(SP_genind_SWLU, freq = TRUE, NA.method = "mean")

mySpca <- spca(SP_genind_SWLU, type=2,ask=FALSE,scannf=FALSE)
barplot(mySpca$eig,main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(1,100)))
barplot(mySpca$eig, main="A variant of the plot\n of sPCA eigenvalues",col=spectral(length(mySpca$eig)))
legend("topright", fill=spectral(2),leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")
screeplot(mySpca)
SP_genind_CO@tab <- tab(SP_genind_SWLU, freq = TRUE, NA.method = "mean")
myGtest <- global.rtest(SP_genind_SWLU$tab, mySpca$lw, nperm=999)
myGtest
myLtest <- local.rtest(SP_genind_SWLU$tab, mySpca$lw, nperm=999)
myLtest

plot(mySpca)
colorplot(mySpca,cex=3,main="colorplot of mySpca, first global score")

library(akima)
x <- other(SP_genind_SWLU)$xy[,1]
y <- other(SP_genind_SWLU)$xy[,2]
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