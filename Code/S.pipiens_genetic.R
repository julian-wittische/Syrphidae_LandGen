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

# Read
SP_df_raw <- as.data.frame(readxl::read_excel("Data/SPipiens_raw.xlsx",
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
# RESULT: Yes! (mean difference is reasonable - similar ton other data sets I have seen)

### Hardy-Weinberg equilibrium
hw.test(SP_genind, B=1000)
HW_LU <- hw.test(SP_genind[SP_genind@pop=="LU"], B=1000)
HW_CO <- hw.test(SP_genind[SP_genind@pop=="CO"], B=1000)
HW_SW <- hw.test(SP_genind[SP_genind@pop=="SW"], B=1000)

which(HW_LU[,4] < 0.05 & HW_CO[,4] < 0.05 & HW_SW[,4] < 0.05)
# RESULT: Spp142, Spp051, Spp108, Spp141, Spp360  are somewhat in disequilibrium
# but let us keep them because Alain has subsampled 30 individuals 10 times and they pass

################################################################################
############### Nonspatial population structure

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

### PCA
# Let us replace those NA values
sum(is.na(SP_genind$tab))
PCAdf_SP <- tab(SP_genind, freq = TRUE, NA.method = "mean")

# Check
class(PCAdf_SP)
dim(PCAdf_SP)
PCAdf_SP[1:5,1:5]


PCA <- dudi.pca(PCAdf_SP, scale = FALSE, scannf = FALSE, nf = 3)
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

colorplot(PCA$li, PCA$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA based on microsatellite genotypes \naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

colorplot(PCA$li[c(1,3)], PCA$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
title("PCA - S. pipiens \naxes 1-3")
abline(v=0,h=0,col="grey", lty=2)

# Weird vertical pattern let's test temporal

dates <- rownames(SP_genind@tab)
#dates <- gsub("_[^_]*$|^[^_]*_", "", dates, perl=T)
dates <- sub("_[[:alpha:]]", "", dates)
dates <- sub("_[[:alpha:]]", "", dates)
# I have no idea where these errors originated but they were not present in the original metadata
dates[c(1172, 1223, 1225, 1226, 1484, 1764)] <- c("LUX_78_1_26-07-21",
                                               "LUX_220_1_10-07-21",
                                               "LUX_220_1_25-07-21",
                                               "LUX_220_1_28-07-21",
                                               "CO_217_1_07-06-21",
                                               "CO_36_1_31-07-21")
library(stringr)
dates <- str_sub(dates, start= -8)
dates <- as.Date(dates, tryFormats = c("%d-%m-%Y"))
class(dates)
dates <- as.numeric(format(dates, "%j"))
hist(dates)
sum(is.na(dates))
cor(PCA$li$Axis1, dates)
cor(PCA$li$Axis2, dates)

colorplot(PCA$li, data.frame(dates, rep(0,length(dates)), rep(0,length(dates))), transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA - S. pipiens \naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

colorplot(PCA$li, data.frame(dates, dates, dates), transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA - S. pipiens \naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

###### Check linkage disequilibrium

SP_LD_CHECK <- poppr::ia(SP_genind, sample=199)
SP_LD_CHECK
# There is statistically significant association among the markers but the overall correlation is very low.
SP_LD_pair <- poppr::pair.ia(SP_genind)
SP_LD_pair
# RESULT: /!\ HIGH /!\ between Spp141 and Spp051 ! Those two are also somewhat in HWE.

# What if we get rid of Spp141 (LD, HWE, high Fis)
SP_genind_noSpp141 <- SP_genind[loc=c("Spp010", "Spp053", "Spp080", "Spp142",
                                "Spp231", "Spp273", "Spp476", "Spp051",
                                "Spp108", "Spp313", "Spp360", "Spp391",
                                "Spp416"), drop=TRUE]

PCAdf_SP_noSpp141 <- tab(SP_genind_noSpp141, freq = TRUE, NA.method = "mean")

PCA_noSpp141  <- dudi.pca(PCAdf_SP_noSpp141 , scale = FALSE, scannf = FALSE, nf = 20)

colorplot(PCA_noSpp141$li, PCA_noSpp141$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA based on microsatellite genotypes (without Spp141) \naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

colorplot(PCA_noSpp141$li[c(1,3)], PCA_noSpp141$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
title("PCA based on microsatellite genotypes (without Spp141) \naxes 1-3")
abline(v=0,h=0,col="grey", lty=2)

# Further LD check (advised by Florian Privé and James Schnable)
loadingscor <- cor(t(as.matrix(PCA_noSpp141$c1)))
sapply(loadingscor, function(x) which(x > 0.9 & x != 1, arr.ind = TRUE))

which(matrix(loadingscor > 0.9 & loadingscor != 1,
             dim(loadingscor)), arr.ind = T)

loadingscor[which(matrix(loadingscor > 0.9 & loadingscor != 1,
                         dim(loadingscor)), arr.ind = T)]

# There is a strong positive correlation between the alleles 190 and 186 of Spp.108
# There is a strong positive correlation between the alleles 240 and 236 of Spp.313

which(matrix(loadingscor  < -0.9 & loadingscor != -1,
             dim(loadingscor)), arr.ind = T)

loadingscor[which(matrix(loadingscor  < -0.9 & loadingscor != -1,
               dim(loadingscor)), arr.ind = T)]

# There is a strong negative correlation between the alleles 161 and 153 of Spp.010
# WARNING : looking at the loading plot for CO, Spp010 could be an issue
# There is a strong negative correlation between the alleles 160 and 156 of Spp.053
# WARNING : looking at the loading plot for the SWLUCO overall, Spp053 could be an issue
# There is a strong negative correlation between the alleles 135 and 133 of Spp.231
# There is a strong negative correlation between the alleles 236 and 238 of Spp.313
# There is a strong negative correlation between the alleles 240 and 238 of Spp.313

# High correlations between some alleles (same locus) but I am not sure what to do

library(factoextra)


fviz_eig(PCA_noSpp141)
fviz_pca_ind(PCA_noSpp141,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_biplot(PCA_noSpp141, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

SP_genind_noSpp141_360 <- SP_genind[loc=c("Spp010", "Spp053", "Spp142",
                                      "Spp231", "Spp273", "Spp476", "Spp051",
                                      "Spp108", "Spp313", "Spp391",
                                      "Spp416", "Spp080"), drop=TRUE]

PCAdf_SP_noSpp141_360 <- tab(SP_genind_noSpp141_360, freq = TRUE, NA.method = "mean")

PCA_noSpp141_360  <- dudi.pca(PCAdf_SP_noSpp141_360 , scale = FALSE, scannf = FALSE, nf = 20)

# colorplot(PCA_noSpp141_360_080$li, PCA_noSpp141_360_080$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
# title("PCA based on microsatellite genotypes (without Spp141 and Spp360) \naxes 1-2")
# abline(v=0,h=0,col="grey", lty=2)
# 
# SP_genind_noSpp141_360_080 <- SP_genind[loc=c("Spp010", "Spp053", "Spp142",
#                                           "Spp231", "Spp273", "Spp476", "Spp051",
#                                           "Spp108", "Spp313", "Spp391",
#                                           "Spp416"), drop=TRUE]
# 
# PCAdf_SP_noSpp141_360_080 <- tab(SP_genind_noSpp141_360_080, freq = TRUE, NA.method = "mean")
# 
# PCA_noSpp141_360_080  <- dudi.pca(PCAdf_SP_noSpp141_360_080 , scale = FALSE, scannf = FALSE, nf = 20)
# 
# colorplot(PCA_noSpp141_360_080$li, PCA_noSpp141_360_080$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
# title("PCA based on microsatellite genotypes (without Spp141 and Spp360) \naxes 1-2")
# abline(v=0,h=0,col="grey", lty=2)


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

##### PCA per study area (without Spp141)

# CO
PCAdf_SP_CO <- tab(SP_genind_noSpp141[SP_genind_noSpp141@pop=="CO"], freq = TRUE, NA.method = "mean")
PCA_noSpp141_CO  <- dudi.pca(PCAdf_SP_CO , scale = FALSE, scannf = FALSE, nf = 3)

colorplot(PCA_noSpp141_CO$li, PCA_noSpp141_CO$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA based on microsatellite genotypes (without Spp141) \naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

colorplot(as.data.frame(SP_genind_noSpp141[SP_genind_noSpp141@pop=="CO"]@other$xy),
          PCA_noSpp141_CO$li, transp=TRUE, cex=3)

#SWLU
PCAdf_SP_SWLU <- tab(SP_genind_noSpp141[SP_genind_noSpp141@pop=="SW"|SP_genind_noSpp141@pop=="LU"], freq = TRUE, NA.method = "mean")
PCA_noSpp141_SWLU  <- dudi.pca(PCAdf_SP_SWLU , scale = FALSE, scannf = FALSE, nf = 3)

colorplot(PCA_noSpp141_SWLU$li, PCA_noSpp141_SWLU$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA based on microsatellite genotypes (without Spp141) \naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)

### Integrity check
plot(SP_genind_noSpp141@other$xy, col=SP_genind_noSpp141@pop)

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
grp <- find.clusters(SP_genind_noSpp141_360_080, max.n.clust=40)
1000
3

dapc1 <- dapc(SP_genind_noSpp141_360_080, grp$grp)
50
2

plot(SP_genind_noSpp141_360_080@other$xy, col=dapc1$grp)

scatter(dapc1)

myCol <- c("darkblue", "purple", "green", "orange", "red", "blue")
scatter(dapc1, ratio.pca=0.3, bg="white", pch=20, cell=0 ,cstar=0, col=myCol,
        solid=.4, cex=3, clab=0,mstree=TRUE, scree.da=FALSE,
        posi.pca="bottomright", leg=TRUE, txt.leg=paste("Cluster",1:3))
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

empir_geo_dist_SP <- as.matrix(dist(as.data.frame(SP_genind_noSpp141_360_080@other$xy)))
empir_geo_dist_SP_SW <- as.matrix(dist(as.data.frame(SP_genind_noSpp141_360_080[SP_genind_noSpp141_360_080@pop=="SW"]@other$xy)))
empir_geo_dist_SP_LU <- as.matrix(dist(as.data.frame(SP_genind_noSpp141_360_080[SP_genind_noSpp141_360_080@pop=="LU"]@other$xy)))
empir_geo_dist_SP_CO <- as.matrix(dist(as.data.frame(SP_genind_noSpp141_360_080[SP_genind_noSpp141_360_080@pop=="CO"]@other$xy)))
empir_geo_dist_SP_SWLU <- as.matrix(dist(as.data.frame(SP_genind_noSpp141_360_080[SP_genind_noSpp141_360_080@pop=="SW"|SP_genind_noSpp141_360_080@pop=="LU"]@other$xy)))

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

SP_genind_noSpp141_360_080_nogeo <- SP_genind_noSpp141_360_080
SP_genind_noSpp141_360_080_nogeo@other <- NULL
SP_genind_noSpp141_360_080_nogeo@pop <- NULL

empirLoiselle_EcoGenetics_SP <- eco.kin.loiselle(genind2ecogen(SP_genind_noSpp141_360_080_nogeo))
empirLoiselle_EcoGenetics_SP_SW <- eco.kin.loiselle(genind2ecogen(SP_genind_noSpp141_360_080_nogeo[SP_genind_noSpp141_360_080@pop=="SW"]))
empirLoiselle_EcoGenetics_SP_LU <- eco.kin.loiselle(genind2ecogen(SP_genind_noSpp141_360_080_nogeo[SP_genind_noSpp141_360_080@pop=="LU"]))
empirLoiselle_EcoGenetics_SP_CO <- eco.kin.loiselle(genind2ecogen(SP_genind_noSpp141_360_080_nogeo[SP_genind_noSpp141_360_080@pop=="CO"]))
empirLoiselle_EcoGenetics_SP_SWLU <- eco.kin.loiselle(genind2ecogen(SP_genind_noSpp141_360_080_nogeo[SP_genind_noSpp141_360_080@pop=="SW"|SP_genind_noSpp141_360_080@pop=="LU"]))

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