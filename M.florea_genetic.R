################################################################################
########################### Julian Wittische - MNHNL ###########################
################################################################################

############### Loading the INCOMPLETE M. florea microsatellite dataset

library(adegenet)
library(pegas)
library(hierfstat)
library(sp)
library(EcoGenetics)

# Read
MF_df_raw <- as.data.frame(readxl::read_excel("Mflorea_raw.xlsx", sheet=1, .name_repair="minimal"))

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
row.names(MF_df)[which(duplicated(row.names(MF_df)))] <- "LUX_98_2_08-09-21_A_bis"

# Transform into genind object
MF_genind <- df2genind(MF_df, ncode=3, ploidy=2, sep="/", type="codom", NA.char=" NA/ NA")

# Add pop information (study area)
MF_genind@pop <- as.factor(substr(row.names(MF_df), 1, 2))

# Add coordinates and reproject them
MF_geo_sp <- SpatialPoints(MF_geo, CRS(SRS_string = "EPSG:4326"))
MF_geo_sp <-spTransform(MF_geo_sp, CRS(SRS_string = "EPSG:3035"))
MF_genind@other$xy <- MF_geo_sp

############### Basic exploration
  
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

### Population structure

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

### Multivariate analyses

sum(is.na(MF_genind$tab))

# Let us replace those NA values
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
title("PCA - M.florea \naxes 1-2")
add.scatter.eig(PCA$eig[1:20], 3, 1, 2)
# Three groups as originally shown by Alain

# Viz 2 - Study area focused
s.class(PCA$li, pop(MF_genind))
title("PCA - M.florea \naxes 1-2")
add.scatter.eig(PCA$eig[1:20], 3, 1, 2)
# The three groups are obviously not geographical

# Viz 3 - Third axis (it does not seem very important)
s.class(PCA$li,pop(MF_genind),xax=1,yax=3,sub="PCA 1-3",csub=2)
title("PCA - M.florea \naxes 1-3")
add.scatter.eig(PCA$eig[1:20], nf=3, xax=1, yax=3)
