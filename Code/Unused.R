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

SP_genind@loc.n.all
SP_genind@all.names
stepcheck <- data.frame(apply(SP_genind@tab, 2, as.numeric))
(lol <- stepcheck[,grep("Spp051", colnames(stepcheck))])
colSums(lol, na.rm=TRUE)
colSums(lol, na.rm=TRUE)[sort(names(colSums(lol, na.rm=TRUE)))]


SP_df_raw[which(stepcheck$Spp273.279>0),2]

stepcheck 


SP_genind@loc.n.all
colSums(stepcheck, na.rm=TRUE)[11:27]

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

# SP_genind@loc.n.all
# SP_genind@all.names
# stepcheck <- data.frame(apply(SP_genind@tab, 2, as.numeric))
# (lol <- stepcheck[,grep("Spp051", colnames(stepcheck))])
# colSums(lol, na.rm=TRUE)
# colSums(lol, na.rm=TRUE)[sort(names(colSums(lol, na.rm=TRUE)))]
# SP_df_raw[which(stepcheck$Spp273.279>0),2]
# stepcheck 
# SP_genind@loc.n.all
# colSums(stepcheck, na.rm=TRUE)[11:27]
