################################################################################
########################### Julian Wittische - MNHNL ###########################
################################################################################

############### Loading the INCOMPLETE M. florea microsatellite dataset

library(adegenet)
library(pegas)
library(hierfstat)

# Read
MF_df_raw <- as.data.frame(readxl::read_excel("dataminus22.xlsx", sheet=1))

# Set first column as row names and remove
MF_df <- MF_df_raw[,-1]

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
row.names(MF_df) <- as.vector(MF_df_raw[,1])

# Transform into genind object
MF_genind <- df2genind(MF_df, ncode=3, ploidy=2, sep="/", type="codom", NA.char=" NA/ NA")

# Add pop information (study area)
MF_genind@pop <- as.factor(substr(row.names(MF_df), 1, 2))

# Add coordinates
meta <- read.csv("NH/Fieldwork_metadata - SW.csv")
meta$w_ <- toupper(meta$w_)
meta$w_ <- gsub("LOGNE", "", meta$w_)
meta$w_ <- gsub("XEMBOURG CITY", "", meta$w_)


#MF_genind@other$xy <- meta$Study.area[row.names(MF_genind@tab) %in% meta$w_ ]

############### Basic exploration
  
MF_genind_summary <- summary(MF_genind)

### Alleles and sample sizes
par(mfrow=c(2,2))
plot(MF_genind_summary$n.by.pop, MF_genind_summary$pop.n.all, xlab="Study area sample size",ylab="Number of alleles",main="Alleles numbers and sample sizes",type="n")
text(MF_genind_summary$n.by.pop,MF_genind_summary$pop.n.all,lab=names(MF_genind_summary$n.by.pop))
barplot(MF_genind_summary$loc.n.all, ylab="Number of alleles",main="Number of alleles per locus")
barplot(MF_genind_summary$Hexp-MF_genind_summary$Hobs, main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")
barplot(MF_genind_summary$n.by.pop, main="Sample sizes per population",ylab="Number of genotypes",las=3)

### Is mean observed H significantly lower than mean expected H ?
par(mfrow = c(1,1))
bartlett.test(list(MF_genind_summary$Hexp,MF_genind_summary$Hobs))
t.test(MF_genind_summary$Hexp,MF_genind_summary$Hobs,pair=T,var.equal=TRUE,alter="greater")
# Yes! (same mean difference as in other data sets I have seen)

### Hardy-Weinberg equilibrium
hw.test(MF_genind, B=1000)

### Classic IBD

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
