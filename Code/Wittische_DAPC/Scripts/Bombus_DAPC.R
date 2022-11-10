# This script performs a Discriminant Analysis of Principal Components (DAPC) using the adegenet package.

# Notes:
# Main advantages of this approach is that is does not come with 'typical' population genetic assumptions 
# (i.e. data does not need to be in HWE, LE). Additionally, by minimizing the variation within groups,
# variation among groups is maximized, resulting in high power to detect even subtle and gradual changing population structure.

# To consider:
# Only samples with non-missing data are incorporated here i.e. if a sample has missing data at a particular locus, the entire locus
# is considered missing. This works for diploids and haploids where allelic dosage is unambiguous. 


# This script is based on the official DAPC tutorial for adegenet v.2.1.5 and uses adegenet v.2.1.5

# Input data is dpds without Billed, as too small populations might not represent the intrapopulation variation correctly.

# Code was checked for correctness on 20211220 by MG.


# Packages, data import and sanity checks ------------------------------------------------

# Specifying packages to be loaded
pckg <- c("here", "adegenet", "poppr", "RColorBrewer", "pals", "png", "grid", "gridExtra")
# Loading packages
lapply(pckg, require, character.only = TRUE)


# Importing genetic data
# Genetic data is provided in Genalex format (saved as .csv file) and converted to a genind-class object using the poppr package v.2.9.3
# No spatial location data are imported. The header information fit the matrix information perfectly (i.e. no mismatch between header and data)

setwd(here("Input"))
gendat <- poppr::read.genalex(
  "Final_5_without_Billed.csv",
  ploidy = 2,
  geo = FALSE, # FALSE = no spatial coordinates are imported
  region = FALSE, # FALSE = no region information is imported
  genclone = FALSE, # FALSE = data imported as genind object
  sep = ",",
  recode = FALSE # only needed if ploidy > 2
)

# Get a summary of the imported genetic data and help on how to access certain parts of the data
gendat


# All information has been imported correctly!

# Finding a meaningful cluster size de novo -------------------------------------

# This is an alternative to using populations as priors and might be suited best for analyses of population structure
# Goal: The most meaningful set of parameters is found de novo
# Background: Performing a search for the most likely number of clusters is preferred over
# simply assuming a certain number of clusters. The find.clusters() function performs a
# PCA, then runs a k-means algorithm for increasing values of K and eventually computes
# a BIC values for each K. The statistically most meaningful K is the one that minimizes the BIC
# value

# Note: Apart from saving some computational time, there is no advantage of taking a subset of axes only.
# Hence, the search for the most likely number of K should be performed using all axes; here achieved by setting n.pca to 100.
# Note, the number of randomly chosen starting centroids (default = 10) was increased to solve model convergence issues
options(nwarnings = 1000000) 
setwd(here("Output"))
png(filename = "BICxclusters_denovo.png", width = 7, height = 5, units = "in",  res = 300)
grp_denovo <- find.clusters(gendat, max.n.clust=40, n.pca = 100, stat ="BIC", n.iter = 1000000, n.start = 700) 
dev.off()

# In 52 out of 700 runs, the model did not converge


# As producing this file is very time-consuming, we may store it as rds file 
setwd(here("Input"))
#saveRDS(grp_denovo, file="grp_denovo.rds")

# If needed, import the rds file again
grp_denovo <- readRDS("grp_denovo.rds")

# Choosing the number of clusters interactively is highly recommended
# Write down the number of clusters that minimize BIC
grp_denovo$Kstat
min(grp_denovo$Kstat) # shows the minimal BIC computed

# The number of clusters that minimizes the BIC value is 7.


# Performing DAPC ---------------------------------------------------------

# Important: Unlike for the find.clusters function performed above, not too many PCs
# should be used in DAPC as retaining too many PCs with respect to the numbers of individuals used can result in
# over-fitting of the data and sub-optimal clustering performance (i.e. data is oversplit).

# Bottomn line
# - Aim to retain a few axes without sacrificing too much information
# - If only few clusters (not tens of clusters) are analysed, one can take all discriminants

# Background: The number of discriminant functions to be retained must be specified. If tens of clusters were specified
# it is likely that only the first e.g. five discriminant functions will explain a lot of variation and hence only these may
# be retained while keeping information loss to a minimum; discriminant functions = linear combinations of variables (principal components of PCA)

# 2 DAPCs are performed: One with populations as priors and the other one with a number of clusters identified de novo (see above):

# Running de novo DAPC:
# n.pca was selected after inspecting the variance x # PCs plot
dapc_denovo_no <- dapc(x= gendat, pop= grp_denovo$grp, scale = FALSE, n.pca = 40, n.da = 100)# setting n.da to a high value ensures that all discriminant functions are retained; no = not optimized
# Running DAPC using population priors:
dapc_priors_no <- dapc(x= gendat, pop= gendat$pop, scale = FALSE, n.pca = 40, n.da = 14)# here, we do not retain all discriminants, but the first 14 which explain most of the variation


# Data visualization ------------------------------------------------------

# Use color-blind friendly colors:

pic_dapc_denovo_no <- scatter(dapc_denovo_no, bg="white",pch=19, cstar=0, col=palette("Paired"), scree.pca=TRUE,
        posi.pca="bottomleft", cex = 2, bg.inset = "white", ratio.da = 0.2, ratio.pca = 0.2, clab=0)

pic_dapc_priors_no <- scatter(dapc_priors_no, bg="white",pch=19, cstar=0, col=palette("Paired"), scree.pca=TRUE,
                              posi.pca="bottomleft", cex = 2, bg.inset = "white", ratio.da = 0.2, ratio.pca = 0.2, clab = 1)


# Avoiding over and under-fitting ------------------------------------------

# Using too many PCs is not beneficial as it results in over-fitting of the data = the model is
# then flexible enough to discriminate between any possible clusters, which results in a loss of predictive capacity.

# Solution: using either a-score or cross-validation to optimize the number of PCs retained

# a-score
# Background: A trade-off between power of discrimination and over-fitting can be measured by the a-score, 
# which is simply the difference between the proportion of successful reassignment of
# the analysis (observed discrimination) and values obtained using random groups (random discrimination)

# Interpretation: An a-score close to one is a sign that the DAPC solution is both strongly 
# discriminating and stable, while low values (toward 0 or lower) indicate either weak discrimination or instability of the results.

# For the de novo approach:
setwd(here("Output"))
ascore_dapc_denovo_no <- a.score(dapc_denovo_no, n.sim = 500)
ascore_dapc_denovo_no$mean # Mean a-score across all populations
write.csv(ascore_dapc_denovo_no$mean, file = "Ascore_mean_denovo_not_optimized.csv")

# For the population prior approach:
ascore_dapc_priors_no <- a.score(dapc_priors_no, n.sim = 500)
ascore_dapc_priors_no$mean # Mean a-score across all populations
write.csv(ascore_dapc_priors_no$mean, file = "Ascore_mean_priors_not_optimized.csv")

# Optimize the a-scores by selecting the optimal number of PCs to retain
# Optimal values are provided in $best - also, provide these graphs in the manuscript
# Because we retained only 40 PCs in the initial DAPC, there is no need to increase n.pca.
png(filename = "DAPC_denovo_opt_PCs.png", width = 6, height = 5, units = "in", res = 300)
score_dapc_denovo_op <- optim.a.score(dapc_denovo_no, n.pca = 1:50, n.sim=1000, smart = FALSE) 
dev.off()

png(filename = "DAPC_priors_opt_PCs.png", width = 6, height = 5, units = "in", res = 300)
score_dapc_priors_op <- optim.a.score(dapc_priors_no, n.pca = 1:50, n.sim=1000, smart = FALSE)
dev.off()


# Perform DAPCs with the optimal number of PCs
dapc_denovo_op <- dapc(x= gendat, pop= grp_denovo$grp, scale = FALSE, n.pca = score_dapc_denovo_op$best, n.da = 100)
dapc_priors_op <- dapc(x= gendat, pop= gendat$pop, scale = FALSE, n.pca = score_dapc_priors_op$best, n.da = 100)

# Save optimized DAPCs to rds files
saveRDS(dapc_denovo_op, file = "dapc_denovo_op.rds")
saveRDS(dapc_priors_op, file = "dapc_priors_op.rds")

# Import them here if needed
dapc_denovo_op <- readRDS(file = "dapc_denovo_op.rds")
dapc_priors_op <- readRDS(file = "dapc_priors_op.rds")

# Display optimized a-scores
ascore_dapc_denovo_op <- a.score(dapc_denovo_op)
ascore_dapc_denovo_op$mean
write.csv(ascore_dapc_denovo_op$mean, file = "Ascore_mean_denovo_optimized.csv")
ascore_dapc_priors_op <- a.score(dapc_priors_op)
ascore_dapc_priors_op$mean
write.csv(ascore_dapc_priors_op$mean, file = "Ascore_mean_priors_optimized.csv")


# Visualizing optimized DAPCs ---------------------------------------------
setwd(here("Output"))
png(filename = "DAPC_optimized_denovo.png", width = 8, height = 6, units = "in", res = 300)
pic_dapc_denovo_op <- scatter(dapc_denovo_op, bg="white", pch =19, cstar=0, col=palette("Dark2"), scree.pca=TRUE,
                              posi.pca="bottomleft", cex = 2, bg.inset = "white", ratio.da = 0.20, ratio.pca = 0.2, clab=1, cellipse = 2.5) #cellipse = 2.5 corresponds to a an alpha threshold of 95% 
dev.off()

# Needs a colour palette with 21 colours; the pals package will be used
png(filename = "DAPC_optimized_priors.png", width = 8, height = 6, units = "in", res = 300)
pic_dapc_priors_op <- scatter(dapc_priors_op, bg="white",pch=19, cstar=0, col=unname(cols25()), scree.pca=TRUE,
                              posi.pca="bottomleft", cex = 2, bg.inset = "white", ratio.da = 0.20, ratio.pca = 0.2, clab=FALSE, cellipse= 2.5)
dev.off()


# Arranging pics
setwd(here("Output"))
pc_sel1 <- readPNG("DAPC_priors_opt_PCs.png")
pc_sel2 <- readPNG("DAPC_denovo_opt_PCs.png")


png("DAPC_PC_selection_graphs.png", height = 7, width = 4, units = "in", res= 300)
PCsel <- grid.arrange(rasterGrob(pc_sel1),rasterGrob(pc_sel2), ncol =1)
dev.off()


# Export of miscellaneous results -------------------------------------------------------
# Export a list showing what individuals were assigned to what cluster
setwd(here("Output"))
denovo_cluster_assign <- as.data.frame(grp_denovo$grp)
write.csv(denovo_cluster_assign, "Denovo_cluster_assignment_individuals.csv")


# Cross-validation of PCs retained (not done here)----------------------------------------
# Important: The number of PCs retained can have substantial impact on the results
# Using cross-validation can help finding the perfect number of PCs to retain and also 

# PROVIDES INFORMATION ON THE PREDICTIVE POWER OF THE DAPC!!!

# Approach: DAPC is carried out on the training set with variable numbers of PCs retained, and the
# degree to which the analysis is able to accurately predict the group membership of excluded
# individuals (those in the validation set) is used to identify the optimal number of PCs to
# retain. At each level of PC retention, the sampling and DAPC procedures are repeated n.rep times

# mat <- tab(gendat, NA.method="mean")
# grp <- pop(gendat)
# 
# png(filename = "Crossval_dpds_without_billed.png", width = 6, height = 5, units = "in", res = 300)
# dapc_xval <- xvalDapc(mat, grp, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 200, xval.plot = TRUE) # Setting n.pca.max to a high value will ensure that all PCS are used in the cross-validation; ; n.pca = NULL: the axes to be retained for the cross validation will be determined automatically
# dev.off()
# 
# 
# # Print results, including the recommended number of PCs_
# # Advice: Select the number of PCs that show the lowest RMSE!
# dapc_xval[2:6]
# print(dapc_xval)