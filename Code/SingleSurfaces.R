################################################################################
########## Julian Wittische - July 2022 - Hoverfly landscape genetics ##########
################################################################################
#
# SCRIPT OBJECTIVE:
#
# - running single surface analyses first to check if everything works
#
# - if it takes forever, I will run a 500m test run 
#_______________________________________________________________________________

source("Code/Sample_Data.R")

################################################################################
############################# LUXEMBOURG SS ANALYSIS ###########################
################################################################################

r.stack <- stack(LUX_DEM, LUX_IMP, LUX_TCD, LUX_GRA, LUX_WAW)

############### Syritta pipiens

GA.inputs <- GA.prep(ASCII.dir = r.stack,
                     Results.dir = "./Results/",
                     parallel = 2)

jl.inputs <- jl.prep(n.Pops = dim(SP_LUX@tab)[1],
                     response = 1-lower(propShared(SP_LUX)),
                     CS_Point.File = SP_LUX@other$xy,
                     JULIA_HOME = JULIA_HOME)

jl.optim <- SS_optim(jl.inputs = jl.inputs,
                     GA.inputs = GA.inputs)

# Test1 : SP; 1 - proportion of shared alleles
# A warning message indicates that IBD is positive and likely to create issues
# -99999 issue
# RESULT: 

# Test2 : SP; 1 - Loiselle's kinship
# RESULT: 

# Test3 : SP; 40 PCA axes distance
# RESULT:

############### Myathropa florea

GA.inputs <- GA.prep(ASCII.dir = r.stack,
                     Results.dir = "./Results/",
                     parallel = 2)

jl.inputs <- jl.prep(n.Pops = dim(MF_LUX@tab)[1],
                     response = 1-lower(propShared(MF_LUX)),
                     CS_Point.File = MF_LUX@other$xy,
                     JULIA_HOME = JULIA_HOME)

jl.optim <- SS_optim(jl.inputs = jl.inputs,
                     GA.inputs = GA.inputs)

# Test1 : MF; 1 - proportion of shared alleles
# RESULT: 

# Test2 : MF; 1 - Loiselle's kinship
# RESULT: 

# Test3 : MF; 40 PCA axes distance
# RESULT:

################################################################################
################################ COLOGNE SS ANALYSIS ###########################
################################################################################


r.stack <- stack(COL_DEM, COL_IMP, COL_TCD, COL_GRA, COL_WAW)

############### Syritta pipiens

GA.inputs <- GA.prep(ASCII.dir = r.stack,
                     Results.dir = "./Results/",
                     parallel = 2)

jl.inputs <- jl.prep(n.Pops = dim(SP_COL@tab)[1],
                     response = 1-lower(propShared(SP_COL)),
                     CS_Point.File = SP_COL@other$xy,
                     JULIA_HOME = JULIA_HOME)

jl.optim <- SS_optim(jl.inputs = jl.inputs,
                     GA.inputs = GA.inputs)

# Test1 : SP; 1 - proportion of shared alleles
# RESULT: 

# Test2 : SP; 1 - Loiselle's kinship
# RESULT: 

# Test3 : SP; 40 PCA axes distance
# RESULT: 

############### Myathropa florea

GA.inputs <- GA.prep(ASCII.dir = r.stack,
                     Results.dir = "./Results/",
                     parallel = 2)

jl.inputs <- jl.prep(n.Pops = dim(MF_COL@tab)[1],
                     response = 1-lower(propShared(MF_COL)),
                     CS_Point.File = MF_COL@other$xy,
                     JULIA_HOME = JULIA_HOME)

jl.optim <- SS_optim(jl.inputs = jl.inputs,
                     GA.inputs = GA.inputs)

# Test1 : MF; 1 - proportion of shared alleles
# RESULT: 

# Test2 : MF; 1 - Loiselle's kinship
# RESULT: 

# Test3 : MF; 40 PCA axes distance
# RESULT: 
