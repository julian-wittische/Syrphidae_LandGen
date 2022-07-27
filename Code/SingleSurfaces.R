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

source("Sample_Data.R")

# Luxembourg SS analysis

r.stack <- stack(LUX_DEM, LUX_IMP, LUX_TCD, LUX_GRA, LUX_WAW)

GA.inputs <- GA.prep(ASCII.dir = r.stack,
                     Results.dir = "./Results/",
                     parallel = 2)

jl.inputs <- jl.prep(n.Pops = dim(SP_LUX@tab)[1],
                     response = lower(propShared(SP_LUX)),
                     CS_Point.File = SP_LUX@other$xy,
                     JULIA_HOME = JULIA_HOME)

jl.optim <- SS_optim(jl.inputs = jl.inputs,
                     GA.inputs = GA.inputs)