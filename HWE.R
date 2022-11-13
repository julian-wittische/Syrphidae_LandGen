################################################################################
######### Julian Wittische - Hoverfly landscape genetics #########
################################################################################

#-------------------------------------------------------------------------------
############################### Loading packages ###############################
#-------------------------------------------------------------------------------

# Source dta loading script (SP_genind is the object with all 14 loci)

source("Code/GenDataPrep.R")

### Hardy-Weinberg equilibrium within study areas

HW_SP_LU <- hw.test(SP_genind[SP_genind@pop=="LU"], B=1000)
HW_SP_CO <- hw.test(SP_genind[SP_genind@pop=="CO"], B=1000)
names(which(HW_SP_LU[,4] < 0.05 & HW_SP_CO[,4] < 0.05))
# RESULT: Spp142, Spp051, Spp108, Spp141, Spp360  are somewhat in disequilibrium
# but let us keep them because Alain has subsampled 30 individuals 10 times and they pass (spring 2022)

### Hardy-Weinberg equilibrium within communes within study areas

SP_genind

HW_SP_LU <- hw.test(SP_genind[SP_genind@pop=="LU"], B=1000)
HW_SP_CO <- hw.test(SP_genind[SP_genind@pop=="CO"], B=1000)
names(which(HW_SP_LU[,4] < 0.05 & HW_SP_CO[,4] < 0.05))
# RESULT: Spp142, Spp051, Spp108, Spp141, Spp360  are somewhat in disequilibrium
# but let us keep them because Alain has subsampled 30 individuals 10 times and they pass
