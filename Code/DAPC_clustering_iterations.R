################################################################################
########## Julian Wittische - July 2022 - Hoverfly landscape genetics ##########
################################################################################
#
# SCRIPT OBJECTIVE:
#
# - check 
#_______________________________________________________________________________

library(adegenet)
library(poppr)

dapc_SP_prior <- dapc(SP, SP@pop)
50
2
optim.a.score(dapc_SWLU, n=100, smart=TRUE)

set.seed(1)
SP_DAPC_xval_first <- xvalDapc(tab(SP, NA.method = "mean"), pop(SP))

set.seed(1)
SP_DAPC_xval_reps <- xvalDapc(tab(SP, NA.method = "mean"), pop(SP),
                              n.pca = 90:110, n.rep = 100,
                              parallel = "snow", ncpus = 4L)

grp_SP <- find.clusters(SP, max.n.clust=40)
1000
10

dapc_SP <- dapc(SP, grp$grp)
50
9

plot(as.data.frame(SP@other$xy), col=dapc$grp, pch=19)

scatter(dapc_SWLU)


contrib <- loadingplot(dapc_SWLU$var.contr, axis=2, thres=.07, lab.jitter=1)

round(head(dapc_SWLU$posterior),3)
summary(dapc_SWLU)
assignplot(dapc_SWLU, subset=1:50)
compoplot(dapc_SWLU, subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:6), lab="",ncol=2, xlab="individuals")


set.seed(1)
set.seed(2)
set.seed(3)

grp_SWLU <- find.clusters(SP_LUX, max.n.clust=40)
1000
10

dapc_SWLU <- dapc(SP_LUX, grp_SWLU$grp)
50
9

plot(as.data.frame(SP_LUX@other$xy), col=dapc_SWLU$grp, pch=19)

scatter(dapc_SWLU)

optim.a.score(dapc_SWLU, n=100, smart=TRUE)

contrib <- loadingplot(dapc_SWLU$var.contr, axis=2, thres=.07, lab.jitter=1)

round(head(dapc_SWLU$posterior),3)
summary(dapc_SWLU)
assignplot(dapc_SWLU, subset=1:50)
compoplot(dapc_SWLU, subset=1:50, posi="bottomright",txt.leg=paste("Cluster", 1:6), lab="",ncol=2, xlab="individuals")
