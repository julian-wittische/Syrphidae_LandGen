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