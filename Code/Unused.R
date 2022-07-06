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