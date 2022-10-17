################################################################################
####### Julian Wittische - September 2022 - Hoverfly landscape genetics ########
################################################################################

library(graph4lg)

genind_to_genepop(SP_genind_noSpp141, output = "Data/SPno141.txt")
pop_gen_index(SP_genind_noSpp141)
a<-gen_graph_indep(SP_genind_noSpp141)
plot(a)
