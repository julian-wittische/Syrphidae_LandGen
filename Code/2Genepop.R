################################################################################
####### Julian Wittische - September 2022 - Hoverfly landscape genetics ########
################################################################################

writeGenPop(SP, "Data/SPno141.gen",
            comment="Coordinates are in the file SPgeo.txt")

write.table(as.data.frame(SP@other$xy), "Data/SPgeo.txt",
            col.names=FALSE, row.names = FALSE)

writeGenPop(MF, "Data/MF.gen",
            comment="Coordinates are in the file MFgeo.txt")

write.table(as.data.frame(MF@other$xy), "Data/MFgeo.txt",
            col.names=FALSE, row.names = FALSE)
