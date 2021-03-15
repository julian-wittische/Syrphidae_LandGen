library(raster)
library(sf)
library(sp)
library(latticeExtra)
library(dplyr)
library(rgeos)

#cologne <- stack("C:/Users/Utilisateur/Desktop/Projects/Syrphidae_LandGen/ascii/cologne_in.asc")

cologne <- raster("C:/Users/Utilisateur/Desktop/Projects/Syrphidae_LandGen/ascii/cologne_in.asc")
crs(cologne) <- CRS('+init=EPSG:4181')
nlayers(cologne)

plot(cologne)

# put at 1km res
cologne1km <- aggregate(cologne, fact=100)
res(cologne1km)
dim(cologne1km)

#checking step
plot(cologne1km, colNA = "black")
cell_number <- cologne1km
sum(getValues(cell_number)!="NaN")



cell_number[which(getValues(cell_number)!="NaN")] <- 1:sum(getValues(cell_number)!="NaN")

plot(cologne1km, colNA = "black")
text(cell_number, cex=0.4)
# 
# cell_sel <- cologne1km
# cell_sel <- cellFromXY(cell_sel, cbind(c(130000,135000), c(220500, 230000)))
# # cell_sel <- cellFromRowColCombine(cell_sel, row=2:8, col=5:15)
# cell_sel <- rasterFromCells(cologne1km, cell_sel, values=TRUE)
# plot(cell_sel, add=TRUE, col="red")

