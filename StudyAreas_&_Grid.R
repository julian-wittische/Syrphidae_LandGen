library(raster)
library(sf)
library(sp)
library(latticeExtra)
library(dplyr)
library(rgeos)
library(rasterVis)

############################### LUXEMBOURG CITY ###############################

# Load, define CRS, check
lux <- raster("ascii/lux_in.asc")
crs(lux) <- CRS('+init=EPSG:2169')
lux
nlayers(lux)
plot(lux)
res(lux)

# Change resolution to 1 kilometer, check
lux1km <- aggregate(lux, fact=1000/res(lux))
res(lux1km) #checking step
dim(lux1km)
ncell(lux1km)
plot(lux1km)

plot(lux1km, colNA = "black")
cell_number_lux <- lux1km
# sum(getValues(cell_number_lux)!="NaN")
cell_number_lux[which(getValues(cell_number_lux)!="NaN")] <- 1:sum(getValues(cell_number_lux)!="NaN")
cell_number_lux[which(getValues(cell_number_lux)=="NaN")] <- NA
plot(cell_number_lux, colNA = "black")

# # FOR EASE:
# cell_number_lux <- projectRaster(cell_number_lux, crs = CRS("+proj=longlat +datum=WGS84"))
# plot(cell_number_lux, colNA = "black")

teams <- cell_number_lux
teams[teams < 40] <- as.integer(1)
teams[teams %in% c(67:73, 84:91, 99:104, 116:127, 133:142, 152:156)] <- as.integer(2)
teams[teams > 177] <- as.integer(3)
teams[!teams %in% c(1:3,NA)] <- as.integer(4)

plot(teams, colNA = "black")
levels(teams)=data.frame(ID=1:4, code=LETTERS[1:4])
rasterVis::levelplot(teams)
text(teams)







# cell_sel <- cologne1km
# cell_sel <- cellFromXY(cell_sel, cbind(c(130000,135000), c(220500, 230000)))
# # cell_sel <- cellFromRowColCombine(cell_sel, row=2:8, col=5:15)
# cell_sel <- rasterFromCells(cologne1km, cell_sel, values=TRUE)
# plot(cell_sel, add=TRUE, col="red")

