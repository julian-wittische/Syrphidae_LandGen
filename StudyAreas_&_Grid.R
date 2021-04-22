library(raster)
library(sf)
library(sp)
library(latticeExtra)
library(dplyr)
library(rgeos)
library(rasterVis)
library(reshape2)
library(RColorBrewer)
library(RgoogleMaps)



############################### COLOGNE ###############################

# Load, define CRS, check
cologne <- raster("ascii/cologne_in.asc")
crs(cologne) <- CRS('+init=EPSG:2169')
cologne
nlayers(cologne)
plot(cologne)
res(cologne)
sum(is.na(values(cologne)))
sum(!is.na(values(cologne)))

# Change resolution to 1 kilometer, check
cologne1km <- aggregate(cologne, fact=1000/res(cologne))
res(cologne1km) #checking step
dim(cologne1km)
ncell(cologne1km)
plot(cologne1km)

cell_number_cologne <- cologne1km
cell_number_cologne[!is.na(cell_number_cologne)] <- 1:length(cell_number_cologne[!is.na(cell_number_cologne)])

pdf("FirstMapDraft.pdf", width=10, height=10)
plot(cologne, xlim=c(121000, 151000), legend=FALSE)
plot(rasterToPolygons(cologne1km, digits=20), add=TRUE, border='black', lwd=0.001)
text(cell_number_cologne, cex=0.55)
dev.off()

latlonwgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
cologne1km_ll <- projectRaster(cologne1km,
                              crs = latlonwgs84,
                              method = "ngb")

centroids <- as.data.frame(cologne1km_ll, xy=TRUE, na.rm=TRUE, centroids=TRUE)
dim(centroids)
#  
# write.csv(centroids, "centroids.csv")
# gmap <- GetMap.bbox(c(6.7, 7.2), c(50.8, 51.2))
# PlotOnStaticMap(gmap)
# image(gmap)

# cell_number_cologne_ll <- projectRaster(cell_number_cologne, crs = CRS("+proj=longlat +datum=WGS84"))
# plot(cell_number_cologne_ll, colNA = "black")
# 
# spplot(rasterToPolygons(cologne1km))
# sp.layoutlist('sp.pointLabel', pts, label=someAuthors,
# spplot(cologne1km)
# 








# cell_sel <- cologne1km
# cell_sel <- cellFromXY(cell_sel, cbind(c(130000,135000), c(220500, 230000)))
# # cell_sel <- cellFromRowColCombine(cell_sel, row=2:8, col=5:15)
# cell_sel <- rasterFromCells(cologne1km, cell_sel, values=TRUE)
# plot(cell_sel, add=TRUE, col="red")
# 
# colr <- colorRampPalette(brewer.pal(11, 'YlOrRd'))
# 
# myPanel <- function(x, y, z, ...) {
#   panel.levelplot(x,y,z,...)
#   panel.text(x, y,  rtp[cbind(x,y)]) ## use handy matrix indexing
# }
# 
# rtp <- as.data.frame(rasterToPolygons(cologne1km, digits=40)
# 
# levelplot(cologne1km,
#           margin = FALSE,
#           colorkey=NULL,
#           col.regions=colr,
#           panel=myPanel)
# 








latlonwgs84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
cologne1kmll <- projectRaster(cologne1km,
                              crs = latlonwgs84,
                              method = "ngb")





# ############################### LUXEMBOURG CITY ###############################
# 
# # Load, define CRS, check
# lux <- raster("ascii/lux_in.asc")
# crs(lux) <- CRS('+init=EPSG:2169')
# lux
# nlayers(lux)
# plot(lux)
# res(lux)
# 
# # Change resolution to 1 kilometer, check
# lux1km <- aggregate(lux, fact=1000/res(lux))
# res(lux1km) #checking step
# dim(lux1km)
# ncell(lux1km)
# plot(lux1km)
# 
# plot(lux1km, colNA = "black")
# cell_number_lux <- lux1km
# # sum(getValues(cell_number_lux)!="NaN")
# cell_number_lux[which(getValues(cell_number_lux)!="NaN")] <- 1:sum(getValues(cell_number_lux)!="NaN")
# cell_number_lux[which(getValues(cell_number_lux)=="NaN")] <- NA
# plot(cell_number_lux, colNA = "black")
# 
# # # FOR EASE:
# # cell_number_lux <- projectRaster(cell_number_lux, crs = CRS("+proj=longlat +datum=WGS84"))
# # plot(cell_number_lux, colNA = "black")
# 
# teams <- cell_number_lux
# teams[teams < 40] <- as.integer(1)
# teams[teams %in% c(67:73, 84:91, 99:104, 116:127, 133:142, 152:156)] <- as.integer(2)
# teams[teams > 177] <- as.integer(3)
# teams[!teams %in% c(1:3,NA)] <- as.integer(4)
# 
# plot(teams, colNA = "black")
# levels(teams)=data.frame(ID=1:4, code=LETTERS[1:4])
# rasterVis::levelplot(teams)
# text(teams)
