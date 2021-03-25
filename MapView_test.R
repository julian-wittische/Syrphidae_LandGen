library(mapview)

mapview(teams,
        method = "ngb",
        na.color = rgb(0, 0, 255, max = 255, alpha = 0),
        query.type = "click",
        trim = TRUE, 
        legend = FALSE)
# 
# 
# 
# grid <- SpatialPixels(SpatialPoints(coordinates(teams)[!is.na(values(teams)),]))
# 
# crs(grid) <- CRS('+init=EPSG:2169')
# 
# plot(grid)
# 
# mapview(grid, method = "ngb",
#         na.color = rgb(0, 0, 255, max = 255, alpha = 0),
#         query.type = "click",
#         trim = TRUE, 
#         legend = FALSE)
# 
# g <- as(teams, 'SpatialGridDataFrame')
# 
# mapview(g)
# g <- spTransform(g, CRS("+proj=longlat +datum=WGS84"))
# 
# teams_latlong <- projectRaster(teams, crs = CRS("+proj=longlat +datum=WGS84"))
# grid_latlong <- as(teams_latlong, 'SpatialGridDataFrame')
# 
# mapview(grid_latlong, method = "ngb",
#         na.color = rgb(0, 0, 255, max = 255, alpha = 0),
#         query.type = "click",
#         trim = TRUE, 
#         legend = FALSE, grid = TRUE)
