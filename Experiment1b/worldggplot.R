library(rgdal)
library(maptools)
library(ggplot2)
library(broom)
lands <- readOGR(dsn = "C:/ZSwork/glbm/Experiment2/ne_110m_land", layer = "ne_110m_land")
lands_data <- tidy(lands) 
library(rgeos)

## The ogriginal data is in [-180, 180] and we want to change it to [0, 360]
## First split the polygons cross the zero degree
## identify which polygons
lands@polygons
spLand <- SpatialPolygons(lands@polygons)
plot(spLand, col = "grey")
summary(spLand)
