# Make ArcGIS shapefiles out of moa (These can be read by R but not by ArcGIS for some reason!)

library(shapefiles)
library(maptools)


SH_Line_to_poly <- function (path,decimate=1) {
  # load shapefiles
  SH_line <- readShapeSpatial(path)
  SH_line_fort <- fortify(SH_line)
  SH_line_fort$long <- SH_line_fort$long
  SH_line_fort$lat <- SH_line_fort$lat
  SH_line_fort <- SH_line_fort[seq(1,nrow(SH_line_fort),decimate),]
  
  
  dd <- data.frame(Id=as.numeric(SH_line_fort$id),
                   X= SH_line_fort$long,
                   Y= SH_line_fort$lat)
  
  ddTable <- data.frame(Id=0,
                        Name=unique(SH_line_fort$id))
  ddShapefile <- convert.to.shapefile(dd, ddTable, "Id", 5) # Make polygon
  return(ddShapefile)
  
}

ddShapefile1 <- SH_Line_to_poly("./AIS_Data/Coastline/moa_groundingline.shp",decimate=1)
write.shapefile(ddShapefile1, "./AIS_Data/ArcGIS_shapefiles/moa_grounding_poly", arcgis=T)

ddShapefile2 <- SH_Line_to_poly("./AIS_Data/Coastline/moa_coastline.shp",decimate=1)
write.shapefile(ddShapefile2, "./AIS_Data/ArcGIS_shapefiles/moa_coast_poly", arcgis=T)

Islands <- read.table("./AIS_Data/Coastline/all_islands.txt",col.names=c("id","x","y"))
dd <- data.frame(Id=Islands$id,X= Islands$x*1000, Y=Islands$y*1000)
ddTable <- data.frame(Id=unique(Islands$id),  Name=unique(Islands$id))
ddShapefile3 <- convert.to.shapefile(dd, ddTable, "Id", 5) # Make polygon
write.shapefile(ddShapefile3, "./AIS_Data/ArcGIS_shapefiles/islands_poly", arcgis=T)






