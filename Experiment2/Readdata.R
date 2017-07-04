## Read in and plot the ssh altimetry data 
library(ncdf4)
alt_nc <- nc_open("glbm/Experiment2/trend_SSH_CCI_200501_201512.nc")
print(alt_nc)
lon <- ncvar_get(alt_nc, "lon")
lat <- ncvar_get(alt_nc, "lat")
trend_ssh <- ncvar_get(alt_nc, "trend_ssh") #note that there re NAs for land datat
err_ssh <- ncvar_get(alt_nc, "err")

alt_data <- data.frame(trend_ssh = as.numeric(trend_ssh), err_ssh = as.numeric(err_ssh))
alt_data$lon <- rep(lon, 180) 
alt_data$lat <- rep(lat, each = 360)

library(ggplot2)

## Create the base of word map
world_map <- map_data("world2")
p <- ggplot() + coord_fixed() + xlab("") + ylab("") 
#Add map to base plot
base_world <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                               colour="black", fill="black")



map_alt <- ggplot(data=alt_data) + geom_raster(aes(x = lon, y = lat, fill = trend_ssh)) + 
  coord_fixed() + xlab("Longitude") + ylab("Latitude") + 
  scale_x_continuous(limits=c(0,360),  expand = c(0, 0)) + scale_y_continuous(limits=c(-90,90),  expand = c(0, 0)) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "orange", midpoint = 0, name = "",
                       limit = c(-20, 20),
                       guide = guide_colorbar(barwidth = 2, barheight = 10, label.position = "right")) 

map_alt2 <- map_alt + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                   colour="black", fill = NA, alpha = 0.3) + ggtitle("Altimetry trend (mm/year)") 

beauty <- 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "white", colour = 'white'), 
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30), 
        axis.text = element_text(size = 30),
        axis.title = element_text(size = 30),
        axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5, size = 30),
        panel.border = element_blank())
map_altf <- map_alt2 + beauty
png("glbm/Experiment2/Altrend_map.png", width = 1200, height = 600, pointsize = 20)
print(map_altf)
dev.off()


map_altE <- ggplot(data=alt_data) + geom_raster(aes(x = lon, y = lat, fill = err_ssh)) + 
  coord_fixed() + xlab("Longitude") + ylab("Latitude") + 
  scale_x_continuous(limits=c(0,360),  expand = c(0, 0)) + scale_y_continuous(limits=c(-90,90),  expand = c(0, 0)) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "orange", midpoint = 2, name = "",
                       limit = c(0,4),
                       guide = guide_colorbar(barwidth = 2, barheight = 10, label.position = "right")) 

map_altE2 <- map_altE + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                   colour="black", fill = NA, alpha = 0.3) + ggtitle("Altimetry error (mm/year)") +beauty


png("glbm/Experiment2/AltErr_map.png", width = 1200, height = 600, pointsize = 20)
print(map_altE2)
dev.off()


### Readin and plot GRACE data
