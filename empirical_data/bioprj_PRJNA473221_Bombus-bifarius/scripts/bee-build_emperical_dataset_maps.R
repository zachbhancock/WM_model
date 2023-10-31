#idea: make map of species range and sampled genetic points - SPECIFIC FOR BEE

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(sp)
library(marmap)
library(stringr)
library(cowplot)
library(ggalt)
library(ggforce)
library(legendMap)
library(rnaturalearthdata)
library(rnaturalearth)
library(smoothr)
library(concaveman)
#devtools::install_github("3wen/legendMap") # to put N arrow and scale bar on map



rm(list=ls())
gc()



#get a world map
world <- ne_countries(scale = "medium", returnclass = "sf")


#load (already fetched) bathymetric maps of world (one -180 to 180 longitude and another 0 - 360 longitude)
#bathymap <- getNOAA.bathy(lat1 = -90, lon1 =-180, lat2 = 90, lon2 = 180, resolution = 4, keep = FALSE, antimeridian = FALSE, path = NULL)
#save(bathymap, file = "/Users/rachel/Desktop/DivDiv/mapping_GBIF_and_genetic_sample_points/world_NOAA_bathy_res4.Robj")
#load bathymap of world
load("empirical_data/world_NOAA_bathy_res4.Robj")
bathmap.simple <- bathymap
bathmap.simple[bathmap.simple > 0] <- 5000
bathmap.simple[bathmap.simple <= 0] <- -5000

#get GBIF range points and genetic sample points
genetic <- read.delim(paste0("empirical_data/bioprj_PRJNA473221_Bombus-bifarius/lat_long_table-bioprj_PRJNA473221_Bombus-bifarius.txt")) %>% 
  dplyr::rename("y"="lat", "x"="long")
gbif <- read.csv(paste0("empirical_data/bioprj_PRJNA473221_Bombus-bifarius/GBIF_Bombus-bifarius_occurrences.csv")) %>% 
  dplyr::rename("y"="decimalLatitude", "x"="decimalLongitude")
genclusts <- read.csv("empirical_data/bioprj_PRJNA473221_Bombus-bifarius/genetic_cluster_for_sample_map.csv")


#quick look at raw point data
ggplot() + 
  geom_sf(data = world) +
  geom_point(data = gbif, aes(x=x, y=y), colour = "black") +
  geom_point(data = genetic, aes(x=x, y=y), colour = "red") +
  scale_y_continuous(limits = c(21, 70), breaks = c(seq(-70,70,10))) +
  scale_x_continuous(limits = c(-170, -80), breaks = c(seq(-170,170,10)))

#build convex hull (can set length of segments to include to make it smoother)
range <- gbif %>% dplyr::select(x,y) %>% as.matrix()
range <- concaveman::concaveman(points=range, concavity=3.25, length_threshold=9) %>% 
  as.data.frame() %>% 
  dplyr::rename("x"="V1", "y"="V2") %>% as.matrix()
#smooth shape out even more
range <- smoothr::smooth_ksmooth(range, wrap = TRUE, smoothness = 2.5) %>% 
  as.data.frame() %>% dplyr::rename("x"="V1", "y"="V2")
#how does it look?
ggplot() + 
  geom_sf(data = world) +
  geom_point(data = gbif, aes(x=x, y=y), colour = "red") +
  geom_polygon(data = range, aes(x=x, y=y), fill = "transparent", color = "black") +
  scale_y_continuous(limits = c(21, 70), breaks = c(seq(-70,70,10))) +
  scale_x_continuous(limits = c(-170, -80), breaks = c(seq(-170,170,10)))

#get zoomed bathymap
maxlat <- max(gbif$y, genetic$y, range$y) + 3
minlat <- min(gbif$y, genetic$y, range$y) - 7
maxlon <- max(gbif$x, genetic$x, range$x) + 3
minlon <- min(gbif$x, genetic$x, range$x) - 3
if (maxlat > 90) {maxlat = 90}
if (minlat < -90) {minlat = -90}
if (maxlon > 180) {maxlon = 180}
if (minlon < -180) {minlon = -180}

bathymapZOOMED <- getNOAA.bathy(lat1 = minlat, lon1 = minlon, lat2 = maxlat, lon2 = maxlon, resolution = 4, keep = FALSE, antimeridian = FALSE, path = NULL)

bathymapZOOMED.simple <- bathymapZOOMED
bathymapZOOMED.simple[bathymapZOOMED.simple > 0] <- 5000
bathymapZOOMED.simple[bathymapZOOMED.simple <= 0] <- -5000



#world with box for inset -----
inset.plot <- ggplot() + 
  #build map
  geom_raster(data = bathmap.simple, aes(x=x, y=y, fill=z)) +
  geom_contour(data = bathymap, aes(x=x, y=y, z=z),
               breaks=0, #contour for land
               colour="gray70", linewidth=0.15) +
  scale_fill_gradient2(low="white", mid="white", high="gray90", midpoint = 0) +
  scale_y_continuous(expand = c(0, 0), limits = c(21, 70), breaks = c(seq(-70,70,10))) +
  scale_x_continuous(expand = c(0, 0), limits = c(-170, -80), breaks = c(seq(-170,170,10))) +
  #add range polygon
  geom_polygon(data = range, aes(x=x, y=y), fill = "#e6b8f5", color = "transparent", alpha = 0.75) +
  #add zoom box
  geom_rect(aes(xmin = minlon, xmax = maxlon, ymin = minlat, ymax = maxlat), colour = "black", fill = "transparent") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  guides(fill = "none") +
  coord_fixed()
print(inset.plot)






#plot zoomed map again with new clusters
ggplot() + 
  
  #build map
  geom_raster(data = bathymapZOOMED.simple, aes(x=x, y=y, fill=z)) +
  geom_contour(data = bathymapZOOMED, aes(x=x, y=y, z=z),
               breaks=0, #contour for land
               colour="gray70", linewidth=0.3) +
  scale_fill_gradient2(low="transparent", mid="transparent", high="gray90", midpoint = 0) +


  
  #add genetic points
  geom_point(data = out, aes(x=x, y=y, size=newn.samps.perpop), shape = 21, color = "black", fill = "darkorange1") +
  scale_y_continuous(limits = c(35, 50), breaks = c(seq(-70,70,2))) +
  scale_x_continuous(limits = c(-125.5, -116), breaks = c(seq(-170,170,2))) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_bw()

#






  
  #manually make a legend for points
  # scale_color_identity(guide = "legend",
  #                      name = "",
  #                      breaks = c("blue", "darkorange1"),
  #                      labels = c("Occurrence reports (GBIF)", "Genetic samples")) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(30,55,10)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(10,30,10)) +
  legendMap::scale_bar(lon = 31.25, lat = 5.75,
                       distance_lon = 400, distance_lat = 50,
                       distance_legend = 125, dist_unit = "km",
                       arrow_length = 290, legend_size = 2.5,
                       arrow_north_size = 3.25, arrow_distance = 260) +

  ggsn::north(x.min = minlat, x.max = maxlat, y.min = minlon, y.max = maxlon, symbol = 3) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.5),
        legend.position = c(0.79, 0.08),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'mm'),
        axis.text = element_text(size = 8.5),
        axis.title = element_text(size = 10)) +
  guides(fill = "none",
         size = "none",
         colour = "none") +
  coord_fixed()
print(main.plot)

# put together and save
#pdf(file = "fig-map_of_sampled_points_and_species_range-clownfish.pdf", width = 3.15, height = 8.5)
plot <- cowplot::ggdraw() +
  draw_plot(main.plot) +
  draw_plot(inset.plot, x = 0.495, y = 0.381, width = 0.5, height = 0.5)
#print(plot)
#dev.off()
cowplot::ggsave2(plot = plot, file = "fig-map_of_sampled_points_and_species_range-bombusbifarius.png", 
                 width = 3.15, height = 8.5, units = c("in"), dpi = 600)



