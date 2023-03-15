#idea: make map of species range and sampled genetic points

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

rm(list=ls())
gc()

#source world map
load(paste0("empirical_data/countries.outlines.Robj"))

#load (already fetched) bathymetric maps of world (one -180 to 180 longitude and another 0 - 360 longitude)
#bathymap <- getNOAA.bathy(lat1 = -90, lon1 =-180, lat2 = 90, lon2 = 180, resolution = 4, keep = FALSE, antimeridian = FALSE, path = NULL)
#save(bathymap, file = "/Users/rachel/Desktop/DivDiv/mapping_GBIF_and_genetic_sample_points/world_NOAA_bathy_res4.Robj")
#load bathymap of world
load("empirical_data/world_NOAA_bathy_res4.Robj")
bathmap.simple <- bathymap
bathmap.simple[bathmap.simple > 0] <- 5000
bathmap.simple[bathmap.simple <= 0] <- -5000

#get GBIF range points
genetic <- read.delim(paste0("empirical_data/bioprj_PRJNA294760_Amphiprion-bicinctus/lat_long_table-bioprj_PRJNA294760_Amphiprion-bicinctus.txt")) %>% 
  dplyr::rename("y"="lat", "x"="long")
gbif <- read.csv(paste0("empirical_data/bioprj_PRJNA294760_Amphiprion-bicinctus/GBIFLocations_Amphiprion bicinctus.csv")) %>% 
  dplyr::rename("y"="decimalLatitude", "x"="decimalLongitude")

#specific to clownfish
#clean GBIF points to match range map in published paper for dataset
gbif <- gbif %>% filter(y > 10)

#add dummy pop variables and N samps per pop
genetic <- genetic %>% group_by(x,y) %>% mutate(dummypop = paste0(y,"_",x))
genetic %>% distinct(dummypop)
genetic <- genetic %>% group_by(dummypop) %>% mutate(n.samps.perpop = n())


#get zoomed bathymap
maxlat <- max(gbif$y, genetic$y) + 3
minlat <- min(gbif$y, genetic$y) - 3
maxlon <- max(gbif$x, genetic$x) + 3
minlon <- min(gbif$x, genetic$x) - 3
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
               colour="gray70", size=0.2) +
  scale_fill_gradient2(low="white", mid="white", high="gray90", midpoint = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #add zoom box
  geom_rect(aes(xmin = minlon, xmax = maxlon, ymin = minlat, ymax = maxlat), colour = "red", fill = "transparent") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)) +
  guides(fill = "none") +
  coord_fixed()

# range polygon and genetic points - zoomed --------
main.plot <- ggplot() + 
  #build map
  geom_raster(data = bathymapZOOMED.simple, aes(x=x, y=y, fill=z)) +
  geom_contour(data = bathymapZOOMED, aes(x=x, y=y, z=z),
               breaks=0, #contour for land
               colour="gray70", size=0.3) +
  scale_fill_gradient2(low="white", mid="white", high="gray90", midpoint = 0) +
  #add range polygon
  #ggalt::geom_encircle(data = gbif, aes(x=x, y=y), alpha = 0.35, expand=0, color = "blue") +
  #geom_point(data = gbif, aes(x=x, y=y, color = "blue"), alpha = 0.35, size = 1) +
  ggforce::geom_mark_hull(data = gbif, aes(x=x, y=y), concavity = 6, expand=0, radius=0, fill = "blue", colour = "transparent", alpha = 0.25) +
  #ggforce::geom_mark_hull(data = gbif, aes(x=x, y=y), concavity = 6, expand=0, radius=0.001, color = "green") +
  #ggforce::geom_mark_hull(data = gbif, aes(x=x, y=y), color = "green", expand=0.0) +
  #ggforce::geom_mark_hull(data = gbif, aes(x=x, y=y), concavity = 5, expand=0, radius=0.1, color = "blue")
  
  #add genetic points
  geom_point(data = genetic, aes(x=x, y=y, size = n.samps.perpop, color = "darkorange1"), alpha = 1) +
  geom_text(data = genetic, aes(x=x, y=y, label = n.samps.perpop), hjust=0.5, vjust=0.45, colour = "black", size = 3) +
  #manually make a legend for points and lines
  scale_color_identity(guide = "legend",
                       name = "",
                       breaks = c("blue", "darkorange1"),
                       labels = c("Occurrence reports (GBIF)", "Genetic samples")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.88, 0.08),
        legend.background = element_blank()) +
  guides(fill = "none",
         size = "none") +
  coord_fixed()
print(main.plot)
  
# put together and save
pdf(file = "fig-map_of_sampled_points_and_species_range-clownfish.pdf", width = 11, height = 8.5)
plot <- cowplot::ggdraw() +
  draw_plot(main.plot) +
  draw_plot(inset.plot, x = 0.04, y = 0.04, width = 0.25, height = 0.25)
print(plot)
dev.off()





#graveyard ---------------

#GBIF points and genetic points
# points zoomed
main.plot <- ggplot() + 
  #build map
  geom_raster(data = bathymapZOOMED.simple, aes(x=x, y=y, fill=z)) +
  geom_contour(data = bathymapZOOMED, aes(x=x, y=y, z=z),
               breaks=0, #contour for land
               colour="gray70", size=0.3) +
  scale_fill_gradient2(low="white", mid="white", high="gray90", midpoint = 0) +
  #add points
  geom_point(data = gbif, aes(x=x, y=y, color = "blue"), alpha = 0.35, size = 1) +
  geom_point(data = genetic, aes(x=x, y=y, color = "orange"), alpha = 1, size = 1) +
  #manually make a legend for points and lines
  scale_color_identity(guide = "legend",
                       name = "",
                       breaks = c("blue", "orange"),
                       labels = c("Occurrence reports (GBIF)", "Genetic samples")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = c(0.88, 0.08),
        legend.background = element_blank()) +
  guides(fill = "none") +
  coord_fixed()

