#idea: make map of species range and sampled genetic points - SPECIFIC FOR CLOWNFISH

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
#devtools::install_github("3wen/legendMap") # to put N arrow and scale bar on map

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

#build polygon for range that will cover area we want
range <- data.frame(y = c(29.966984540136707,29.59510655999322,15.594107911593486,18.805084570680975,14.990062473658979,13.680852758692668,11.579769771626555,10.025826209455921,11.799981379838396,20.471441358647393,17.20610054229364,12.53016144228723,11.954411104179366,13.963804784857466), 
                    x = c(31.8084972701229,34.98629239993202,52.23208084432265,36.70839937740826,39.75161228440387,41.9488778666583,42.52016688186535,44.695459909731525,51.300645259909096,46.45593469754394,48.097139470917966,55.25110899588169,55.04069812750039,53.567822048831395)
                    ) %>% 
  mutate(sampid = 1:nrow(.))


#get zoomed bathymap
maxlat <- max(gbif$y, genetic$y) + 3
minlat <- min(gbif$y, genetic$y) - 7
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
               colour="gray70", size=0.15) +
  scale_fill_gradient2(low="white", mid="white", high="gray90", midpoint = 0) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
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
#print(inset.plot)



# range polygon and genetic points - zoomed --------
#google maps ocean color: #89b4f9
main.plot <- ggplot() + 
  #add range polygon
  #ggforce::geom_mark_hull(data = gbif, aes(x=x, y=y), concavity = 6, expand=0, radius=0, fill = "blue", colour = "transparent", alpha = 0.25) +
  #geom_point(data = range, aes(x=x, y=y), alpha = 1, colour = "magenta") +
  #geom_text(data = range, aes(x=x, y=y, label = sampid), hjust=0.5, vjust=0.45, colour = "black", size = 3) +
  ggforce::geom_mark_hull(data = range, aes(x=x, y=y), concavity = 0, expand=0, radius=0, fill = "#e6b8f5", colour = "transparent", alpha = 0.75) +

  #build map
  geom_raster(data = bathymapZOOMED.simple, aes(x=x, y=y, fill=z)) +
  geom_contour(data = bathymapZOOMED, aes(x=x, y=y, z=z),
               breaks=0, #contour for land
               colour="gray70", size=0.3) +
  scale_fill_gradient2(low="transparent", mid="transparent", high="gray90", midpoint = 0) +

  #add genetic points
  geom_point(data = genetic, aes(x=x, y=y, size = n.samps.perpop), alpha = 1, color = "darkorange1") +
  geom_text(data = genetic, aes(x=x, y=y, label = n.samps.perpop), hjust=0.5, vjust=0.45, colour = "black", size = 3) +
  
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
  labs(x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  ggsn::north(x.min = minlat, x.max = maxlat, y.min = minlon, y.max = maxlon, symbol = 3) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", size = 1.5),
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
cowplot::ggsave2(plot = plot, file = "fig-map_of_sampled_points_and_species_range-clownfish.png", 
                 width = 3.15, height = 8.5, units = c("in"), dpi = 600)



