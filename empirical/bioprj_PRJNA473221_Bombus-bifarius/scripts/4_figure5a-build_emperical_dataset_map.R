#idea: make map of species range and sampled genetic points for empirical dataset

#load libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(concaveman)
library(smoothr)
library(legendMap)
#devtools::install_github("3wen/legendMap") # to put N arrow and scale bar on map


rm(list=ls())
gc()

#path to WM_model git repo
headir = "/Users/rachel/WM_model/"

#get a world map
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

#get Great Lakes outlines
lakes <- rnaturalearth::ne_download(scale="large", category = 'physical', type = "lakes", returnclass = "sf")
#subset world lake outlines to just GLs
gls <- c("Lake Superior", "Lake Huron",
         "Lake Michigan", "Lake Erie", "Lake Ontario")
gls <- subset(lakes, name %in% gls)
#how does it look?
ggplot() +
  geom_sf(data = gls)

#get GBIF range points and genetic sample points
genetic <- read.delim(paste0(headir,"empirical/bioprj_PRJNA473221_Bombus-bifarius/PRJNA473221_SraRunTable.txt"), sep  = ",") %>%
  filter(Organism == "Bombus bifarius") %>% 
  separate(lat_lon, into = c("decimalLatitude","dir1","decimalLongitude","dir2"), sep = " ", remove = F) %>% 
  mutate(decimalLongitude = as.numeric(decimalLongitude)*(-1)) %>% 
  mutate(decimalLatitude = as.numeric(decimalLatitude)) %>% 
  dplyr::rename("y"="decimalLatitude", "x"="decimalLongitude")
gbif <- read.csv(paste0(headir,"empirical/bioprj_PRJNA473221_Bombus-bifarius/figure5a/GBIF_Bombus-bifarius_occurrences.csv")) %>% 
  dplyr::rename("y"="decimalLatitude", "x"="decimalLongitude")
genclusts <- read.csv(paste0(headir,"empirical/bioprj_PRJNA473221_Bombus-bifarius/figure5a/genetic_clusters_for_sample_map.csv"))

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

#set bounds for inset map aka zoomed in map of sample locations
minlat=35
maxlat=50
minlon=-125.5
maxlon=-116

#world with box for inset -----
inset.plot <- ggplot() + 
  #build map
  geom_sf(data = world) +
  geom_sf(data = gls, fill = "white") +
  scale_y_continuous(expand = c(0, 0), limits = c(21.5, 77), breaks = c(seq(-70,70,30))) +
  scale_x_continuous(expand = c(0, 0), limits = c(-170, -59), breaks = c(seq(-170,170,30))) +
  #add range polygon
  geom_polygon(data = range, aes(x=x, y=y), fill = "gray40", color = "transparent", alpha = 0.75) +
  #add zoom box
  geom_rect(aes(xmin = minlon, xmax = maxlon, ymin = minlat, ymax = maxlat), colour = "black", fill = "transparent") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        #axis.title = element_blank(),
        #axis.text = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  guides(fill = "none") +
  labs(x = "Longitude",
       y = "Latitude") +
  coord_sf()
print(inset.plot)
ggsave(file = paste0(headir,"empirical/bioprj_PRJNA473221_Bombus-bifarius/figure5a/fig5a-map_of_range-bombusbifarius.pdf"), 
       width = 3, height = 3, units = c("in"), dpi = 600)



#plot zoomed map again with new clusters ----------
main.plot <- ggplot() + 
  #build map
  geom_sf(data = world, fill = "gray40") +
  #add genetic points
  geom_point(data = genclusts, aes(x=x, y=y, size=newn.samps.perpop), shape = 21, color = "black", fill = "goldenrod2") +
  scale_y_continuous(expand = c(0, 0), limits = c(35, 50), breaks = c(seq(-70,70,5))) +
  scale_x_continuous(expand = c(0, 0), limits = c(-125.5, -116), breaks = c(seq(-170,170,5))) +
  scale_size_continuous(breaks= c(5,10,20,40), range = c(1,15)) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 1.5),
        legend.position = c(1.25, 0.11),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, 'mm'),
        legend.title = element_text(size = 10),
        axis.text = element_text(size = 8.5),
        axis.title = element_text(size = 10)) +
  guides(size=guide_legend(title="N. indivs.")) +
  coord_sf() +
  legendMap::scale_bar(lon = -125, lat = 35.25,
                       distance_lon = 150, distance_lat = 20,
                       distance_legend = 40, dist_unit = "km",
                       arrow_length = 100, legend_size = 2.25,
                       arrow_north_size = 4, arrow_distance = 90) 
print(main.plot)
ggsave(file = paste0(headir,"empirical/bioprj_PRJNA473221_Bombus-bifarius/figure5a/fig5a-map_of_sampled_points-bombusbifarius.pdf"), 
        width = 7.25, height = 7.25, units = c("in"), dpi = 600)


#put two plots together in Affinity to make final file: fig5a-map-bombusbifarius.pdf

