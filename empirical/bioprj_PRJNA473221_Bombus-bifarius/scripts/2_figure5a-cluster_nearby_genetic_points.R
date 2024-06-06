#idea: look at bee genetic points and aggregate into some larger clusters for mapping

#load libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(rnaturalearth)
library(fields)

rm(list=ls())
gc()



#define path to WM_model git repo
headir = "/Users/rachel/WM_model/"

#get a world map
world <- ne_countries(scale = "medium", returnclass = "sf")

#get genetic sample points
genetic <- read.delim(paste0(headir,"empirical/bioprj_PRJNA473221_Bombus-bifarius/PRJNA473221_SraRunTable.txt"), sep  = ",") %>%
  filter(Organism == "Bombus bifarius") %>% 
  separate(lat_lon, into = c("decimalLatitude","dir1","decimalLongitude","dir2"), sep = " ", remove = F) %>% 
  mutate(decimalLongitude = as.numeric(decimalLongitude)*(-1)) %>% 
  mutate(decimalLatitude = as.numeric(decimalLatitude)) %>% 
  dplyr::rename("y"="decimalLatitude", "x"="decimalLongitude")

#add dummy pop names and N samps per pop
genetic <- genetic %>% group_by(x,y) %>% mutate(dummypop = paste0(y,"_",x))
genetic %>% group_by(dummypop) %>% summarise(n=n()) %>% as.data.frame()
genetic <- genetic %>% group_by(dummypop) %>% mutate(n.samps.perpop = n()) %>% as.data.frame()
genetic.locales <- genetic %>% dplyr::select(x,y,dummypop,n.samps.perpop) %>% distinct() %>% 
  arrange(y,x) %>% mutate(clusterid = 1:nrow(.))

#see where we are at viz-wise
ggplot() + 
  #build map
  geom_sf(data = world) +
  #add genetic points
  geom_point(data = genetic.locales, aes(x=x, y=y), shape = 21, color = "black", size = 3, fill = "transparent") +
  scale_y_continuous(limits = c(35, 50), breaks = c(seq(-70,70,2))) +
  scale_x_continuous(limits = c(-127, -115), breaks = c(seq(-170,170,2))) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_bw()

#figure out which clusters are super close to each other and lump for viz purposes
points <- genetic.locales
row.names(points) <- points$clusterid
points <- points %>% dplyr::select(x,y)
gcd <- fields::rdist.earth(points,miles=FALSE) #output in kms
dists <- gcd %>% as.data.frame() %>% mutate(samp = row.names(.)) %>% 
  dplyr::select(samp,everything()) %>% 
  pivot_longer(., 2:ncol(.), names_to = "samp2", values_to = "dist_km") %>% 
  mutate(samp2 = gsub("V","",samp2)) %>% mutate(self = ifelse(samp == samp2, "self","notself")) %>% 
  filter(self == "notself") %>% 
  arrange(dist_km)

#Lump clusters (less than 20 km to each other)
out <- genetic.locales %>% mutate(newclusterid = clusterid)
out$newclusterid[out$clusterid == 14] = 13
out$newclusterid[out$clusterid == 24] = 23
out$newclusterid[out$clusterid == 31] = 30
out$newclusterid[out$clusterid == 7] = 6
out$newclusterid[out$clusterid == 3] = 2
out$newclusterid[out$clusterid == 39] = 37
out$newclusterid[out$clusterid == 16] = 15
out$newclusterid[out$clusterid == 10] = 9
out$newclusterid[out$clusterid == 38] = 37
out$newclusterid[out$clusterid == 5] = 6
out$newclusterid[out$clusterid == 25] = 23
out$newclusterid[out$clusterid == 38] = 37
out$newclusterid[out$clusterid == 32] = 29
out$newclusterid[out$clusterid == 28] = 27
out$newclusterid[out$clusterid == 8] = 6
out$newclusterid[out$clusterid == 4] = 2
out$newclusterid[out$clusterid == 35] = 34
out$newclusterid[out$clusterid == 26] = 27
out$newclusterid[out$clusterid == 18] = 17
out <- out %>% group_by(newclusterid) %>% 
  mutate(x = mean(x), y = mean(y)) %>% 
  mutate(newn.samps.perpop = sum(n.samps.perpop)) %>% 
  as.data.frame()

#visualize how we did with building clusters
ggplot() + 
  #build map
  geom_sf(data = world) +
  #add genetic points
  geom_text(data = out, aes(x=x, y=y, label=newclusterid), hjust=0.5, vjust=0.45, colour = "black", size = 3) +
  scale_y_continuous(limits = c(35, 50), breaks = c(seq(-70,70,2))) +
  scale_x_continuous(limits = c(-127, -115), breaks = c(seq(-170,170,2))) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_bw()

ggplot() + 
  #build map
  geom_sf(data = world) +
  #add genetic points
  geom_point(data = out %>% dplyr::select(x,y,newclusterid,newn.samps.perpop) %>% distinct(), 
             aes(x=x, y=y, size=newn.samps.perpop), 
             shape=21, colour = "black", fill = "transparent") +
  scale_y_continuous(limits = c(35, 50), breaks = c(seq(-70,70,2))) +
  scale_x_continuous(limits = c(-127, -115), breaks = c(seq(-170,170,2))) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme_bw()

#check distances between new clusters
points <- out
points <- points %>% dplyr::select(x,y,newclusterid) %>% distinct()
row.names(points) <- points$newclusterid
dists <- points %>% dplyr::select(x,y)
dists <- fields::rdist.earth(dists, miles=FALSE) %>% as.data.frame() #output in kms
dists <- cbind(points$newclusterid, dists) %>% rename("samp" = 1) %>% 
  dplyr::select(samp, everything())
colnames(dists) <- c("samp", points$newclusterid)
dists <- dists %>% 
  pivot_longer(., 2:ncol(.), names_to = "samp2", values_to = "dist_km") %>% 
  mutate(self = ifelse(samp == samp2, "self","notself")) %>% 
  filter(self == "notself") %>% 
  arrange(dist_km) %>% mutate(samp = as.numeric(samp))
  #good, all further apart than 20km as expected

#save
write.csv(out %>% dplyr::select(x,y,newclusterid,newn.samps.perpop) %>% distinct(), 
          paste0(headir,"empirical/bioprj_PRJNA473221_Bombus-bifarius/figure5a/newgenetic_clusters_for_sample_map.csv"), 
          row.names = FALSE)



