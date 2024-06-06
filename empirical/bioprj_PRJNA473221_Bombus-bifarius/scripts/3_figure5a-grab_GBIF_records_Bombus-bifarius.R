#idea: pull all occurrences records for species X from GBIF

#from http://ropensci.org/tutorials/rgbif_tutorial.html
#if on linux and this fails because rgeos fails, see https://stackoverflow.com/questions/38924767/error-installing-r-package-for-linux
#also, make sure ggplot2 is installed, because gbifmap depends on it

#load libraries
library(rgbif)
library(ggplot2)
library(dplyr)
library(rnaturalearthdata)
library(rnaturalearth)


rm(list = ls())
gc()



# define some variables ---------------
#path to WM_model git repo
headir = "/Users/rachel/WM_model/"

#replace XXXXX with GBIF account info 
user<-'XXXXXXX' 
password<- 'XXXXXXX'
email<-'XXXXXXX'

#get taxon key for species of interest aka little brown bats
taxonKey <- rgbif::name_backbone(name="Bombus bifarius")$speciesKey



#build search query for GBIF ---------------
searchRequest <- rgbif::occ_download(type='and',pred('taxonKey',taxonKey),
                                     pred("hasGeospatialIssue", FALSE),
                                     pred("hasCoordinate", TRUE),
                                     pred_not(pred("basisOfRecord", "FOSSIL_SPECIMEN")),
                                     pred_or(
                                       pred_not(pred("establishmentMeans","MANAGED")),
                                       pred_not(pred_notnull("establishmentMeans"))
                                     ),
                                     pred_or(
                                       pred("occurrenceStatus","PRESENT"),
                                       pred_not(pred_notnull("occurrenceStatus"))
                                     ),
                                     format='SIMPLE_CSV',
                                     user=user,pwd=password,email=email)

#run to find status and download data
metaStatus = rgbif::occ_download_meta(searchRequest)
print(metaStatus)
print('waiting')
rgbif::occ_download_wait(searchRequest)
print('done waiting')
whatGot = rgbif::occ_download_get(searchRequest,overwrite=TRUE)
datGBIF = rgbif::occ_download_import(x=whatGot)
datGBIF <- datGBIF %>% mutate(gbifID = as.character(gbifID))



# do some data clean-up --------------
#view on map quick
world <- ne_countries(scale = "medium", returnclass = "sf")

datGBIF %>% 
  ggplot() + 
  geom_sf(data = world) +
  geom_point(aes(x=decimalLongitude, y=decimalLatitude)) +
  scale_y_continuous(limits = c(-79, 79), breaks = c(seq(-79,79,10))) +
  scale_x_continuous(limits = c(-179, 0))

#remove outliers
keep <- datGBIF %>% filter(decimalLongitude > -150) %>% filter(decimalLongitude < -102.5) %>% 
  filter(!gbifID %in% c("42586830", "42585999", "42585857", 
                        "42585846", "42585835", "42585831", 
                        "42585827", "42585823","658149114", 
                        "657672348", "42587995", "42585999", 
                        "1039184721"))

keep %>% 
  ggplot() + 
  geom_sf(data = world, fill = "transparent") +
  geom_point(aes(x=decimalLongitude, y=decimalLatitude)) +
  scale_y_continuous(limits = c(21, 70), breaks = c(seq(-70,70,10))) +
  scale_x_continuous(limits = c(-170, -80), breaks = c(seq(-170,170,10))) +
  geom_vline(xintercept = -105, colour = "red") +
  geom_hline(yintercept = 48, colour = "red")

#save
write.csv(keep, paste0(headir,"empirical/bioprj_PRJNA473221_Bombus-bifarius/figure5a/GBIF_Bombus-bifarius_occurrences.csv"), row.names = FALSE)


