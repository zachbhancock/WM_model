#idea: build pairwise geographic distance matrix between all sequenced individuals from lat/long coordinates to use in WM model

#load libraries
library(dplyr)
library(tidyr)
library(fields)

rm(list = ls())
gc()


#path to WM_model git repo
headir = "/Users/rachel/WM_model/"


#calculate and save pairwise geographic distances

#get NCBI metadata table with lat/longs
df <- read.delim(paste0(headir,"empirical/data/bioprj_PRJNA473221_Bombus-bifarius/PRJNA473221_SraRunTable.txt"), sep = ",")

#pull out lat/longs for just B. bifarius
df <- df %>% separate(lat_lon, into = c("decimalLatitude","dir1","decimalLongitude","dir2"), sep = " ", remove = F)
df %>% dplyr::select(lat_lon, decimalLatitude, dir1, decimalLongitude, dir2) %>% group_by(dir1,dir2) %>% summarise(n=n())
df <- df %>% 
  mutate(decimalLongitude = as.numeric(decimalLongitude)*(-1)) %>% 
  mutate(decimalLatitude = as.numeric(decimalLatitude)) %>% 
  rename("run_acc_sra"="Run") %>% 
  filter(Organism == "Bombus bifarius")

#data wrangling for rdist.earth()
geodist <- df %>% dplyr::select(decimalLongitude, decimalLatitude, run_acc_sra) %>%
  dplyr::rename("x"="decimalLongitude", "y"="decimalLatitude") %>% as.data.frame()
rownames(geodist) <- geodist$run_acc_sra
geodist <- geodist %>% dplyr::select(y,x)
geodist <- as.matrix(geodist)

#calc pw great circle distances
gcd <- fields::rdist.earth(geodist, miles=FALSE)

#save
max_and_pw_dists <- list("pw.gcd.genetic" = gcd)

save(max_and_pw_dists, file=paste0(headir,"empirical/data/bioprj_PRJNA473221_Bombus-bifarius/inputs_for_WM/max_and_pw_distsbioprj_PRJNA473221_Bombus-bifarius.Robj"))


