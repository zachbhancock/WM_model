#idea: build pairwise geographic distance matrix between all sequenced individuals from lat/long coordinates to use in WM model

#load libraries
library(dplyr)
library(tidyr)
library(readxl)

rm(list = ls())
gc()


#path to WM_model git repo
headir = "/Users/rachel/WM_model/"


#calculate and save pairwise geographic distances

#get lat/long table
df <- readxl::read_xlsx(paste0(headir,"empirical_data/bioprj_PRJNA473221_Bombus-bifarius/E1590_PRJNA473221_validated_QC.xlsx"),
                        sheet = 2)

#data wrangling for rdist.earth()
geodist <- df %>% dplyr::select(decimalLongitude, decimalLatitude, run_acc_sra) %>%
  dplyr::rename("y"="decimalLongitude", "x"="decimalLatitude") %>% as.data.frame()
rownames(geodist) <- geodist$run_acc_sra
geodist <- geodist %>% dplyr::select(y,x)
geodist <- as.matrix(geodist)

#calc pw great circle distances
gcd <- fields::rdist.earth(geodist, miles=FALSE)

#save
max_and_pw_dists <- list("pw.gcd.genetic" = gcd)

save(max_and_pw_dists, file=paste0(headir,"empirical_data/bioprj_PRJNA473221_Bombus-bifarius/max_and_pw_distsbioprj_PRJNA473221_Bombus-bifarius.Robj"))


