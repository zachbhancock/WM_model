#idea:	running inference of wright-malecot model on empirical data

# load libraries ---------------
library(rstan)
library(dplyr)
library(tidyr)

rm(list=ls())
gc()



#define some variables for empirical dataset --------------
functiondir = "/Users/rachel/WM_model/" #path to where functions live

popgenfile = "/Users/rachel/WM_model/empirical/bioprj_PRJNA473221_Bombus-bifarius/inputs_for_WM/popgenstats.0.5.bioprj_PRJNA473221_Bombus-bifarius_stacks_littlem_3_bigm_3_n3_nequalM_popgenstats.Robj" #pwp matrix in this file
sampkeyfile = "/Users/rachel/WM_model/empirical/bioprj_PRJNA473221_Bombus-bifarius/inputs_for_WM/samplenamekey.txt" #sample name matching btwn bioinf output and NCBI
nlocifile = "/Users/rachel/WM_model/empirical/bioprj_PRJNA473221_Bombus-bifarius/inputs_for_WM/bpstats.0.5.bioprj_PRJNA473221_Bombus-bifarius_stacks_littlem_3_bigm_3_n3_nequalM_BPstats.Robj" #Nloci in this file
geogdistfile = "/Users/rachel/WM_model/empirical/bioprj_PRJNA473221_Bombus-bifarius/inputs_for_WM/max_and_pw_distsbioprj_PRJNA473221_Bombus-bifarius.Robj"

outdir = "/Users/rachel/WM_model/empirical/bioprj_PRJNA473221_Bombus-bifarius/outputs_from_WM" #outdir

#create outdir, source functions, load models
dir.create(outdir)

source(paste0(functiondir,"/wm_lib.R"))
source(paste0(functiondir,"/models/wm_hom_cmpPar_mod_block_scaled.R"))

ibsMod <- stan_model(model_code=stanBlock)



#get model inputs (pwp, geog dist, Nloci) --------------
sampkey <- read.delim(sampkeyfile)
  
#get number of loci (L)
load(nlocifile)
Nloci = BPstats$nLoci

#get pw geographic distance matrix
load(geogdistfile)
geoDist <- max_and_pw_dists$pw.gcd.genetic

#rename geoDist from SRR IDs to sampleX IDs
geoDistnames <- colnames(geoDist) %>% as.data.frame() %>% 
  dplyr::rename("run_acc_sra" = ".") %>% 
  dplyr::mutate(order = 1:dplyr::n()) %>% 
  base::merge(., sampkey %>% dplyr::select(run_acc_sra,sampid_assigned_for_bioinf), 
              by = "run_acc_sra", all.x = T) %>% 
  arrange(order)
geoDistnames <- geoDistnames$sampid_assigned_for_bioinf
colnames(geoDist) <- geoDistnames
rownames(geoDist) <- geoDistnames

#get pw pi aka genetic matrix
load(popgenfile)
pwp <- popgenstats$pwp

#get list of just samples in both geog and pwp matrices
sampstouse <- merge(rownames(geoDist) %>% as.data.frame() %>% mutate(x = "geo") %>% dplyr::rename("samp" = "."),
                    rownames(pwp) %>% as.data.frame() %>% mutate(y = "geno") %>% dplyr::rename("samp" = "."),
                    by = "samp", all.x = T, all.y = T) %>% 
  filter(is.na(x)==F & is.na(y)==F) %>% 
  arrange(samp)
sampstouse <- as.character(sampstouse$samp)
#keep just samps in both geog and pwp matrices
geoDist <- geoDist[rownames(geoDist) %in% sampstouse, colnames(geoDist) %in% sampstouse]
pwp <- pwp[rownames(pwp) %in% sampstouse, colnames(pwp) %in% sampstouse]
#and order pwp and geogDist the same/to match
geoDist.ordered <- geoDist[sampstouse, sampstouse]
geoDist <- geoDist.ordered
pwp.ordered <- pwp[sampstouse, sampstouse]
pwp <- pwp.ordered

#check ordered same now
if (identical(row.names(pwp),row.names(geoDist)) == TRUE) {
    print("names and order of rownames for pwp and geoDist match") 
  } else {
    print("ERROR ! - names and/or order of rownames for pwp and geoDist do not match!")
  }
  
#convert pwp to hom
hom <- 1-pwp
diag(hom) <- 1 # add inbreeding
flag <- 0
scl_min <- min(hom)
scl_max <- max(hom-scl_min)
tmphom <- (hom-scl_min)/scl_max

#check if hom is positive definite
posdefcheck <- any(eigen(tmphom)$values < 0)
while(posdefcheck & flag < 100){
    diag(hom) <- diag(hom) + 1e-3
    tmphom <- (hom-scl_min)/scl_max
    posdefcheck <- any(eigen(tmphom)$values < 0)
    flag <- flag + 1
    }
if(posdefcheck){
    cat("ERROR: I couldn't get the dataset to be positive definite\n")
  } else {
    print("dataset is positive definite")
  }
  


#run model -----------------

# make dataBlock for stan
# N = number of samples
# L = number of loci
# hom = pairwise homozygosity
# k = geographic distance w/in which W-M breaks down, (should be ~2*sigma)
# geoDist = pairwise geographic distance
dataBlock <- list("N"=nrow(pwp),
                  "L" = Nloci,
                  "hom"=hom,
                  "k" = 1,
                  "geoDist"=geoDist)
  
# run inference
runWM(stanMod = ibsMod,
      dataBlock = dataBlock,
      nChains = 5,
      nIter = 5e3,
      prefix = paste0("WMfit-wishart-bioprj_PRJNA473221_Bombus-bifarius_stacks_littlem_3_bigm_3_n3_nequalM"),
      MLjumpstart = FALSE
      )




