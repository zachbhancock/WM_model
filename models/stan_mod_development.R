library(rstan)

source("models/wm_hom_cmpPar_mod_block_scaled.R")

mod <- stan_model(model_code=stanBlock)

source('~/Dropbox/spatial_models/WM_model/wm_lib.R')
coords <- data.matrix(read.table("~/Dropbox/bedassle-paper/sims/slimulations/ibd/ibd1_coords.txt"))
pwp <- data.matrix(read.table("~/Dropbox/bedassle-paper/sims/slimulations/ibd/ibd1_pwp.txt"))
geoDist <- fields::rdist(coords)
hom <- 1-pwp
diag(hom) <- 1
db <- list("N" = nrow(coords),
		   "L" = 1e4,
		   "hom" = hom,
		   "k" = 1,
		   "geoDist" = geoDist)


#ml2init(db=db,mod=mod,nRuns=1e1)

fit <- runWM(stanMod=mod,dataBlock=db,nChains=2,nIter=1e3,prefix="test")

fit2 <- runWM(stanMod=mod,dataBlock=db,nChains=2,nIter=1e3,prefix="testMLjumpstart",MLjumpstart=TRUE,nMLruns=20)