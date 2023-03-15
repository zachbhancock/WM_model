################################################################
################################################################
#	running inference of wright-malecot model on SLiMulated data
################################################################
################################################################


#This code is for running the WM model locally. It uses the example data from the repo - just swap out the "sampled" and "sampled.coal" objects with your own matrices!

#source our functions 
source(wm_lib.R)

stanFile <- "/models/wm_hom_cmpPar_mod_block_scaled.R"
source(stanFile))
ibsMod <- stan_model(model_code=stanBlock)

# here, we read in the metadata and convert geographic locations to pairwise geographic distance
sampled <- read.table("example_locs.txt"), header=TRUE)
coords <- sampled[,c("x","y")]
geoDist <- fields::rdist(coords)
sampled.coal <- data.matrix(read.table("example_pwp.csv"), header=TRUE)  
  
#transform pairwise pi into pairwise homozygosity
pwp <- sampled.coal
hom <- 1-pwp
# add inbreeding, this standardizes the diagonal to 1
diag(hom) <- 1
  
# make dataBlock for stan
# N = number of samples
# L = number of loci
# hom = pairwise homozygosity
# k = geographic distance w/in which W-M breaks down 
# (should be ~2*sigma)
# geoDist = pairwise geographic distance
  
dataBlock <- list("N"=nrow(pwp),
                      "L" = 1e4,
                      "hom"=hom,
                      "k" = 0.25,
                      "geoDist"=geoDist)
# run inference
try(runWM(stanMod = ibsMod,
	dataBlock = dataBlock,
	nChains = 4,
	nIter = 4e3,
	prefix = "est_wishart"))