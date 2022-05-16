################################################################
################################################################
#	running inference of wright-malecot model on SLiMulated data
################################################################
################################################################


# load libraries and source relevant functions
library(rstan, warn.conflicts = FALSE, quietly = TRUE)
library(plyr, quietly = TRUE, warn.conflicts = FALSE)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)
library(tidyr, quietly = TRUE, warn.conflicts = FALSE)
#print session info and args for HPCC metadata/log files
print("STARTING TO RUN exe_WM.R SCRIPT")
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("we just opened R and loaded all libraries")
print("session info:")
sessionInfo()
#get and print a vector of all arguments passed from "-export" in submit file
args <- commandArgs(trailingOnly = TRUE)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("we are printing the contents of the vector args to show all variables passed from htcondor to R enviro")
print("arguments passed:")
cat(args, sep = "\n")

#define some variables (pulled in from bash script)
workingdir = args[1]
treefile = args[2] %>% gsub("/", "", .)
print(paste0("workingdir is ", workingdir))
print(paste0("treefile is ", treefile))

#for local testing, example of formatting for treefile variable
#treefile = "/linear_46352997_slimIter_11_sigma_0.2"

#source our functions 
source(paste0(workingdir,"/wm_lib.R"))

stanFile <- "wm_hom_cmpPar_mod_block_scaled.R"
source(paste0(workingdir,"/",stanFile))
ibsMod <- stan_model(model_code=stanBlock)

#print objects loaded in R (so we can see functions and variables were correctly loaded)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("so far we have loaded libraries and printed args")
print("objects currently loaded in R:")
ls()

#print all warnings/errors as they occur
options(warn=1)

#under model
#N <- 30
#coords <- matrix(runif(2*N,0,30),nrow=N,ncol=2)
#k <- 2
#s <- 0.4
#m <- 1e-3
#nbhd <- 15
#inDeme <- 0.5
#nugget <- 0.01
#geoDist <- fields::rdist(coords)
#parIBS <- s + besselK(nu=0,x=m * geoDist) / nbhd
#parIBS[which(geoDist < k)] <- inDeme + s
#diag(parIBS) <- diag(parIBS) + nugget
#plot(geoDist,parIBS)
#db <- list("N" = N,
#           "L" = 5e3,
#           "hom" = parIBS,
#           "k" = k,
#           "geoDist" = geoDist)
#
#ibsMod <- stan_model(model_code=stanBlock)
#fit <- sampling(ibsMod,data=db,iter=1e3,chains=1)


#get list of all square iterations and square sizes to process
list_of_prefixes <- list.files(path = workingdir, pattern = paste0(treefile,"*"), full.names = FALSE) %>% 
  as.data.frame() %>% dplyr::rename("file" = ".") %>% separate(., file, into = c("prefix"), sep = "-", extra = "drop") %>% 
  distinct() %>% filter(stringr::str_detect(prefix, "K") %>% filter(stringr::str_detect(prefix, "sigma"))

#process for each K and sigma combo
for (prefix in list_of_prefixes$prefix){
  
  print(paste0("prefix is: ", prefix)) #prints to .out file
  #also print this to .err so we can tell when errors are thrown if they are
  cat(paste0("prefix is: ", prefix, " now\n"), file = stderr())
  
  # read in slim metadata and make dataBlock for inference
  sampled <- read.table(file=paste0(workingdir, "/", prefix, "-pi_locs.txt"), header=TRUE)
  coords <- sampled[,c("x","y")]
  geoDist <- fields::rdist(coords)
  sampled.coal <- data.matrix(read.table(file=paste0(workingdir, "/", prefix, "-pi.csv"), header=TRUE))
  #for torus
  #geoDist <- dist.torus(coords)
  #geoDist <- as.matrix(geoDist)
  
  
  # translate coalescent times into pairwise pi
  #	note that the coalescent time matrix output by SLiM is 
  #	actually 2*TMRCA, so multiplying by a mutation rate 
  #	should give pairwise pi
  pwp <- sampled.coal
  hom <- 1-pwp
  # jitter it
  #for(i in 1:nrow(pwp)){
  #  for(j in i:nrow(pwp)){
  #    hom[i,j] <- hom[i,j] + rnorm(1,0,1e-7)
  #    hom[j,i] <- hom[i,j]
  #  }
  #}
  
  # add inbreeding
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
        prefix = paste0(prefix, "-est")))
        
    
}

#source our functions 
source(paste0(workingdir,"/wm_lib.R"))

#now run compLik
stanFile <- "wm_hom_cmpPar_cmpLnL_mod_block_scaled.R"

#process each K and sigma combo
for (prefix in list_of_prefixes$prefix){
  
  print(paste0("prefix is: ", prefix)) #prints to .out file
  #also print this to .err so we can tell when errors are thrown if they are
  cat(paste0("prefix is: ", prefix, " now\n"), file = stderr())
  
  # read in slim metadata and make dataBlock for inference
  sampled <- read.table(file=paste0(workingdir, "/", prefix, "-pi_locs.txt"), header=TRUE)
  coords <- sampled[,c("x","y")]
  geoDist <- fields::rdist(coords)
  sampled.coal <- data.matrix(read.table(file=paste0(workingdir, "/", prefix, "-pi.csv"), header=TRUE))
  #for torus
  #geoDist <- dist.torus(coords)
  #geoDist <- as.matrix(geoDist)
  
  
  # translate coalescent times into pairwise pi
  #	note that the coalescent time matrix output by SLiM is 
  #	actually 2*TMRCA, so multiplying by a mutation rate 
  #	should give pairwise pi
  pwp <- sampled.coal
  hom <- 1-pwp
  # jitter it
  #for(i in 1:nrow(pwp)){
  #  for(j in i:nrow(pwp)){
  #    hom[i,j] <- hom[i,j] + rnorm(1,0,1e-7)
  #    hom[j,i] <- hom[i,j]
  #  }
  #}
  
  # add inbreeding
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
  

source(paste0(workingdir,"/",stanFile))
ibsMod <- stan_model(model_code=stanBlock)
    try(runWM(stanMod = ibsMod,
        dataBlock = dataBlock,
        nChains = 4,
        nIter = 4e3,
        prefix = paste0(prefix, "-cmpLnL_est")))
  
}
