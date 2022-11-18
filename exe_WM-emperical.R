################################################################
################################################################
#	running inference of wright-malecot model on empirical data
################################################################
################################################################


# load libraries and source relevant functions
suppressMessages(library(rstan, warn.conflicts = FALSE, quietly = TRUE))
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
run_name = args[1] #dataset name
indir = args[3] #indir
outdir = args[4] #outdir
minPropIndivsScoredin = as.numeric(args[5]) #percent of indivs that locus must be scored in to save
model_flavor = args[6] #flavor of WM model to run

#for local testing
# indir="../scripts/"
# outdir="../troubleshooting"
# run_name="bioprj_PRJNA294760_Amphiprion-bicinctus"
# minPropIndivsScoredin = 0.5

#source our functions/load models
source(paste0(indir,"/wm_lib.R"))

if (model_flavor == "wishart") { stanFile <- "wm_hom_cmpPar_mod_block_scaled.R" }
if (model_flavor == "cmplnl") { stanFile <- "wm_hom_cmpPar_cmpLnL_mod_block.R" }
print(paste0("using stanFile: ", stanFile))
source(paste0(indir,"/",stanFile))
ibsMod <- stan_model(model_code=stanBlock)

#print objects loaded in R (so we can see functions and variables were correctly loaded)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("so far we have loaded libraries and printed args")
print("objects currently loaded in R:")
ls()

#print all warnings/errors as they occur
options(warn=1)


#get pwp and distance matrices --------------

popgenfiles_list <- list.files(path = paste0(indir), pattern = paste0("popgenstats.",minPropIndivsScoredin), full.names = TRUE)
print("popgenfiles_list contains:")
print(popgenfiles_list)

for (loop.iter in 1:length(popgenfiles_list)) {
  
  #get bookkeeping variables
  popgenfile = popgenfiles_list[loop.iter]   
  stacksparams = popgenfile %>% strsplit(., split = "/") %>% as.data.frame() %>% .[nrow(.),] %>% 
    gsub("_popgenstats|Robj","",.) %>% gsub(run_name,"",.) %>% gsub(minPropIndivsScoredin,"",.) %>% 
    gsub("\\.","",.) %>% gsub("popgenstats_stacks_","",.)
  print(paste0("working on run_name ", run_name, " and Stacks params ", stacksparams))
  #also print this to .err so we can tell when errors are thrown if they are
  cat(paste0("working on run_name ", run_name, " and Stacks params ", stacksparams, " now\n"), file = stderr())
  
  #get sample names
  sampkey <- read.delim(paste0(indir,"/samplenamekey.txt"))
  
  #get number of polymorphic loci
  load(paste0(indir, "/bpstats.", minPropIndivsScoredin, ".", run_name, "_stacks_", stacksparams, "_BPstats.Robj"))
  Npolyloci = BPstats$nLoci
  cat("here1\n", file = stderr())

  #get pw geographic distance matrix
  load(paste0(indir, "/max_and_pw_dists", run_name, ".Robj", sep = ""))
  #geoDist <- max_and_pw_dists$pw.gcd.genetic
  geoDist <- max_and_pw_dists$pw.seadist.genetic
  if ( is.null(geoDist)==T ) {
    geoDist <- max_and_pw_dists$pw.gcd.genetic
    print("using great circle distance")
  } else {
    print("using sea distance")
  }
  # rename geoDist from SRR IDs to sampleX IDs
  geoDistnames <- colnames(geoDist) %>% as.data.frame() %>% dplyr::rename("run_acc_sra" = ".") %>% 
    dplyr::mutate(order = 1:dplyr::n()) %>% 
    base::merge(., sampkey %>% dplyr::select(run_acc_sra,sampid_assigned_for_bioinf), by = "run_acc_sra", all.x = T) %>% 
    arrange(order)
  geoDistnames <- geoDistnames$sampid_assigned_for_bioinf
  colnames(geoDist) <- geoDistnames
  rownames(geoDist) <- geoDistnames
  cat("here2\n", file = stderr()) 

  #get pw pi aka genetic matrix
  load(popgenfile)
  #pwp <- popgenstats$pwp$pwp
  pwp <- popgenstats$pwp
  cat("here3\n", file = stderr())

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
  cat("here4\n", file = stderr())
  #check ordered same now
  if (identical(row.names(pwp),row.names(geoDist)) == TRUE) {
    print("names and order of rownames for pwp and geoDist match")
  } else {
    print("ERROR ! - names and/or order of rownames for pwp and geoDist do not match!")
  }
  cat("here5\n", file = stderr())
  
  #convert pwp to hom
  hom <- 1-pwp
  cat("here6\n", file = stderr())
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
  cat("here7\n", file = stderr())
  #se <- popgenstats$pwp$se
  
  #run model
  if (model_flavor == "wishart") {
    
    # make dataBlock for stan
    # N = number of samples
    # L = number of loci
    # hom = pairwise homozygosity
    # k = geographic distance w/in which W-M breaks down, (should be ~2*sigma)
    # geoDist = pairwise geographic distance
    dataBlock <- list("N"=nrow(pwp),
                      "L" = Npolyloci,
                      "hom"=hom,
                      "k" = 1,
                      "geoDist"=geoDist)
    # run inference
    try(runWM(stanMod = ibsMod,
              dataBlock = dataBlock,
              nChains = 5,
              nIter = 5e3,
              prefix = paste0("WMfit",model_flavor,"-",run_name,"_stacks_",stacksparams),
              MLjumpstart = FALSE))
    
  }
  
  if (model_flavor == "cmplnl") {
    
    ut = upper.tri(hom, diag=FALSE)
    
    dataBlock <- list("lut" = length(hom[ut]),
                      "hom" = hom[ut],
                      "k" = 0.25,
                      "geoDist" = geoDist[ut],
                      "se" = se[ut])
    # run inference
    try(runWM_cmpLnl(stanMod = ibsMod,
                     dataBlock = dataBlock,
                     nChains = 3,
                     nIter = 4e3,
                     prefix = paste0("WMfit",model_flavor,"-",run_name,"_stacks_",stacksparams)))
    
  }

  #clean up at end of loop iteration in case something goes awry
  rm(popgenfile,stacksparams,sampkey,Npolyloci,
     geoDist,geoDistnames,geoDist.ordered,
     pwp.ordered,hom,pwp,sampstouse)

}


