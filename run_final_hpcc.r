#Load Stan output and make plots
# load libraries and source relevant functions
library(ggplot2, warn.conflicts = FALSE, quietly = TRUE)
library(plyr, warn.conflicts = FALSE, quietly = TRUE)
library(tidyr, warn.conflicts = FALSE, quietly = TRUE)
library(dplyr, warn.conflicts = FALSE, quietly = TRUE)

#print session info and args for HPCC metadata/log files
print("STARTING TO RUN collectpi_plots.R SCRIPT")
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
workingdir = compiled_output
model_flavor = Wishart
#treefile = args[2] %>% gsub("/", "", .)
#note to Zach for testing, run: treefile="linear_42014183_slimIter_0_sigma_0.5"
#working dir is just the directory where "all our stuff" for one job (aka one value of sigma and one main slim iteration number) lives

#print objects loaded in R (so we can see functions and variables were correctly loaded)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("so far we have loaded libraries and printed args")
print("objects currently loaded in R:")
ls()

#get list of all K / sigma combos to process
list_of_prefixes <- list.files(path = workingdir, pattern = paste0("*_out.Robj"), full.names = FALSE) %>% 
  as.data.frame() %>% dplyr::rename("Robjfile" = ".") %>% 
  tidyr::separate(., Robjfile, into = c("prefix","temp"), sep = "-", extra = "drop", remove = FALSE) %>% 
  distinct() %>% 
  tidyr::separate(., temp, into = c("garbage","model_flavor"), sep = "_", extra = "drop") %>% dplyr::select(-garbage) %>%  distinct()

#process each square iteration / size combo
wmModel_out <- data.frame("slurm_job_id"=NA, "slimIter"=NA, "sigma"=NA, "K"=NA,
                          "model_flavor"=NA, "chain"=NA, "iteration"=NA,
                          "col_pi"=NA, "nbhd"=NA, "inDeme"=NA, "m"=NA, "posterior"=NA,
                          "theo_nbhd"=NA, "delta_nbhd"=NA)

for (loop.iter in 1:length(list_of_prefixes$Robjfile)) {
  
  #define some variables
  Robjfile = list_of_prefixes[loop.iter, ]$Robjfile
  prefix = list_of_prefixes[loop.iter, ]$prefix
  model_flavor = list_of_prefixes[loop.iter, ]$model_flavor
  
  print(paste0("starting file: ", Robjfile))
  
  labels <- prefix %>% as.data.frame() %>% dplyr::rename("prefix"=".") %>%
    tidyr::separate(., prefix, into = c("garbage","slurm_job_id","garbage2","slimIter","garbage3","sigma","garbage4","K"), sep = "_") %>%
    mutate(model_flavor = model_flavor) %>% 
    mutate(K = as.numeric(K), sigma = as.numeric(sigma)) %>% 
    dplyr::select(slurm_job_id, slimIter, sigma, K, model_flavor)

  #read dist files
  K.density <- read.table(file=paste0(prefix, "_density"), header = FALSE)
  K.density <- K.density$V1
  MK_eff_sigma <- read.table(file=paste0(prefix, "-eff_disp.txt"), header = FALSE)
  MK_eff_sigma <- MK_eff_sigma$V1
  eff_density <- read.table(file=paste0(prefix, "-eff_dens.txt"), header=FALSE)
  eff_density <- eff_density$V1
  gen.dist <- data.matrix(read.table(file=paste0(prefix, "-pi.csv"), header=TRUE))
  gen.dist <- as.matrix(gen.dist)
  gen.dist <- gen.dist[lower.tri(gen.dist)]
  gen.dist <- as.data.frame(gen.dist)
  fst <- data.matrix(read.table(file=paste0(prefix, "-Fst.csv"), header=TRUE))
  fst <- as.matrix(fst)
  fst <- fst[lower.tri(fst)]
  fst <- as.data.frame(fst)
  geo.dist <- read.table(file=paste0(prefix, "-pi_locs.txt"), header = TRUE)
  coords.dist <- geo.dist[,c("x","y")]
  geo.dist <- fields::rdist(coords.dist)
  geo.dist <- as.matrix(geo.dist)
  geo.dist <- geo.dist[lower.tri(geo.dist)]
  geo.dist <- as.data.frame(geo.dist)
  dist.all <- cbind(gen.dist, geo.dist)
  dist.fst <- cbind(fst, geo.dist)
  names(dist.all)[1] <- "gen.dist"
  names(dist.all)[2] <- "geo.dist"
  names(dist.fst)[1] <- "fst"
  names(dist.fst)[2] <- "geo.dist"
  
  #Rousset method for estimating Nb
  dist.fst$y <- dist.fst$fst / (1 - dist.fst$fst)
  reg <- lm(dist.fst$y ~ dist.fst$geo.dist)
  beta <- summary(reg)$coefficients["x", 1]
  Rousset_Nb <- 1 / beta
  
  #read WM_model output
  load(paste0(workingdir, "/", Robjfile), verbose=TRUE)
  col_pi <- rstan::extract(out$fit, "s", inc_warmup=TRUE, permute=FALSE)
  nbhd_pi <- rstan::extract(out$fit, "nbhd", inc_warmup=TRUE, permute=FALSE)
  inDeme_pi <- rstan::extract(out$fit, "inDeme", inc_warmup=TRUE, permute=FALSE)
  m_pi <- rstan::extract(out$fit, "m", inc_warmup=TRUE, permute=FALSE)
  post_pi <- rstan::get_logposterior(out$fit, inc_warmup=FALSE)
  #make long
  col_pi <- plyr::adply(col_pi, c(1,2,3)) %>% 
    dplyr::rename("iteration"="iterations", "chain"="chains") %>% 
    mutate(chain = gsub("chain:","",chain)) %>%
    mutate(col_pi = 1 - V1) %>% dplyr::select(-parameters,-V1)
  nbhd_pi <- plyr::adply(nbhd_pi, c(1,2,3)) %>% 
    dplyr::rename("iteration"="iterations", "chain"="chains", "nbhd"="V1") %>% 
    mutate(chain = gsub("chain:","",chain)) %>% dplyr::select(-parameters)
  inDeme_pi <- plyr::adply(inDeme_pi, c(1,2,3)) %>% 
    dplyr::rename("iteration"="iterations", "chain"="chains", "inDeme"="V1") %>% 
    mutate(chain = gsub("chain:","",chain)) %>% dplyr::select(-parameters)
  m_pi <- plyr::adply(m_pi, c(1,2,3)) %>% 
    dplyr::rename("iteration"="iterations", "chain"="chains",  "m"="V1") %>% 
    mutate(chain = gsub("chain:","",chain)) %>% dplyr::select(-parameters)
  post_pi <- data.frame(matrix(unlist(post_pi), nrow=length(post_pi), byrow=FALSE)) %>% t(.) %>% as.data.frame() %>% 
    mutate(iteration = gsub("X","",row.names(.))) %>%
    pivot_longer(., names_to = "chain", values_to = "posterior", cols = 1:4) %>% 
    mutate(chain = gsub("V","",chain))
  labels_pi <- do.call("rbind", replicate(nrow(col_pi), labels, simplify = FALSE))
  #bind
  wmModel_pi <- cbind(labels_pi, col_pi)
  wmModel_pi <- wmModel_pi %>% 
    merge(., nbhd_pi, by = c("iteration","chain")) %>% 
    merge(., inDeme_pi, by = c("iteration","chain")) %>% 
    merge(., m_pi, by = c("iteration","chain")) %>% 
    merge(., post_pi, by = c("iteration","chain"))
  wmModel_pi <- wmModel_pi %>% 
    mutate(theo_nbhd = 4*pi*(K.density)*(sigma^2)) %>% 
    mutate(delta_nbhd = theo_nbhd - nbhd)
  
  wmModel_pi <- wmModel_pi %>% dplyr::select(slurm_job_id, slimIter, sigma, K,
                                             model_flavor, chain, iteration,
                                             col_pi, nbhd, inDeme, m, posterior,
                                             theo_nbhd, delta_nbhd)
  
  wmModel_out <- rbind(wmModel_out, wmModel_pi)
  wmModel_out$Rousset_Nb <- Rousset_Nb
  wmModel_out$MK_eff_nbhd <- 4*pi*(K.density)*MK_eff_sigma
  wmModel_out$theo_eff_nbhd <- 4*pi*(K.density)*(sigma*sqrt(3/2))^2
  wmModel_out$MK_eff_dens_nbhd <- 4*pi*(eff_density)*MK_eff_sigma
  wmModel_out$theo_eff_dens_nbhd <- 4*pi*(eff_density)*(sigma*sqrt(3/2))^2
  
  
}

wmModel_pi <- wmModel_out %>% filter(is.na(slurm_job_id)==FALSE)

#save text output
write.table(wmModel_pi, file=paste0(prefix, "-est.txt"))

#end



