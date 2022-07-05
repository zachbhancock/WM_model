#Load Stan output and make plots
# load libraries and source relevant functions
library(ggplot2, warn.conflicts = FALSE, quietly = TRUE)
library(plyr, warn.conflicts = FALSE, quietly = TRUE)
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
workingdir = args[1]
model_flavor = args[2]
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
list_of_prefixes <- list.files(path = workingdir, pattern = paste0("wmModel_", "*"), full.names = FALSE) %>% 
  as.data.frame() %>% dplyr::rename("file" = ".") %>% separate(., file, into = c("prefix"), sep = "-", extra = "drop") %>% 
  distinct() %>% filter(stringr::str_detect(prefix, "K"))

#process each square iteration / size combo
for (prefix in list_of_prefixes$prefix) {
  
  print(paste0("prefix is: ", prefix))  
  
  #get label info (added by Rachel / modified by Zach)
  labels <- prefix %>% as.data.frame() %>% dplyr::rename("prefix"=".") %>% 
    separate(., prefix, into = c("slurm_job_id","slimIter","sigma","K"), sep = "_") %>% 
    dplyr::select(slurm_job_id, slimIter, sigma, K)
  
  #read dist files
  gen.dist <- data.matrix(read.table(file=paste0(prefix, "-pi.csv"), header=TRUE))
  gen.dist <- as.matrix(gen.dist)
  gen.dist <- gen.dist[lower.tri(gen.dist)]
  gen.dist <- as.data.frame(gen.dist)
  geo.dist <- read.table(file=paste0(prefix, "-pi_locs.txt"), header = TRUE)
  coords.dist <- geo.dist[,c("x","y")]
  geo.dist <- fields::rdist(coords.dist)
  geo.dist <- as.matrix(geo.dist)
  geo.dist <- geo.dist[lower.tri(geo.dist)]
  geo.dist <- as.data.frame(geo.dist)
  dist.all <- cbind(gen.dist, geo.dist)
  names(dist.all)[1] <- "gen.dist"
  names(dist.all[2] <- "geo.dist"
  
  load(paste(prefix, "-est_"[wishart, cmpLnL]"_out.Robj", sep=""), verbose=TRUE)
  col_pi <- rstan::extract(out$fit, "s", inc_warmup=TRUE, permute=FALSE)
  nbhd_pi <- rstan::extract(out$fit, "nbhd", inc_warmup=TRUE, permute=FALSE)
  inDeme_pi <- rstan::extract(out$fit, "inDeme", inc_warmup=TRUE, permute=FALSE)
  m_pi <- rstan::extract(out$fit, "m", inc_warmup=TRUE, permute=FALSE)
  post_pi <- rstan::get_logposterior(out$fit,inc_warmup=FALSE)
  post_pi <- data.frame(matrix(unlist(post_pi), nrow=length(post_pi), byrow=FALSE))
  post_pi <- cbind(stack(post_pi[1:4]))
  names(post_pi)[1] <- "posterior"
  col_pi <- plyr::adply(col_pi, c(1,2,3))
  wmModel_pi <- cbind(col_pi, post_pi, nbhd_pi, inDeme_pi, m_pi) %>% 
    dplyr::mutate(col_pi = 1 - V1, iterations = as.numeric(iterations)) %>% 
    dplyr::select(-chains)
  wmModel_pi$theo_nbhd <- 4*pi*wmModel_pi$K*(wmModel_pi*sigma)^2
  wmModel_pi$delta_nbhd <- theo.nbhd - wmModel_pi$nbhd_pi
  wmModel_pi$model <- "Wischart"
  
  #load(paste(prefix, "-est_cmpLnL_out.Robj", sep=""), verbose=TRUE)
  #col_cmpLnL <- rstan::extract(out$fit, "s", inc_warmup=TRUE, permute=FALSE)
  #nbhd_cmpLnL <- rstan::extract(out$fit, "nbhd", inc_warmup=TRUE, permute=FALSE)
  #inDeme_cmpLnL <- rstan::extract(out$fit, "inDeme", inc_warmup=TRUE, permute=FALSE)
  #m_cmpLnL <- rstan::extract(out$fit, "m", inc_warmup=TRUE, permute=FALSE)
  #post_cmpLnL <- rstan::get_logposterior(out$fit,inc_warmup=FALSE)
  #post_cmpLnL <- data.frame(matrix(unlist(post_cmpLnL), nrow=length(post_cmpLnL), byrow=FALSE))
  #post_cmpLnL <- cbind(stack(post_cmpLnL[1:4]))
  #names(post_LnL)[1] <- "posterior"
  #col_cmpLnL <- plyr::adply(col_cmpLnL, c(1,2,3))
  #wmModel_cmpLnL <- cbind(col_cmpLnL, post_cmpLnL, nbhd_cmpLnL, inDeme_cmpLnL, m_cmpLnL) %>% 
  #  dplyr::mutate(col_cmpLnL = 1 - V1, iterations = as.numeric(iterations)) %>% 
  #  dplyr::select(-chains)
  #wmModel_cmpLnL$theo_nbhd <- 4*pi*wmModel_cmpLnL$K*(wmModel_cmpLnL*sigma)^2
  #wmModel_cmpLnL$delta_nbhd <- theo.nbhd - wmModel_pi$nbhd_pi
  #wmModel_cmpLnL$model <- "cmpLnL"
  
  wmModel_all <- cbind(wmModel_pi, wmModel_cmpLnL)
  
}

write.table(wmModel_pi, file=paste(treefile, "-wmModel_est.txt", sep=""))

pdf(file=paste(prefix, "-nhbd_K.pdf", sep=""))
est.plot <- ggplot(wmModel_all, aes(x=K, y=nbhd_pi, fill=sigma, shape=model)) + geom_point(size=2)
print(est.plot)
dev.off()

pdf(file=paste(prefix, "-nhbd_delta_K.pdf", sep=""))
delta.plot <- ggplot(wmModel_all, aes(x=K, y=delta_nbhd, fill=sigma, shape=model)) + geom_point(size=2)
print(delta.plot)
dev.off()

pdf(file=paste(prefix, "-col_K.pdf", sep=""))
col_K.plot <- ggplot(wmModel_all, aes(x=K, y=col_pi, fill=sigma, shape=model)) + geom_point(size=2)
print(col_K.plot)
dev.off()

pdf(file=paste(prefix, "-inDeme_K.pdf", sep=""))
col_K.plot <- ggplot(wmModel_all, aes(x=K, y=inDeme_pi, fill=sigma, shape=model)) + geom_point(size=2)
print(col_K.plot)
dev.off()
