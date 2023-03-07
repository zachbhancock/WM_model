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
#treefile = args[2] %>% gsub("/", "", .)
print(paste0("workingdir is ", workingdir))
#print(paste0("treefile is ", treefile))

#for local testing, example of formatting for treefile variable
#treefile = "/linear_46352997_slimIter_11_sigma_0.2"

#source our functions 

#print objects loaded in R (so we can see functions and variables were correctly loaded)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("so far we have loaded libraries and printed args")
print("objects currently loaded in R:")
ls()

#print all warnings/errors as they occur
options(warn=1)

#get list of all square iterations and square sizes to process
list_of_prefixes <- list.files(path = workingdir, pattern = paste0("wmModel_", "*"), full.names = FALSE) %>% 
  as.data.frame() %>% dplyr::rename("file" = ".") %>% separate(., file, into = c("prefix"), sep = "-", extra = "drop") %>% 
  distinct() %>% filter(stringr::str_detect(prefix, "K")) %>% filter(stringr::str_detect(prefix, "sigma"))

print(list_of_prefixes)

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
  
  edge <- sampled.coal[lower.tri(sampled.coal)]
  edge <- as.data.frame(edge)
  geoDist.edge <- as.matrix(geoDist)
  geoDist.edge <- geoDist.edge[lower.tri(geoDist.edge)]
  geoDist.edge <- as.data.frame(geoDist.edge)
  edge.dist <- cbind(edge, geoDist.edge)
  names(edge.dist)[1] <- "gen.dist"
  names(edge.dist)[2] <- "geo.dist"
  edge <- filter(edge.dist, geo.dist >= 20)
  names(edge)[1] <- "gen.dist"

  load(paste(prefix, "-est_wishart_out.Robj", sep=""), verbose=TRUE)
  col_pi <- rstan::extract(out$fit, "s", inc_warmup=TRUE, permute=FALSE)
  post_pi <- rstan::get_logposterior(out$fit,inc_warmup=FALSE)
  post_pi <- data.frame(matrix(unlist(post_pi), nrow=length(post_pi), byrow=FALSE))
  post_pi <- data.frame(matrix(unlist(post_pi), nrow=length(post_pi), byrow=FALSE))
  post_pi <- cbind(stack(post_pi[1:4]))
  names(post_pi)[1] <- "posterior"
  col_pi <- plyr::adply(col_pi, c(1,2,3))
  col_pi <- cbind(col_pi, post_pi) %>% 
    dplyr::mutate(col_pi = 1 - V1, iterations = as.numeric(iterations)) %>% 
    dplyr::select(-chains)
  
  col_pi <- col_pi %>% 
    #get mean post. for each chain
    dplyr::group_by(ind) %>% dplyr::mutate(mean.posterior = mean(posterior)) %>%
    #get max of those 4 means
    dplyr::ungroup() %>% dplyr::mutate(max.mean.posterior = max(mean.posterior)) %>% 
    as.data.frame() %>%
    #take diff of mean and max to ID which is actually the max
    dplyr::mutate(ismaxchain = ifelse(abs(mean.posterior - max.mean.posterior)==0, "max.post.chain","drop")) %>% 
    #keep just (250) values from max chain
    dplyr::filter(ismaxchain == "max.post.chain")
  
  max.div <- edge$gen.dist
  dist.max <- as.data.frame(max.div)
  names(dist.max)[1] <- "gen.dist"
  dist.max$value <- "true"
  est.col.pi <- as.data.frame(col_pi$col_pi)
  names(est.col.pi)[1] <- "gen.dist"
  est.col.pi$value <- "estimated"
  est.true.pi <- rbind(dist.max, est.col.pi)
  
    
}

write.table(est.true.pi, file=paste(prefix, "-est_edge_values.txt", sep=""))

    
}