#Load Stan outputs and make master .csv of outputs

# load libraries and source relevant functions
library(ggplot2, warn.conflicts = FALSE, quietly = TRUE)
library(plyr, warn.conflicts = FALSE, quietly = TRUE)
library(dplyr, warn.conflicts = FALSE, quietly = TRUE)
library(tidyr, warn.conflicts = FALSE, quietly = TRUE)


#rm(list=ls())
#gc()

#print session info and args for HPCC metadata/log files
print("STARTING TO RUN collate_all_txt_outputs.R SCRIPT")
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

#print objects loaded in R (so we can see functions and variables were correctly loaded)
print("----------------------------------------------------------------------------------------")
print("----------------------------------------------------------------------------------------")
print("so far we have loaded libraries and printed args")
print("objects currently loaded in R:")
ls()

#local testing
#list_of_txtfiles <- list.files(path = "txt_outputs", pattern = "true_est_values.txt", full.names = TRUE)

#get list of all square iteration / square size combos to process
list_of_txtfiles <- list.files(path = workingdir, pattern = "-est.txt", full.names = TRUE)

#process each square iteration / size combo
#make blank df to store results of loop in
df <- data.frame("slurm_job_id"=NA, "slimIter"=NA, "sigma"=NA, "K"=NA,
                          "model_flavor"=NA, "chain"=NA, "iteration"=NA,
                          "col_pi"=NA, "nbhd"=NA, "inDeme"=NA, "m"=NA, "posterior"=NA,
                          "theo_nbhd"=NA, "delta_nbhd"=NA)

for (loop.iter in 1:length(list_of_txtfiles)) {
  
  print(paste0("on loop iteration ", loop.iter, " of ", length(list_of_txtfiles)))  
  file <- list_of_txtfiles[loop.iter]
  out.sub <- read.table(file=paste0(file), header=TRUE)
  df <- rbind(df, out.sub)

}
#remove empty first row of df that we set up above
df <- df %>% dplyr::filter(is.na(slurm_job_id)==F)
write.table(df, file=paste(workingdir,"alloutputs-",Sys.Date(),".txt", sep=""))



# !!! START HERE !!! #

#summarize to mean and SD for 250 estimated values from chain with best posterior, 1000 true values, and then keep one single best posterior estimate
df.summary <- df %>% group_by(slurm_job_id, slimIter, sigma, model_flavor, value) %>% 
  dplyr::summarise(mean.pi = mean(gen.dist), sd.pi = sd(gen.dist)) %>% 
  ungroup()

#do a few quick checks
# df.summary %>% distinct(slurm_job_id) %>% nrow() #should = N sigma values
# df.summary %>% distinct(slimIter) %>% nrow() #should = N slim iterations
# df.summary %>% distinct(sigma) %>% nrow() #should = N sigma values

# df.summary %>% group_by(slimIter,sigma) %>% summarise(n.values = n()) %>% ungroup() %>% distinct(n.values) #should be 3
# df.summary %>% dplyr::select(slimIter,sigma,squareRep) %>% distinct() %>% 
#   group_by(slimIter,sigma) %>% summarise(n=n()) %>% group_by(sigma) %>% summarise(n=n()) #should be N sigma values long and n should = N slim iters

#save full compiled output
write.table(df.summary, file=paste0(workingdir,"summary_alloutputs-",Sys.Date(),".txt"))


#graveyard ----------
