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
workingdir = "~/Desktop/ZBH_stuff/compiled_output"
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
  tidyr::separate(., temp, into = c("garbage"), sep = "_", extra = "drop") %>% dplyr::select(-garbage) %>%  distinct()

#process each square iteration / size combo
wmModel_out <- data.frame("slurm_job_id"=NA, "slimIter"=NA, "sigma"=NA, "K"=NA,
                          "chain"=NA, "iteration"=NA,
                          "col_pi"=NA, "nbhd"=NA, "inDeme"=NA, "m"=NA, "posterior"=NA,
                          "theo_nbhd"=NA, "delta_nbhd"=NA, "Rousset_Nb"=NA, "MK_eff_nbhd"=NA,
                          "K_density"=NA, "eff_density"=NA, "theo_eff_nbhd"=NA, "MK_eff_sigma"=NA,
                          "MK_eff_dens_nbhd"=NA, "theo_eff_dens_nbhd"=NA, "eff_N"=NA, "true_s"=NA, "num_roots"=NA)

for (loop.iter in 1:length(list_of_prefixes$Robjfile)) {
  
  #define some variables
  Robjfile = list_of_prefixes[loop.iter, ]$Robjfile
  prefix = list_of_prefixes[loop.iter, ]$prefix
  #model_flavor = list_of_prefixes[loop.iter, ]$model_flavor
  
  print(paste0("starting file: ", Robjfile))
  
  labels <- prefix %>% as.data.frame() %>% dplyr::rename("prefix"=".") %>%
    tidyr::separate(., prefix, into = c("garbage","slurm_job_id","garbage2","slimIter","garbage3","sigma","garbage4","K"), sep = "_") %>%
    mutate(K = as.numeric(K), sigma = as.numeric(sigma)) %>% 
    dplyr::select(slurm_job_id, slimIter, sigma, K)
  
  #read dist files
  num_roots <- read.table(file=paste0(prefix, "-roots.txt"), header=FALSE)
  num_roots <- num_roots$V1
  K_density <- read.table(file=paste0(prefix, "_density"), header = FALSE)
  K_density <- K_density$V1
  MK_eff_sigma <- read.table(file=paste0(prefix, "-eff_sigma.txt"), header = FALSE)
  MK_eff_sigma <- MK_eff_sigma$V1
  eff_density <- read.table(file=paste0(prefix, "-eff_dens.txt"), header=FALSE)
  eff_density <- eff_density$V1
  eff_N <- read.table(file=paste0(prefix, "-eff_N.txt"), header = FALSE)
  eff_N <- eff_N$V1
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
  
  edge <- filter(dist.all, geo.dist >= 20)
  true_s <- mean(edge$gen.dist)
  
  #Rousset method for estimating Nb
  dist.fst$y <- dist.fst$fst / (1 - dist.fst$fst)
  reg <- lm(dist.fst$y ~ dist.fst$geo.dist)
  beta <- summary(reg)$coefficients[2, 1]
  Rousset_Nb <- 1 / beta
  
  #read WM_model output
  load(Robjfile, verbose=TRUE)
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
    mutate(theo_nbhd = 4*pi*(K_density)*(sigma^2)) %>% 
    mutate(delta_nbhd = theo_nbhd - nbhd)
  
  wmModel_pi <- wmModel_pi %>% dplyr::select(slurm_job_id, slimIter, sigma, K,
                                             chain, iteration,
                                             col_pi, nbhd, inDeme, m, posterior,
                                             theo_nbhd, delta_nbhd)
  
  wmModel_pi$Rousset_Nb <- Rousset_Nb
  wmModel_pi$MK_eff_nbhd <- 4*pi*(K_density)*MK_eff_sigma
  wmModel_pi$K_density <- K_density
  wmModel_pi$eff_density <- eff_density
  wmModel_pi$theo_eff_nbhd <- 4*pi*(wmModel_pi$K_density)*(wmModel_pi$sigma[1]*sqrt(3/2))^2
  wmModel_pi$MK_eff_sigma <- MK_eff_sigma
  wmModel_pi$MK_eff_dens_nbhd <- 4*pi*(wmModel_pi$eff_density)*MK_eff_sigma
  wmModel_pi$theo_eff_dens_nbhd <- 4*pi*(wmModel_pi$eff_density)*(wmModel_pi$sigma[1]*sqrt(3/2))^2
  wmModel_pi$eff_N <- eff_N
  wmModel_pi$true_s <- true_s
  wmModel_pi$num_roots <- num_roots
  wmModel_out <- rbind(wmModel_out, wmModel_pi)
  
  
}

wmModel_pi <- wmModel_out %>% filter(is.na(slurm_job_id)==FALSE)

####expcon
expcon_01 <- wmModel_pi
expcon_01$fraction <- 0.1
expcon_1 <- wmModel_pi
expcon_1$fraction <- 1
expcon_10 <- wmModel_pi
expcon_10$fraction <- 10

expcon_all <- rbind(expcon_01, expcon_1, expcon_10)
expcon_all$diff <- expcon_all$eff_N - expcon_all$nbhd
expcon_all$fraction <- factor(expcon_all$fraction, levels=c("0.1", "1", "10"))
expcon_all$K <- factor(expcon_all$K, levels=c("5", "10"), labels=c("expansion", "contraction"))
ggplot(expcon_all, aes(x=fraction, y=diff, fill=fraction)) + geom_boxplot() +
  facet_wrap(~K, labeller=label_parsed) + geom_hline(yintercept=0, linetype="dashed") +
  theme_bw() + theme(axis.text=element_text(size=12), axis.title=element_text(size=14),
                     strip.text=element_text(size=14), legend.position="none") + xlab("time (% of post-change N)") +
  ylab("true nbhd - estimated nbhd")


wmModel_pi$K <- factor(wmModel_pi$K, levels=c("2", "5", "10", "25"))

#save text output
write.table(wmModel_pi, file="all_sims_output.txt")

#end


ggplot(wmModel_pi, aes(x=eff_N, y=nbhd, color=K)) + geom_point() + 
  facet_wrap(~sigma, scales="free") + geom_abline(slope=1)


data_nbhd <- wmModel_pi %>% 
  group_by(sigma, K) %>% 
  summarise_at(vars(nbhd),
               list(min=min, Q1=~quantile(., probs = 0.05),
                    median=median, Q3=~quantile(., probs = 0.95),
                    max=max))


data_theo_nbhd <- wmModel_pi %>% 
  group_by(sigma, K) %>% 
  summarise_at(vars(eff_N),
               list(min=min, Q1=~quantile(., probs = 0.05),
                    median=median, Q3=~quantile(., probs = 0.95),
                    max=max))


nbhd_df <- as.data.frame(data_nbhd)
nbhd_df$est <- "nbhd"
theo_nbhd_df <- as.data.frame(data_theo_nbhd)
theo_nbhd_df$est <- "theo"
data_medians <- rbind(nbhd_df, theo_nbhd_df)

p <- ggplot(data_medians, aes(x=K, y=median)) + facet_wrap(~sigma, scales="free", labeller=label_parsed) + 
  geom_point(data=nbhd_df, size=3) + geom_errorbar(data=nbhd_df, aes(ymin=Q1, ymax=Q3, width=1.5)) +
  geom_point(data=theo_nbhd_df, color="#0072B2", size=3) + theme(legend.position="top", 
                                                                 axis.title = element_text(size=16), 
                                                                 axis.text=element_text(size=12), 
                                                                 strip.text=element_text(size=12), 
                                                                 legend.text=element_text(size=12)) +
  ylab("Neighborhood size") + theme_bw()

data_s <- wmModel_pi %>% 
  group_by(sigma, K) %>% 
  summarise_at(vars(col_pi),
               list(min=min, Q1=~quantile(., probs = 0.05),
                    median=median, Q3=~quantile(., probs = 0.95),
                    max=max))


data_true_s <- wmModel_pi %>% 
  group_by(sigma, K) %>% 
  summarise_at(vars(true_s),
               list(min=min, Q1=~quantile(., probs = 0.05),
                    median=median, Q3=~quantile(., probs = 0.95),
                    max=max))


s_df <- as.data.frame(data_s)
s_df$est <- "pi"
s_theo_df <- as.data.frame(data_true_s)
s_theo_df$est <- "true"
data_medians_s <- rbind(s_df, s_theo_df)

q <- ggplot(data_medians_s, aes(x=K, y=median)) + facet_wrap(~sigma, scales="free", labeller=label_parsed) + 
  geom_point(data=s_df, size=3) + geom_errorbar(data=s_df, aes(ymin=Q1, ymax=Q3, width=1.5)) +
  geom_point(data=s_theo_df, color="#0072B2", size=3) + theme(legend.position="top", 
                                                                 axis.title = element_text(size=16), 
                                                                 axis.text=element_text(size=12), 
                                                                 strip.text=element_text(size=12), 
                                                                 legend.text=element_text(size=12)) +
  ylab(expression(paste(pi[c]))) + theme_bw()


ggarrange(p, q, labels=c("A", "B"), ncol=1, nrow=2)

wmModel_pi$sigma <- factor(wmModel_pi$sigma, levels=c("0.5", "0.75", "1", "1.25", "1.5", 
                                                      "2"), 
                           labels=c("sigma==0.5", "sigma==0.75", "sigma==1", "sigma==1.25",
                                    "sigma==1.5", "sigma==2"))

ggplot(wmModel_pi, aes(x=theo_nbhd, y=eff_N, color=K)) + geom_point() + 
  facet_wrap(~sigma, scales="free", labeller=label_parsed) + geom_abline(slope=1) + 
  xlab(expression(paste(4*pi*(N[c]/w^2)*sigma^2))) + 
  ylab(expression(paste(4*pi*rho[eff]*sigma[eff]^2))) + theme_bw() +
  theme(axis.title = element_text(size=16), strip.text=element_text(size=12))

ggplot(wmModel_pi, aes(x=in_deme_rate, y=inDeme, color=K)) + geom_point() + facet_wrap(~sigma)

#####################
size_25_torus <- wmModel_pi
size_25_torus$size <- "25"
size_25_torus$range <- "torus"
size_50_torus <- wmModel_pi
size_50_torus$size <- "50"
size_50_torus$range <- "torus"
size_50_edges <- wmModel_pi
size_50_edges$size <- "50"
size_50_edges$range <- "edge"
size_25_edges <- wmModel_pi
size_25_edges$size <- "25"
size_25_edges$range <- "edge"

torus_runs <- rbind(size_25_torus, size_25_edges, size_50_edges, size_50_torus)
torus_runs$sigma <- factor(torus_runs$sigma, levels=c("0.5", "1", "2"))
torus_runs$diff_nbhd <- torus_runs$eff_N - torus_runs$nbhd
ggplot(torus_runs, aes(x=sigma, y=diff_nbhd, fill=range)) + geom_boxplot() + facet_wrap(~size)
ggplot(torus_runs, aes(x=sigma, y=diff_nbhd, fill=range)) + geom_boxplot() + ylim(-400,400) + facet_wrap(~size) +
  theme_bw() + theme(axis.title = element_text(size=14), axis.text=element_text(size=12),
                     strip.text = element_text(size=14), legend.text=element_text(size=12)) + xlab(expression(paste(sigma))) +
  ylab("true nbhd - estimated nbhd")

ggplot(torus_runs, aes(x=eff_N, y=nbhd, color=range)) + geom_point() +
  geom_abline(slope=1) + facet_wrap(~size, scales="free") + ylab("estimated nbhd") +
  xlab("true nbhd")

ggplot(torus_runs, aes(x=sigma, y=MK_eff_sigma, fill=range)) + geom_boxplot() + facet_wrap(~size, scales="free")

torus_runs$dens_eff <- torus_runs$eff_N/(4*pi*torus_runs$MK_eff_sigma)

ggplot(torus_runs, aes(x=K_density, y=dens_eff, color=range)) + geom_point() +
  facet_wrap(~size, scales="free") + xlab("census density") + ylab("effective density")

