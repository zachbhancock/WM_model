#idea: prep empirical data from NCBI and dryad to use in WM model


#load libraries
library(vcfR)
library(dplyr)
library(tidyr)
library(ggplot2)

rm(list = ls())
gc()


#get SRA metadata table saved from NCBI ------------
metadf <- read.delim("/Users/rachel/WM_model/empirical_data/bioprj_PRJNA473221_Bombus-vosnesenskii/SraRunTable.txt", sep  = ",")
metadf <- metadf %>% filter(Organism == "Bombus vosnesenskii")
metadf <- metadf %>% dplyr::rename("run_acc_sra"="Run",
                             "read_type_sra"="LibraryLayout") %>% 
  mutate(organism = gsub(" ","-",Organism)) %>% 
  mutate(link = paste0(BioProject,"_",organism)) %>%
  arrange(run_acc_sra) %>% 
  #mutate(sampid_assigned_for_bioinf = paste0("sample", 0:(nrow(.)-1))) %>% 
  mutate(sampid_assigned_for_bioinf = Sample.Name) %>% 
  mutate(pop = "dummy") %>% 
  mutate(identifiers_biosamp = "null", infraspecies_biosamp = "null") %>% 
  mutate(total_number_reads_sra = Bases/AvgSpotLen) %>% 
  mutate(total_reads_written_to_final_fastq = total_number_reads_sra)
key <- metadf %>% dplyr::select(link, run_acc_sra, sampid_assigned_for_bioinf, 
                             pop, identifiers_biosamp, infraspecies_biosamp,
                             read_type_sra, total_number_reads_sra, 
                             total_reads_written_to_final_fastq)
key %>% arrange(total_reads_written_to_final_fastq) %>% dplyr::mutate(order = 1:n()) %>% 
  ggplot() + 
  geom_point(aes(x=total_reads_written_to_final_fastq, y=order)) +
  geom_vline(xintercept = 100000, colour = "orange")


#generate and save samplerenamekey and popmap ----------
readr::write_delim(key,
                   paste0("empirical_data/bioprj_PRJNA473221_Bombus-vosnesenskii/samplenamekey.txt"),
                   delim = "\t")
popmap <- key %>% dplyr::select(sampid_assigned_for_bioinf, pop)
readr::write_delim(popmap,
                   paste0("empirical_data/bioprj_PRJNA473221_Bombus-vosnesenskii/popmap"),
                   delim = "\t")


#calculate and save pairwise geographic distances ------------
geodist <- metadf %>% dplyr::select(BioProject, run_acc_sra, organism, Lat_Lon)
geodist <- geodist %>% separate(Lat_Lon, into = c("lat","latdir","lon","londir"), sep = " ", remove = FALSE) %>% 
  mutate(lat = as.numeric(lat), lon = as.numeric(lon))
geodist %>% group_by(latdir,londir) %>% dplyr::summarise(n=n())
#S and W are negative
geodist <- geodist %>% mutate(lon = ifelse(londir=="W",lon*(-1),lon),
                              lat = ifelse(latdir=="S",lat*(-1),lat))
genetic <- geodist %>% dplyr::select(lon, lat, run_acc_sra) %>% dplyr::rename("y"="lon", "x"="lat")

#calculate pairwise great circle distances
#returns lower triangle distance matrix in km btw points from decimal degrees long,lat df
#will work if longitudes are bewteen -180 to 180 or 0 to 360
pwGCD <- function(points, quantile) {
  
  #calc quantile values
  lower = ((100 - quantile)/2)/100
  upper = (100 - (lower*100))/100
  #calc GCD
  #quantile(obj_of_class_dist) gives same result as long list of pairwise dists with self-self 0's excluded
  gcd <- fields::rdist.earth(points,miles=FALSE)
  rownames(gcd) <- points$run_acc_sra
  colnames(gcd) <- points$run_acc_sra 
  #find the max for X quantile
  max.quant.gcd = stats::quantile(gcd, c(lower, upper)) %>% .[2]
  #select the first pair that represent the X% quantile max dist so we can plot it on map
  max.quant.gcdpr <- gcd %>% as.data.frame() %>% mutate(samp = row.names(.)) %>% 
    dplyr::select(samp,everything()) %>% 
    pivot_longer(., 2:ncol(.), names_to = "samp2", values_to = "dist_km") %>% 
    mutate(distfrommax = abs(dist_km - max.quant.gcd)) %>% 
    filter(distfrommax == min(distfrommax)) %>% 
    .[1,]
  samp1 <- points %>% filter(run_acc_sra == max.quant.gcdpr$samp)
  samp2 <- points %>% filter(run_acc_sra == max.quant.gcdpr$samp2)
  max.quant.gcdpr  <- rbind(samp1,samp2)
  #save
  gcd_dists <- list("max.quant.gcd" = max.quant.gcd,
                    "pw.gcd" = gcd,
                    "max.quant.gcdpr" = max.quant.gcdpr)
  return(gcd_dists)
  
}

quantile = 100
gcd.genetic <- pwGCD(points = genetic, quantile = quantile)
names(gcd.genetic)[names(gcd.genetic) == "max.quant.gcd"] <- paste("max.",quantile,".gcd",".genetic", sep="")
names(gcd.genetic)[names(gcd.genetic) == "pw.gcd"] <- paste("pw.gcd",".genetic", sep="")
names(gcd.genetic)[names(gcd.genetic) == "max.quant.gcdpr"] <- paste("max.",quantile,".gcdpr",".genetic", sep="")

max_and_pw_dists <- c(gcd.genetic)

save(max_and_pw_dists, file=paste("empirical_data/bioprj_PRJNA473221_Bombus-vosnesenskii/max_and_pw_distsbioprj_PRJNA473221_Bombus-vosnesenskii.Robj",sep = ""))


#get VCF file from dyrad ----------
vcf <- vcfR::read.vcfR("empirical_data/bioprj_PRJNA473221_Bombus-vosnesenskii/Samtools_vos_popsover5.vcf")
#extract genotypes
gt <- vcfR::extract.gt(vcf)
# transpose and make into data.frame
#NA here means SNP not called
gt <- as.data.frame(t(gt), stringsAsFactors = FALSE) %>% as.matrix()
#check which syntax is used to denote genotypes
genos <- if(any(grepl("/",gt))){
    c("0/0","0/1","1/0","1/1")
  } else if(any(grepl("|",gt))){
    c("0\\|0","0\\|1","1\\|0","1\\|1")
}
gt <- gsub(genos[1],0,gt)	
gt <- gsub(genos[2],1,gt)
gt <- gsub(genos[3],1,gt)
gt <- gsub(genos[4],2,gt)
class(gt) <- "numeric"




#calc pwp matrix -----------------

#this code from Gideon from divdiv

#' Calculate pairwise pi between two diploid individuals
#'
#' @param ind1 A genotype vector for a single diploid individual
#'			 for which the \emph{i}th element 
#'			 gives the frequency of the counted 
#'			 allele in that individual at the \emph{i}th locus.
#' @param ind2 A genotype vector for a single diploid individual
#'			 for which the \emph{i}th element 
#'			 gives the frequency of the counted 
#'			 allele in that individual at the \emph{i}th locus.
#' @param L The number of loci genotyped in \emph{both} individuals 1 
#'			and 2.
#' @return An estimate of pairwise pi between individuals 1 and 2.
#' @details The genotype data should consist of {0,0.5,1}
#'			with missing data indicated with \code{NA}.
#'			If the genotype vectors \code{ind1} and \code{ind2} 
#'			consist only of polymorphic loci, and \code{L} is the 
#'			number of co-genotyped polymorphic loci, 
#'			this function will return pairwise pi at 
#'			polymorphic sites, which may not be 
#'			comparable across datasets.
#' @export
calcPWP <- function(ind1,ind2,L){
  nMD <- ifelse((any(is.na(ind1)) | any(is.na(ind2))),
                length(unique(which(is.na(ind1)),which(is.na(ind2)))),
                0)
  if(is.null(L)){
    L <- length(ind1) - nMD
  }
  diff.homs = sum(abs(ind1-ind2)==1,na.rm=TRUE)
  hets = length(
    unique(
      c(which(ind1==0.5 & !is.na(ind2)),
        which(ind2==0.5 & !is.na(ind1)))
    )
  )
  return((diff.homs + hets/2)/L)
}

#' Calculate pairwise pi between all pairs of individuals in a dataset
#'
#' @param freqs  A genotype matrix with \emph{N} rows and 
#'			 \emph{L} columns, where \emph{N} is the number 
#'			 of individuals and \emph{L} is the number of 
#'			 loci, for which the \emph{i},\emph{j}th element 
#'			 gives the frequency of the counted allele in each 
#'			 individual at each locus.
#' @param coGeno A symmetric matrix (dimensions \emph{N} x \emph{N}) 
#'			for which the \emph{i},\emph{j}th element gives 
#'			the number of base pairs genotyped in \emph{both} 
#'			samples \emph{i} and \emph{j}. Default is \code{NULL}.
#'			If left unspecified, each cell will be assigned a value of 
#'			\code{L} (i.e., indicating no missing data).
#' @param quiet A \code{logical} value that controls whether 
#'			a progress par is displayed. Default is \code{FALSE}.
#' @return An estimate of pairwise pi between all individuals 
#'			in the data matrix.
#' @details The genotype data should consist of {0,0.5,1}
#'			with missing data indicated with \code{NA}.
#'			If the genotype matrix consists only of 
#'			polymorphic loci, and \code{coGeno} either 
#'			contains the number of co-genotyped polymorphic 
#'			loci for each sample pair or is not specified 
#'			this function will return pairwise pi at 
#'			polymorphic sites, which may not be 
#'			comparable across datasets.
#' @export
freqs2pairwisePi <- function(freqs,coGeno=NULL,quiet=FALSE){
  if(is.null(coGeno)){
    coGeno <- matrix(ncol(freqs),nrow(freqs),nrow(freqs))
  }
  n <- nrow(freqs)
  pwp <- matrix(NA,n,n)
  if(!quiet){
    prog <- utils::txtProgressBar(min=0,max=n+(n*(n-1))/2,char="*",style=3)
  }
  for(i in 1:n){
    for(j in i:n){
      if(!quiet){
        utils::setTxtProgressBar(prog,(i-1)*n+j-(i*(i-1)/2))
      }
      pwp[i,j] <- calcPWP(freqs[i,],freqs[j,],coGeno[i,j])
      if(i == j){
        pwp[i,i] <- 2 * pwp[i,i]
      } else if(i != j){
        pwp[j,i] <- pwp[i,j]
      }
    }
  }
  row.names(pwp) <- row.names(coGeno)
  colnames(pwp) <- row.names(coGeno)
  return(pwp)
}

pwp <- freqs2pairwisePi(freqs=gt/2,quiet=FALSE)
globalPi <- mean(pwp[upper.tri(pwp,diag=TRUE)])

popgenstats <- list("thetaW" = NULL,
                    "pwp" = pwp,
                    "globalPi" = globalPi,
                    "pcs" = NULL,
                    "het" = NULL)

save(popgenstats,file=paste0("empirical_data/bioprj_PRJNA473221_Bombus-vosnesenskii/popgenstats.0.5.bioprj_PRJNA473221_Bombus-vosnesenskii_stacks_fromdryad_popgenstats.Robj"))

Nloci=ncol(gt)
BPstats <- list("nLoci" = Nloci)

save(BPstats,file=paste0("empirical_data/bioprj_PRJNA473221_Bombus-vosnesenskii/bpstats.0.5.bioprj_PRJNA473221_Bombus-vosnesenskii_stacks_fromdryad_BPstats.Robj"))


