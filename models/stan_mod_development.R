library(rstan)

source("wm_hom_cmpPar_mod_block_scaled_Gpr.R")

IBG <- stan_model(model_code=stanBlock)

source("wm_hom_cmpPar_mod_block_scaled.R")

IBD <- stan_model(model_code=stanBlock)


source('~/Dropbox/spatial_models/WM_model/wm_lib.R')

coords <- data.matrix(read.table("~/Dropbox/bedassle-paper/sims/slimulations/ibd/ibd1_coords.txt"))
pwp <- data.matrix(read.table("~/Dropbox/bedassle-paper/sims/slimulations/ibd/ibd1_pwp.txt"))
geoDist <- fields::rdist(coords)
hom <- 1-pwp
diag(hom) <- 1
k <- 1

# testing normal model
db <- list("N" = nrow(coords),
		   "L" = 1e4,
		   "hom" = hom,
		   "k" = k,
		   "geoDist" = geoDist)

fit <- runWM(stanMod=IBD,dataBlock=db,nChains=4,nIter=5e2,prefix="test",MLjumpstart=TRUE,nMLruns=3)


# testing G model
ut <- upper.tri(geoDist,diag=TRUE)
nGpar <- length(geoDist[ut])-length(which(geoDist[ut] < k))

idxsG <- getIdxsG(geoDist,k)

db <- list("N" = nrow(coords),
		   "L" = 1e4,
		   "hom" = hom,
		   "k" = k,
		   "geoDist" = geoDist,
		   "nGpar" = nGpar,
		   "idxsG" = idxsG)

fit <- runWM(stanMod=IBG,dataBlock=db,nChains=2,nIter=250,prefix="test",MLjumpstart=TRUE,nMLruns=3,Gmodel=TRUE)




if(FALSE){
	
tmp <- ml2init(db=db,mod=IBG,nRuns=5,Gmodel=TRUE,prefix="testG_mod")

library(doParallel)
library(foreach)

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

foreach(i =1:4, .packages="rstan") %dopar% {
fit <- runWM(stanMod=IBG,dataBlock=db,nChains=1,nIter=2e4,prefix=sprintf("test%s",i))
}


#fit2 <- runWM(stanMod=mod,dataBlock=db,nChains=2,nIter=1e3,prefix="testMLjumpstart",MLjumpstart=TRUE,nMLruns=20)

parallel::stopCluster(cl)

generateInitPars <- function(dataBlock,breakLimit=1e4,nChains,prefix){
#	recover()
	scl_min <- min(dataBlock$hom)
	scl_max <- max(dataBlock$hom-scl_min)
	posdef <- FALSE
	counter <- 0
	inPrBounds <- FALSE
	while((!posdef | !inPrBounds) & counter < breakLimit){
		k <- dataBlock$k
		s <- rbeta(1,0.8,0.2)
			while(s==1){
				s <- rbeta(1,0.8,0.2)			
			}
		logm <- runif(1,min=-30,max=0)
		m <- exp(logm)
		nbhd <- abs(rnorm(1,1,10))
		inDeme <- rbeta(1,0.9,0.5)
		nugget <- abs(rnorm(1,0.05,0.01))
		loginDeme <- log(inDeme)
		lognugget <- log(nugget)
		Gvec <- db$geoDist[db$idxsG]
		G <- db$geoDist
		nbhd <- abs(rnorm(1,1,10))
		parHom <- makeParaHom(s,m,k,nbhd,inDeme,nugget,G)
		parHom <- (parHom-scl_min)/scl_max
		posdef <- all(eigen(parHom)$values > 0)
		inPrBounds <- checkPrBounds(nbhd,logm,loginDeme,lognugget)
		counter <- counter + 1
	}
	if(counter == breakLimit){
		stop("\nunable to generate initial parameters that generate a positive-definite covariance matrix\n")
	}
	if(posdef){
		cat("\ninitial parameters produce a positive definite matrix\n")
	}
	initPars <- list("s"=s,
				 	 "logm"=logm,
				 	 "Gvec"=Gvec,
				 	 "nbhd"=nbhd,
				 	 "loginDeme"=loginDeme,
				 	 "lognugget"=lognugget,
				 	 "parHom" = parHom)
	save(initPars,file=paste0(prefix,"_initPars.Robj"))
	return(initPars)
}


	tv <- unlist(lapply(1:1219,function(i){var(extract(out$fit,sprintf("Gvec[%s]",i),inc_warmup=FALSE,permute=FALSE)[,1,1])}))
plot(db$geoDist[db$idxsG],tv)

extractGamPars <- function(fit,N){
	A <- matrix(NA,N,N)
	B <- matrix(NA,N,N)
	for(i in 1:N){
		for(j in 1:N){
			x <- sprintf("A[%s,%s]",i,j)
			A[i,j] <- rstan::extract(fit,pars=x,inc_warmup=FALSE,permuted=FALSE)[1]
			x <- sprintf("B[%s,%s]",i,j)
			B[i,j] <- rstan::extract(fit,pars=x,inc_warmup=FALSE,permuted=FALSE)[1]
		}
	}
	return(
		list("A" = A,
			 "B" = B))
}

extractG <- function(fit,chainNo,N){
	G <- array(NA,dim=c(fit@sim$n_save[chainNo],N,N))
	for(i in 1:N){
		for(j in 1:N){
			x <- sprintf("G[%s,%s]",i,j)
			G[,i,j] <- rstan::extract(fit,pars=x,inc_warmup=FALSE,permuted=FALSE)[,chainNo,1]
		}
	}
	return(G)
}

parG <- extractG(out$fit,1,db$N)
ut <- upper.tri(db$geoDist)

plot(db$geoDist,db$geoDist,type='n')
	invisible(
		lapply(sample(1:250,30),
			function(i){
				points(db$geoDist[ut],parG[i,,][ut],pch=19,col=adjustcolor(1,0.1))
			}))




gp <- extractGamPars(out$fit,db$N)

plot(db$geoDist,gp$A/gp$B^2)


for(i in 1:(db$N-1)){
	for(j in (i+1):db$N){
		cat(dgamma(db$geoDist[i,j],shape=gp$A[i,j],rate=gp$B[i,j],log=TRUE),"\n")
	}
}
	
getGammaPars <- function(d,maxD){
	ordMag <- floor(log(maxD,base=10))
	slope <- max(c(-0.1*ordMag,-0.9))
	a <- 0.1*d^(3+slope)
	b <- 0.1*d^(2+slope)
	return(a/b^2)
#	return(
#		list("a" = a,
#				"b" = b))
}


plot(0,xlim=c(0,maxD),ylim=c(0,10),type='n')

maxD <- 2e8
testD <- seq(1,maxD,length.out=1000)
	lines(testD,getGammaPars(testD,maxD=maxD),col="blue")




par(mfrow=c(1,2))
g1 <- rgamma(1e3,shape=5,scale=1)
g2 <- rgamma(1e3,shape=50,scale=10)
hist(g1)
	abline(v=mean(g1),col="red")
hist(g2)
	abline(v=mean(g2),col="red")


maxD <- 2e2
testD <- seq(1,maxD,length.out=1000)
getGammaPars <- function(d,maxD){
	a <- (d/maxD)^(2.1)
	b <- (d/maxD)^(1.1)
	return(a/b)
}
plot(testD,getGammaPars(testD,maxD=maxD),col="blue")
}







