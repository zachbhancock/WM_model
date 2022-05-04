################################################################
################################################################
#	testing DivDiv stan models
################################################################
################################################################

################################
# load libraries and define 
#	useful functions
################################
library(rstan)

# get pairwise homozygosity from 
#	Wright-Malecot model
getWMhom <- function(s,m,k,nbhd,inDeme,geoDist){
	pIBD <- besselK(nu=0,x=sqrt(m) * geoDist) / nbhd
	pIBD[which(geoDist < k)] <- inDeme
	pHom <- pIBD + (1-pIBD)*s
	return(pHom)
}

ml2init <- function(db,mod,nRuns){
	mlRuns <- lapply(1:nRuns,
					function(i){
						optimizing(object=mod,data=db)
					})
	bestRun <- which.max(unlist(lapply(mlRuns,"[[","value")))
	# lnls <- unlist(lapply(mlRuns,"[[","value"))
	# nbhds <- unlist(lapply(mlRuns,function(x){x$par[1]}))
	# ms <- unlist(lapply(mlRuns,function(x){x$par[2]}))
	# inDemes <- unlist(lapply(mlRuns,function(x){x$par[3]}))
	# ss <- unlist(lapply(mlRuns,function(x){x$par[4]}))
	inits <- list("nbhd" = as.numeric(mlRuns[[bestRun]]$par[which(names(mlRuns[[bestRun]]$par)=="nbhd")]),
				  "inDeme" = as.numeric(mlRuns[[bestRun]]$par[which(names(mlRuns[[bestRun]]$par)=="inDeme")]),
				  "s" = as.numeric(mlRuns[[bestRun]]$par[which(names(mlRuns[[bestRun]]$par)=="s")]),
				  "m" = as.numeric(mlRuns[[bestRun]]$par[which(names(mlRuns[[bestRun]]$par)=="m")]))
	return(inits)
}

################################
# compile rstan model
################################
source("../../models/wm_hom_cmpPar_cmpLnL_mod_block.R")
wm <- stan_model(model_code=stanBlock)


################################
# simulate dataset
################################

N <- 40
coords <- matrix(runif(2*N,min=0,50),nrow=N,ncol=2)
geoDist <- fields::rdist(coords)
ut <- upper.tri(geoDist,diag=FALSE)

bad <- TRUE
while(bad){
	pars <- list("s" = max(0.95,rbeta(1,0.8,0.2)),
				 "k" = 1,
				 "m" = abs(rnorm(1,0.1,0.1)),
				 "nbhd" = abs(rnorm(1,1,10)),
				 "inDeme" = rbeta(1,0.9,0.5))
	obsHom <- getWMhom(pars$s,pars$m,pars$k,pars$nbhd,pars$inDeme,geoDist)
	if(any(obsHom > 1) | diff(range(obsHom)) > 0.2){
		bad <- TRUE
	} else {
		bad <- FALSE
	}
}
#plot(geoDist[ut],obsHom[ut])
pars[["N"]] <- N
pars[["coords"]] <- coords
pars[["geoDist"]] <- geoDist


db <- list("lut" = N*(N-1)/2,
		   "k" = 1,
		   "obsHom" = obsHom[ut],
		   "se" = abs(rnorm(N*(N-1)/2,mean=0,sd=0.05)),
		   "geoDist" = geoDist[ut])

nChains <- 2
inits <- lapply(1:nChains,function(i){ml2init(db=db,mod=wm,nRuns=1e3)})

fit <- sampling(object = wm,
				data = db,
				iter = 3e3,
				chains = 2,
				init = inits,
				save_warmup=FALSE)

save(db,file="db.Robj")
save(fit,file="fit.Robj")
save(pars,file="pars.Robj")

################################
# visualize output
################################

if(FALSE){
s <- rstan::extract(fit,"s",permute=FALSE)[,,1]
m <- rstan::extract(fit,"m",permute=FALSE)[,,1]
nbhd <- rstan::extract(fit,"nbhd",permute=FALSE)[,,1]
inDeme <- rstan::extract(fit,"inDeme",permute=FALSE)[,,1]
pHom <- rstan::extract(fit,"pHom",permute=FALSE)[,1,]
lpd <- rstan::get_logposterior(fit)

pdf(file="sim_test.pdf",width=10,height=8)
plot(db$geoDist,db$obsHom,col="red",pch=19,ylim=range(pHom))
	invisible(
		lapply(sample(1:nrow(pHom),25),
			function(i){
				points(db$geoDist,pHom[i,],pch=20,col=adjustcolor(1,0.1))
			}))
plot(lpd[[1]],type='n',ylim=range(lpd),xlim=c(0,length(lpd[[1]])),xlab="iterations",ylab="lpd")
	lines(lpd[[1]],col="blue")
	lines(lpd[[2]],col="goldenrod")
par(mfrow=c(2,2))
	matplot(s,type='l',ylim=range(c(s,pars$s)),col=c("blue","goldenrod"),lty=1)
		abline(h=pars$s,col="red")
	matplot(m,type='l',ylim=range(c(m,pars$m)),col=c("blue","goldenrod"),lty=1)
		abline(h=pars$m,col="red")
	matplot(nbhd,type='l',ylim=range(c(nbhd,pars$nbhd)),col=c("blue","goldenrod"),lty=1)
		abline(h=pars$nbhd,col="red")
	matplot(inDeme,type='l',ylim=range(c(inDeme,pars$inDeme)),col=c("blue","goldenrod"),lty=1)
		abline(h=pars$inDeme,col="red")
dev.off()
}



