################################################################
################################################################
# wright-malecot IBD model fitting function library
################################################################
################################################################

################################
# running the model
################################

runWM <- function(stanMod,dataBlock,nChains,nIter,prefix){
	initPars <- lapply(1:nChains,function(i){generateInitPars(dataBlock=dataBlock)})
	fit <- sampling(object = stanMod,
				 data = dataBlock,
				 iter = nIter,
				 chains = nChains,
				 save_warmup = FALSE,
				 init=initPars,
				 thin = ifelse(nIter/500 > 
        					   1, 
            				   floor(nIter/500), 1))
	out <- list("dataBlock" = dataBlock,
				"fit" = fit)
    saveOut(fit=fit,outPrefix=prefix)
	save(out,file=paste0(prefix,"_out.Robj"))
	vizWMout(wmOutfile=paste0(prefix,"_out.Robj"),outPrefix=prefix)
}

runWM_cmpLnl <- function(stanMod,dataBlock,nChains,nIter,prefix){
	message("running ML analyses to generate initial parameter estimates")
	initPars <- lapply(1:nChains,function(i){ml2init(db=dataBlock,mod=stanMod,nRuns=5e2)})
	message("running Bayesian analyses")
	fit <- sampling(object = stanMod,
				 data = dataBlock,
				 iter = nIter,
				 chains = nChains,
				 save_warmup = FALSE,
				 init=initPars,
				 thin = ifelse(nIter/500 > 
        					   1, 
            				   floor(nIter/500), 1))
	out <- list("dataBlock" = dataBlock,
				"fit" = fit)
    saveOut(fit=fit,outPrefix=prefix)
	save(out,file=paste0(prefix,"_out.Robj"))
	vizWMout_cmpLnl(wmOutfile=paste0(prefix,"_out.Robj"),outPrefix=prefix)
}

ml2init <- function(db,mod,nRuns){
	mlRuns <- lapply(1:nRuns,
					function(i){
						optimizing(object=mod,data=db)
					})
	bestRun <- which.max(unlist(lapply(mlRuns,"[[","value")))
	inits <- list("nbhd" = as.numeric(mlRuns[[bestRun]]$par[which(names(mlRuns[[bestRun]]$par)=="nbhd")]),
				  "inDeme" = as.numeric(mlRuns[[bestRun]]$par[which(names(mlRuns[[bestRun]]$par)=="inDeme")]),
				  "s" = as.numeric(mlRuns[[bestRun]]$par[which(names(mlRuns[[bestRun]]$par)=="s")]),
				  "m" = as.numeric(mlRuns[[bestRun]]$par[which(names(mlRuns[[bestRun]]$par)=="m")]))
	return(inits)
}

makeParaHom <- function(s,m,k,nbhd,inDeme,nugget,geoDist){
	pIBD <- besselK(nu=0,x=sqrt(m) * geoDist) / nbhd
	pIBD[which(geoDist < k)] <- inDeme
	pHom <- pIBD + (1-pIBD)*s
	diag(pHom) <- diag(pHom) + nugget
	return(pHom)
}

generateInitPars <- function(dataBlock,breakLimit=1e4,nChains){
	scl_min <- min(dataBlock$hom)
	scl_max <- max(dataBlock$hom-scl_min)
	posdef <- FALSE
	counter <- 0
	while(!posdef & counter < breakLimit){
		k <- dataBlock$k
		s <- rbeta(1,0.8,0.2)
			while(s==1){
				s <- rbeta(1,0.8,0.2)			
			}
		m <- abs(rnorm(1,0.1,0.1))
		nbhd <- abs(rnorm(1,1,10))
		inDeme <- rbeta(1,0.9,0.5)
		nugget <- abs(rnorm(1,0.05,0.01))
		parHom <- makeParaHom(s,m,k,nbhd,inDeme,nugget,dataBlock$geoDist)
		parHom <- (parHom-scl_min)/scl_max
		posdef <- all(eigen(parHom)$values > 0)
		counter <- counter + 1
	}
	if(counter == breakLimit){
		stop("\nunable to generate initial parameters that generate a positive-definite covariance matrix\n")
	}
	initPars <- list("s"=s,
				 	 "m"=m,
				 	 "nbhd"=nbhd,
				 	 "inDeme"=inDeme,
				 	 "nugget"=nugget)
	save(initPars,file="initPars.Robj")
	return(initPars)
}

################################
# processing output
################################

getPhom <- function(model.fit,chain.no,N){
	par.cov <- array(NA,dim=c(model.fit@sim$n_save[chain.no],N,N))
	for(i in 1:N){
		for(j in 1:N){
			my.par <- sprintf("pHom[%s,%s]",i,j)
			par.cov[,i,j] <- rstan::extract(model.fit,pars=my.par,inc_warmup=FALSE,permuted=FALSE)[,chain.no,]
		}
	}
	return(par.cov)
}

saveOut <- function(fit,outPrefix){
	s <- rstan::extract(fit,"s",inc_warmup=FALSE)
	if("nugget" %in% names(fit)){
		nugget <- rstan::extract(fit,"nugget",inc_warmup=FALSE)
	}
	m <- rstan::extract(fit,"m",inc_warmup=FALSE)
	nbhd <- rstan::extract(fit,"nbhd",inc_warmup=FALSE)
	inDeme <- rstan::extract(fit,"inDeme",inc_warmup=FALSE)
	ptEsts <- list("s" = mean(s[[1]]),
				   "m" = mean(m[[1]]),
				   "nbhd" = mean(nbhd[[1]]),
				   "inDeme" = mean(inDeme[[1]]))
	postDist <- list("s" = s,
		 			 "nugget" = nugget,
					 "m" = m,
					 "nbhd" = nbhd,
					 "inDeme" = inDeme)
	if("nugget" %in% names(fit)){
		ptEsts[["nugget"]] = mean(nugget[[1]])
		postDist[["nugget"]] = nugget[[1]]
	}
	outPars <- list("pt" = ptEsts,
				  "post" = postDist)
	save(outPars,file=paste0(outPrefix,"_pars.Robj"))
	return(invisible("saved"))
}

################################
# visualizing output
################################


vizWMout <- function(wmOutfile,outPrefix){
	load(wmOutfile)
	post <- rstan::get_logposterior(out$fit,inc_warmup=FALSE)
	s <- rstan::extract(out$fit,"s",inc_warmup=FALSE,permute=FALSE)
	nugget <- rstan::extract(out$fit,"nugget",inc_warmup=FALSE,permute=FALSE)
	m <- rstan::extract(out$fit,"m",inc_warmup=FALSE,permute=FALSE)
	nbhd <- rstan::extract(out$fit,"nbhd",inc_warmup=FALSE,permute=FALSE)
	inDeme <- rstan::extract(out$fit,"inDeme",inc_warmup=FALSE,permute=FALSE)
	pHom <- invisible(lapply(1:length(post),
				function(i){
					getPhom(out$fit,i,out$dataBlock$N)}))
	nChains <- length(post)
	chainCols <- c("blue","goldenrod1","red","forestgreen","purple","black")[1:nChains]
	pdf(file=paste0(outPrefix,"plots.pdf"),width=12,heigh=8)
		makeCmpParPlots(post,m,nbhd,s,nugget,inDeme,chainCols)
		for(i in 1:length(post)){
			par(mfrow=c(1,1))
				plotFit(out,pHom[[i]],chainCols[i])
		}
	dev.off()
}

vizWMout_cmpLnl <- function(wmOutfile,outPrefix){
	load(wmOutfile)
	post <- rstan::get_logposterior(out$fit,inc_warmup=FALSE)
	s <- rstan::extract(out$fit,"s",inc_warmup=FALSE,permute=FALSE)
	m <- rstan::extract(out$fit,"m",inc_warmup=FALSE,permute=FALSE)
	nbhd <- rstan::extract(out$fit,"nbhd",inc_warmup=FALSE,permute=FALSE)
	inDeme <- rstan::extract(out$fit,"inDeme",inc_warmup=FALSE,permute=FALSE)
	pHom <- rstan::extract(out$fit,"pHom",inc_warmup=FALSE,permute=FALSE)
	nChains <- length(post)
	chainCols <- c("blue","goldenrod1","red","forestgreen","purple","black")[1:nChains]
	pdf(file=paste0(outPrefix,"plots.pdf"),width=12,heigh=8)
		makeCmpParPlots(post,m,nbhd,s,nugget=NULL,inDeme,chainCols)
		for(i in 1:length(post)){
			par(mfrow=c(1,1))
				plotFit_cmpLnl(out,pHom[,i,],chainCols[i])
		}
	dev.off()
}

makeCmpParPlots <- function(post,m,nbhd,s,nugget=NULL,inDeme,chainCols){
	par(mfrow=c(2,3))
		matplot(Reduce("cbind",post),ylab="posterior probability",type='l',lwd=1.5,lty=1,col=chainCols)
		matplot(m[,,1],type='l',lwd=1.5,lty=1,col=chainCols,xlab="",ylab="m")
		matplot(nbhd[,,1],type='l',lwd=1.5,lty=1,col=chainCols,xlab="",ylab="nbhd")
		matplot(s[,,1],type='l',lwd=1.5,lty=1,col=chainCols,xlab="",ylab="s")
		matplot(inDeme[,,1],type='l',lwd=1.5,lty=1,col=chainCols,xlab="sampled mcmc iterations",ylab="inDeme")
		if(!is.null(nugget)){
			matplot(nugget[,,1],type='l',lwd=1.5,lty=1,col=chainCols,xlab="",ylab="nugget")
		}
	for(i in 1:length(chainCols)){
		par(mfrow=c(2,5))
			plot(m[,i,1],nbhd[,i,1],pch=20,col=adjustcolor(chainCols[i],0.5),xlab="m",ylab="nbhd")
			plot(m[,i,1],s[,i,1],pch=20,col=adjustcolor(chainCols[i],0.5),xlab="m",ylab="s")
			plot(m[,i,1],inDeme[,i,1],pch=20,col=adjustcolor(chainCols[i],0.5),xlab="m",ylab="inDeme")
			plot(nbhd[,i,1],s[,i,1],pch=20,col=adjustcolor(chainCols[i],0.5),xlab="nbhd",ylab="s")
			plot(nbhd[,i,1],inDeme[,i,1],pch=20,col=adjustcolor(chainCols[i],0.5),xlab="nbhd",ylab="inDeme")
			plot(s[,i,1],inDeme[,i,1],pch=20,col=adjustcolor(chainCols[i],0.5),xlab="s",ylab="inDeme")
		if(!is.null(nugget)){
			plot(m[,i,1],nugget[,i,1],pch=20,col=adjustcolor(chainCols[i],0.5),xlab="m",ylab="nugget")
			plot(nbhd[,i,1],nugget[,i,1],pch=20,col=adjustcolor(chainCols[i],0.5),xlab="nbhd",ylab="nugget")
			plot(s[,i,1],nugget[,i,1],pch=20,col=adjustcolor(chainCols[i],0.5),xlab="s",ylab="nugget")
			plot(inDeme[,i,1],nugget[,i,1],pch=20,col=adjustcolor(chainCols[i],0.5),xlab="inDeme",ylab="nugget")
		}
	}
}

plotFit <- function(out,pHom,chainCol){
	plot(out$dataBlock$geoDist,out$myData$hom,
			ylim=range(c(pHom,out$dataBlock$hom)),type='n',
			xlab="pairwise distance",ylab="pairwise homozygosity")
		invisible(
			lapply(seq(1,250,length.out=25),function(i){
				points(out$dataBlock$geoDist,pHom[i,,],pch=20,col=adjustcolor(chainCol,0.05))
			}))
	points(out$dataBlock$geoDist,out$dataBlock$hom)
}

plotFit_cmpLnl <- function(out,pHom,chainCol){
	plot(out$dataBlock$geoDist,out$myData$hom,
			ylim=range(c(pHom,out$dataBlock$hom)),xlim=range(out$dataBlock$geoDist),type='n',
			xlab="pairwise distance",ylab="pairwise homozygosity")
		invisible(
			lapply(seq(1,250,length.out=25),function(i){
				points(out$dataBlock$geoDist,pHom[i,],pch=20,col=adjustcolor(chainCol,0.05))
			}))
	points(out$dataBlock$geoDist,out$dataBlock$hom)
}