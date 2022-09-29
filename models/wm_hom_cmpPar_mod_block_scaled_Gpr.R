stanBlock <- "
functions {
	matrix probHom(int N, real k, real s, real m, real nbhd, real inDeme, real nugget, matrix G) {
		matrix[N,N] pHom;
		matrix[N,N] nugMat;
		real pIBD;
		pHom = rep_matrix(0,N,N);
		nugMat = diag_matrix(rep_vector(nugget,N));
		for(i in 1:N){
			for(j in i:N){
				if(G[i,j] > k){
					pIBD = modified_bessel_second_kind(0, sqrt(m) * G[i,j]) / nbhd;
					pHom[i,j] = pIBD + (1-pIBD)*s;
				} else {
					pHom[i,j] = inDeme + (1-inDeme)*s;
				}
				pHom[j,i] = pHom[i,j];
			}
		}
		return pHom + nugMat;	
	}
	matrix makeG(int N, int nGpar, matrix geoDist, int[,] idxsG, real[] Gvec){
		matrix[N,N] G = geoDist;
		for(i in 1:nGpar){
			G[idxsG[i,1],idxsG[i,2]] = Gvec[i];
			G[idxsG[i,2],idxsG[i,1]] = Gvec[i];
		}
		return G;
	}
	real[] getGammaPars(int nGpar, matrix geoDist, int[,] idxsG, int b){
		real d;
		real gammaPar[nGpar];
		real e = 2.1;
		if(b){
			e = 1.1;
		}
		for(i in 1:nGpar){
			d = geoDist[idxsG[i,1],idxsG[i,2]];
			gammaPar[i] = d^e;
		}
		return gammaPar;
	}
}
data {
	int<lower=2> N; 	  					// number of samples
	int<lower=N+1> L;	    			// number of independent loci
	cov_matrix[N] hom; 				// observed pairwise homozygosity
	real<lower=0> k;					// distance threshold within which Wright-Malecot breaks down
	matrix[N, N] geoDist; 				// matrix of pairwise geographic distance
	int<lower=N> nGpar;				// number of entries in the upper triangle (including the diagonal) of the geoDist matrix that are greater than k
	int idxsG[nGpar,2];			// indices of the entries in the upper triangle (including the diagonal) of the geoDist matrix that are greater than k
}
transformed data {
	real<lower=0> scl_min = min(hom);
	real<lower=0> scl_max = max(hom-scl_min);
	matrix[N,N] hom_scl = (hom-scl_min)/scl_max;				// n.loci multiplied by the sample covariance
	matrix[N,N] Lhom_scl  = L * hom_scl;				// n.loci multiplied by the sample covariance
	real a[nGpar] = getGammaPars(nGpar, geoDist, idxsG, 0);
	real b[nGpar] = getGammaPars(nGpar, geoDist, idxsG, 1);
}
parameters {
	real<lower=0> nbhd;							// Wright's neighborhood size
	real<lower=-25,upper=0> logm;				// scaled migration rate
	real<lower=-5,upper=0> loginDeme;			// within deme p(IBD)
	real<lower=0,upper=1> s;				// minimum IBD
	real<lower=-25,upper=-1> lognugget;			// nugget
	real<lower=0> Gvec[nGpar];
}
transformed parameters {
	matrix<lower=0>[N, N] G = makeG(N, nGpar, geoDist, idxsG, Gvec);					// pairwise effective distance
	real m = exp(logm);
	real inDeme = exp(loginDeme);
	real nugget = exp(lognugget);
	matrix[N,N] pHom = probHom(N, k, s, m, nbhd, inDeme, nugget, G);					// probability of being homozygous
	matrix[N,N] pHom_scl = (pHom-scl_min)/scl_max;											// probability of being homozygous
}
model {
	for(i in 1:nGpar){
		Gvec[i] ~ gamma(a[i],b[i]);
	}
	s ~ beta(1,0.1);					// prior on minimum relatedness
	logm ~ normal(-5,1);					// prior on scaled migration rate
	nbhd ~ normal(100,1000);				// prior on neighborhood size
	loginDeme ~ normal(0,1);				// prior on within-deme p(IBD)
	lognugget ~ normal(0,1);				// prior on nugget
	Lhom_scl ~ wishart(L, pHom_scl);			// likelihood function
}
"