stanBlock <- "
functions {
	matrix probHom(int N, real k, real s, real m, real nbhd, real inDeme, real nugget, matrix geoDist) {
		matrix[N,N] pHom = rep_matrix(0,N,N);
		matrix[N,N] nugMat = diag_matrix(rep_vector(nugget,N));
		real pIBD;
		for(i in 1:N){
			for(j in i:N){
				if(geoDist[i,j] > k){
					pIBD = modified_bessel_second_kind(0, sqrt(m) * geoDist[i,j]) / nbhd;
					pHom[i,j] = pIBD + (1-pIBD)*s;
				} else {
					pHom[i,j] = inDeme + (1-inDeme)*s;
				}
				pHom[j,i] = pHom[i,j];
			}
		}
		return pHom + nugMat;	
	}
}
data {
	int<lower=2> N; 	  					// number of samples
	int<lower=N+1> L;	    			// number of independent loci
	cov_matrix[N] hom; 				// observed pairwise homozygosity
	real<lower=0> k;					// distance threshold within which Wright-Malecot breaks down
	matrix[N, N] geoDist; 				// matrix of pairwise geographic distance
}
transformed data {
	real<lower=0> scl_min = min(hom);
	real<lower=0> scl_max = max(hom-scl_min);
	matrix[N,N] hom_scl = (hom-scl_min)/scl_max;				// n.loci multiplied by the sample covariance
	matrix[N,N] Lhom_scl = L * hom_scl;				// n.loci multiplied by the sample covariance
}
parameters {
	real<lower=0> nbhd;				// Wright's neighborhood size
	real<lower=0> m;				// scaled migration rate
	real<lower=0> inDeme;			// within deme p(IBD)
	real<lower=0,upper=1> s;				// minimum rate of IBD
	real<lower=0> nugget;			// nugget
}
transformed parameters {
	matrix[N,N] pHom = probHom(N, k, s, m, nbhd, inDeme, nugget, geoDist);					// probability of being homozygous
	matrix[N,N] pHom_scl = (pHom-scl_min)/scl_max;					// probability of being homozygous
}
model {
	s ~ beta(1,0.1);					// prior on minimum relatedness
	m ~ beta(0.000001,1);					// prior on scaled migration rate
	nbhd ~ normal(100,1000);				// prior on neighborhood size
	inDeme ~ beta(1,0.01);				// prior on within-deme p(IBD)
	nugget ~ normal(0,1);				// prior on nugget
	Lhom_scl ~ wishart(L, pHom_scl);			// likelihood function
}
"