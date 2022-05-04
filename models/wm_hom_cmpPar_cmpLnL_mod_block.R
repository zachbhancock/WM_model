stanBlock <- "
functions {
	vector probHom(int lut, real k, real s, real m, real nbhd, real inDeme, vector geoDist) {
		vector[lut] pHom;
		real pIBD;
		for(i in 1:lut){
			pIBD = geoDist[i] > k ? modified_bessel_second_kind(0, m * geoDist[i]) / nbhd : inDeme;
			pHom[i] = pIBD + (1-pIBD)*s;
		}
		return pHom;	
	}
}
data {
	int<lower=0> lut;				// length of the upper triangle, excluding the diagonal
	real<lower=0> k;					// distance threshold within which Wright-Malecot breaks down
	vector[lut] obsHom; 			// observed pairwise homozygosity (upper triangle of pairwise matrix)
	vector[lut] se; 					// standard errors of pairwise homozygosity measures
	vector[lut] geoDist; 			// upper triangle of matrix of pairwise geographic distance 
}
transformed data {
	vector[lut] obsHomScl;
	real sclMn;
	real sclMx;
	vector[lut]seScl;
	sclMn = min(obsHom);
	sclMx = max(obsHom-sclMn);
	obsHomScl = (obsHom - sclMn)/sclMx;
	seScl = se/sclMx;
}
parameters {
	real<lower=0> nbhd;				// Wright's neighborhood size
	real<lower=0> m;					// scaled migration rate
	real<lower=0,upper=1> inDeme;	// within deme p(IBD)
	real<lower=0,upper=1> s;		// minimum rate of IBD
}
transformed parameters {
	vector[lut] pHom;
	vector[lut] pHomScl;
	pHom = probHom(lut, k, s, m, nbhd, inDeme, geoDist);
	pHomScl = (pHom-sclMn)/sclMx;
}
model {
	s ~ beta(1,0.1);					// prior on minimum relatedness
	m ~ normal(0,0.1);				// prior on scaled migration rate
	nbhd ~ normal(10,100);			// prior on neighborhood size
	inDeme ~ beta(1,0.01);			// prior on within-deme p(IBD)
//	for(i in 1:lut) obsHomScl[i] ~ normal(pHomScl[i],seScl[i]);
	obsHomScl ~ normal(pHomScl,seScl);
}
"
