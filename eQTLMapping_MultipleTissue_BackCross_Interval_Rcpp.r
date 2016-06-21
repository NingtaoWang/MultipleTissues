REPEAT_LIMIT <<- 1001
RELATIVE_DIFF <<- 0.000001
Method_H0 <<- 1
Method_H1TissueInteraction<<- 2
Method_H1ClusterInteraction <<- 3
Method_H1Equal <<- 4

library(mvtnorm)
library(abind)
library(Rcpp)
library(RcppArmadillo)


cppFunction('
Rcpp::NumericVector MyOuterMvnorm( 	Rcpp::NumericVector DataR,
						Rcpp::NumericVector MuR,
						Rcpp::NumericVector SigmaR) {
	Rcpp::IntegerVector DimsDataR = DataR.attr("dim");
	arma::cube DataC(DataR.begin(), DimsDataR[0], DimsDataR[1], DimsDataR[2]);
	int TC = DimsDataR[2] - 1;

	Rcpp::IntegerVector DimsMuR = MuR.attr("dim");
	arma::cube MuC(MuR.begin(), DimsMuR[0], DimsMuR[1], DimsMuR[2]);

	Rcpp::IntegerVector DimsSigmaR = SigmaR.attr("dim");
	arma::cube SigmaC(SigmaR.begin(), DimsSigmaR[0], DimsSigmaR[1], DimsSigmaR[2]);
	int T2C = DimsSigmaR[2] - 1;

	int n_iter = DimsDataR[0] * DimsDataR[1] * DimsMuR[0] * DimsMuR[1], 
	n_iter_ijk = DimsDataR[0] * DimsDataR[1] * DimsMuR[0], 
	n_iter_ij = DimsDataR[0] * DimsDataR[1], 
	n_iter_i = DimsDataR[0];
	Rcpp::NumericVector Output(n_iter);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		arma::vec DataVec = DataC.subcube(i, j, 0, i, j, TC);
		arma::vec MuVec = MuC.subcube(k, l, 0, k, l, TC);
		arma::vec SigmaVec = SigmaC.subcube(k, l, 0, k, l, T2C);
		arma::mat SigmaMat(SigmaVec.begin(), DimsMuR[2], DimsMuR[2]);
		arma::vec DataMuVec = DataVec - MuVec;
		double tmp1 = pow(2*PI, 0.5 * DimsDataR[2]);
		double tmp2 = pow(arma::det(SigmaMat), 0.5);
		Output[h] = exp(arma::as_scalar(- arma::trans(DataMuVec) * arma::inv(SigmaMat) * DataMuVec / 2)) / (tmp1 * tmp2);
	}
	return Output;      
}', depends=c("RcppArmadillo"))

cppFunction('
Rcpp::NumericVector MyOuterPsiPsi( 	Rcpp::NumericMatrix Psiik,
						Rcpp::NumericVector Psiijkl) {
	Rcpp::IntegerVector DimsPsiijkl = Psiijkl.attr("dim");
	
	int n_iter = DimsPsiijkl[0] * DimsPsiijkl[1] * DimsPsiijkl[2] * DimsPsiijkl[3], 
	n_iter_ijk = DimsPsiijkl[0] * DimsPsiijkl[1] * DimsPsiijkl[2], 
	n_iter_ij = DimsPsiijkl[0] * DimsPsiijkl[1], 
	n_iter_i = DimsPsiijkl[0];
	Rcpp::NumericVector Output(n_iter);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		Output[h] = Psiik(i, k) * Psiijkl(h);
	}
	return Output;      
}', depends=c("RcppArmadillo"))

cppFunction('
Rcpp::NumericVector MyOuterOmegaOmega( 	Rcpp::NumericMatrix Omegajl,
							Rcpp::NumericVector Omegaijkl) {
	Rcpp::IntegerVector DimsOmegaijkl = Omegaijkl.attr("dim");
	
	int n_iter = DimsOmegaijkl[0] * DimsOmegaijkl[1] * DimsOmegaijkl[2] * DimsOmegaijkl[3], 
	n_iter_ijk = DimsOmegaijkl[0] * DimsOmegaijkl[1] * DimsOmegaijkl[2], 
	n_iter_ij = DimsOmegaijkl[0] * DimsOmegaijkl[1], 
	n_iter_i = DimsOmegaijkl[0];
	Rcpp::NumericVector Output(n_iter);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		Output[h] = Omegajl(j, l) * Omegaijkl(h);
	}
	return Output;      
}', depends=c("RcppArmadillo"))

cppFunction('
arma::mat MyOuterMu( 	Rcpp::NumericVector DataR,
				Rcpp::NumericVector PsiPsiOmegaOmega) {
	Rcpp::IntegerVector DimsDataR = DataR.attr("dim");
	arma::cube DataC(DataR.begin(), DimsDataR[0], DimsDataR[1], DimsDataR[2]);

	Rcpp::IntegerVector DimsPsiPsiOmegaOmega = PsiPsiOmegaOmega.attr("dim");
	
	int n_iter = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1] * DimsPsiPsiOmegaOmega[2] * DimsPsiPsiOmegaOmega[3], 
	n_iter_ijk = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1] * DimsPsiPsiOmegaOmega[2], 
	n_iter_ij = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1], 
	n_iter_i = DimsPsiPsiOmegaOmega[0],
	n_iter_t = DimsDataR[2];
	arma::mat Output(n_iter, n_iter_t);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		for(int t = 0; t < n_iter_t; t++){
			Output(h,t) = DataC(i,j,t) * PsiPsiOmegaOmega(h);
		}
	}
	return Output;      
}', depends=c("RcppArmadillo"))

cppFunction('
arma::mat MyOuterSigma( Rcpp::NumericVector DataR,
				Rcpp::NumericVector MuR,
				Rcpp::NumericVector PsiPsiOmegaOmega) {
	Rcpp::IntegerVector DimsDataR = DataR.attr("dim");
	arma::cube DataC(DataR.begin(), DimsDataR[0], DimsDataR[1], DimsDataR[2]);

	Rcpp::IntegerVector DimsMuR = MuR.attr("dim");
	arma::cube MuC(MuR.begin(), DimsMuR[0], DimsMuR[1], DimsMuR[2]);

	Rcpp::IntegerVector DimsPsiPsiOmegaOmega = PsiPsiOmegaOmega.attr("dim");
	
	int n_iter = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1] * DimsPsiPsiOmegaOmega[2] * DimsPsiPsiOmegaOmega[3], 
	n_iter_ijk = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1] * DimsPsiPsiOmegaOmega[2], 
	n_iter_ij = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1], 
	n_iter_i = DimsPsiPsiOmegaOmega[0],
	n_iter_ts = DimsDataR[2] * DimsDataR[2],
	n_iter_t = DimsDataR[2];
	arma::mat Output(n_iter, n_iter_ts);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		for(int ts = 0; ts < n_iter_ts; ts++){
			int t = floor(ts / n_iter_t);
			int s = ts - t * n_iter_t;
			Output(h,ts) = PsiPsiOmegaOmega(h) * (DataC(i, j, t) - MuC(k, l, t)) * (DataC(i, j, s) - MuC(k, l, s));
		}
	}
	return Output;      
}', depends=c("RcppArmadillo"))

cppFunction('
arma::mat MyOuterPsiPsiOmegaOmegaSigmaInverse(		Rcpp::NumericVector SigmaR,
							Rcpp::NumericVector PsiPsiOmegaOmega,
							int T) {
	Rcpp::IntegerVector DimsSigmaR = SigmaR.attr("dim");
	arma::cube SigmaC(SigmaR.begin(), DimsSigmaR[0], DimsSigmaR[1], DimsSigmaR[2]);
	int T2C = DimsSigmaR[2] - 1;

	Rcpp::IntegerVector DimsPsiPsiOmegaOmega = PsiPsiOmegaOmega.attr("dim");
	
	int n_iter = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1] * DimsPsiPsiOmegaOmega[2] * DimsPsiPsiOmegaOmega[3], 
	n_iter_ijk = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1] * DimsPsiPsiOmegaOmega[2], 
	n_iter_ij = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1], 
	n_iter_i = DimsPsiPsiOmegaOmega[0],
	n_iter_ts = T * T,
	n_iter_t = T;
	arma::mat Output(n_iter, n_iter_ts);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		arma::vec SigmaVec = SigmaC.subcube(k, l, 0, k, l, T2C);
		arma::mat SigmaMat(SigmaVec.begin(), T, T);
		arma::mat SigmaInverse = arma::inv(SigmaMat);

		for(int ts = 0; ts < n_iter_ts; ts++){
			int t = floor(ts / n_iter_t);
			int s = ts - t * n_iter_t;
			Output(h,ts) = PsiPsiOmegaOmega(h) * SigmaInverse(t, s);
		}
	}
	return Output;      
}', depends=c("RcppArmadillo"))

cppFunction('
arma::mat MyOuterPsiPsiOmegaOmegaSigmaInverseData(	Rcpp::NumericVector DataR,
							Rcpp::NumericVector SigmaR,
							Rcpp::NumericVector PsiPsiOmegaOmega) {
	Rcpp::IntegerVector DimsDataR = DataR.attr("dim");
	arma::cube DataC(DataR.begin(), DimsDataR[0], DimsDataR[1], DimsDataR[2]);
	int TC = DimsDataR[2] - 1;

	Rcpp::IntegerVector DimsSigmaR = SigmaR.attr("dim");
	arma::cube SigmaC(SigmaR.begin(), DimsSigmaR[0], DimsSigmaR[1], DimsSigmaR[2]);
	int T2C = DimsSigmaR[2] - 1;

	Rcpp::IntegerVector DimsPsiPsiOmegaOmega = PsiPsiOmegaOmega.attr("dim");
	
	int n_iter = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1] * DimsPsiPsiOmegaOmega[2] * DimsPsiPsiOmegaOmega[3], 
	n_iter_ijk = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1] * DimsPsiPsiOmegaOmega[2], 
	n_iter_ij = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1], 
	n_iter_i = DimsPsiPsiOmegaOmega[0],
	n_iter_t = DimsDataR[2];
	arma::mat Output(n_iter, n_iter_t);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		arma::vec DataVec = DataC.subcube(i, j, 0, i, j, TC);
		arma::vec SigmaVec = SigmaC.subcube(k, l, 0, k, l, T2C);
		arma::mat SigmaMat(SigmaVec.begin(), DimsDataR[2], DimsDataR[2]);
		arma::vec SigmaInverseDataVec = arma::inv(SigmaMat) * DataVec;

		for(int t = 0; t < n_iter_t; t++){
			Output(h,t) = PsiPsiOmegaOmega(h) * SigmaInverseDataVec(t);
		}
	}
	return Output;      
}', depends=c("RcppArmadillo"))

cppFunction('
arma::mat MyOuterPsiPsiOmegaOmegaSigmaInverseDataMu(	Rcpp::NumericVector DataR,
							Rcpp::NumericVector MuR,
							Rcpp::NumericVector SigmaR,
							Rcpp::NumericVector PsiPsiOmegaOmega) {
	Rcpp::IntegerVector DimsDataR = DataR.attr("dim");
	arma::cube DataC(DataR.begin(), DimsDataR[0], DimsDataR[1], DimsDataR[2]);
	int TC = DimsDataR[2] - 1;

	Rcpp::IntegerVector DimsMuR = MuR.attr("dim");
	arma::cube MuC(MuR.begin(), DimsMuR[0], DimsMuR[1], DimsMuR[2]);

	Rcpp::IntegerVector DimsSigmaR = SigmaR.attr("dim");
	arma::cube SigmaC(SigmaR.begin(), DimsSigmaR[0], DimsSigmaR[1], DimsSigmaR[2]);
	int T2C = DimsSigmaR[2] - 1;

	Rcpp::IntegerVector DimsPsiPsiOmegaOmega = PsiPsiOmegaOmega.attr("dim");
	
	int n_iter = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1] * DimsPsiPsiOmegaOmega[2] * DimsPsiPsiOmegaOmega[3], 
	n_iter_ijk = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1] * DimsPsiPsiOmegaOmega[2], 
	n_iter_ij = DimsPsiPsiOmegaOmega[0] * DimsPsiPsiOmegaOmega[1], 
	n_iter_i = DimsPsiPsiOmegaOmega[0],
	n_iter_t = DimsDataR[2];
	arma::mat Output(n_iter, n_iter_t);
	for(int h = 0; h < n_iter; h++){
		int l = floor(h / n_iter_ijk);
		int ijk = h - l * n_iter_ijk;
		int k = floor(ijk / n_iter_ij);
		int ij = ijk - k * n_iter_ij;
		int j = floor(ij / n_iter_i);
		int i = ij - j * n_iter_i;

		arma::vec DataVec = DataC.subcube(i, j, 0, i, j, TC);
		arma::vec MuVec = MuC.subcube(k, l, 0, k, l, TC);
		arma::vec DataMuVec = DataVec - MuVec;
		arma::vec SigmaVec = SigmaC.subcube(k, l, 0, k, l, T2C);
		arma::mat SigmaMat(SigmaVec.begin(), DimsDataR[2], DimsDataR[2]);
		arma::vec SigmaInverseDataMuVec = arma::inv(SigmaMat) * DataMuVec;

		for(int t = 0; t < n_iter_t; t++){
			Output(h,t) = PsiPsiOmegaOmega(h) * SigmaInverseDataMuVec(t);
		}
	}
	return Output;      
}', depends=c("RcppArmadillo"))

Maximize.H0 <- function(data, PsiPsiOmegaOmega, PsiOmega.kl, K, L, n, m, T)
{
	PsiOmegaData <- MyOuterMu(data, PsiPsiOmegaOmega)
	PsiOmegaData <- array(PsiOmegaData, dim=c(n, m, K, L, T))
	PsiOmegaData.klt <- apply(PsiOmegaData, 3:5, sum)
	new.mu <- apply(PsiOmegaData.klt, 3, FUN = function(x){x / PsiOmega.kl})
	new.mu <- array(new.mu, dim = c(K, L, T))

	return(new.mu)
}

Newton.Equal <- function(data, sigma, PsiPsiOmegaOmega, K, L, n, m, T)
{
	sigma.trans <- array(sigma, dim = c(K, L, T^2))
	new.mu <- array(NA, dim = c(K, L, T))

	PsiPsiOmegaOmega.SigmaInverse <- MyOuterPsiPsiOmegaOmegaSigmaInverse(sigma.trans, PsiPsiOmegaOmega, T)
	PsiPsiOmegaOmega.SigmaInverse <- array(PsiPsiOmegaOmega.SigmaInverse, dim = c(n, m, K, L, T, T))
	PsiPsiOmegaOmega.SigmaInverse.ltt <- apply(PsiPsiOmegaOmega.SigmaInverse, 4:6, sum)

	PsiPsiOmegaOmega.SigmaInverse.data <- MyOuterPsiPsiOmegaOmegaSigmaInverseData(data, sigma.trans, PsiPsiOmegaOmega)
	PsiPsiOmegaOmega.SigmaInverse.data <- array(PsiPsiOmegaOmega.SigmaInverse.data, dim = c(n, m, K, L, T))
	PsiPsiOmegaOmega.SigmaInverse.data.lt <- apply(PsiPsiOmegaOmega.SigmaInverse.data, 4:5, sum)

	for(l in 1:L)
	{
		new.mu[1, l, ] <- solve(PsiPsiOmegaOmega.SigmaInverse.ltt[l, , ]) %*% PsiPsiOmegaOmega.SigmaInverse.data.lt[l, ]
		new.mu[2, l, ] <- new.mu[1, l, ]
	}

	return(new.mu)
}

Newton.TissueInteraction <- function(data, mu, sigma, PsiPsiOmegaOmega, K, L, n, m, T)
{
	cl <- apply(mu[1, , ] - mu[2, , ], 1, mean)
	mu[2, , ] <- mu[1, , ] - cl %o% rep(1,T)

	sigma.trans <- array(sigma, dim = c(K, L, T^2))
	new.mu <- array(NA, dim = c(K, L, T))

	PsiPsiOmegaOmega.SigmaInverse <- MyOuterPsiPsiOmegaOmegaSigmaInverse(sigma.trans, PsiPsiOmegaOmega, T)
	PsiPsiOmegaOmega.SigmaInverse <- array(PsiPsiOmegaOmega.SigmaInverse, dim = c(n, m, K, L, T, T))
	PsiPsiOmegaOmega.SigmaInverse.ltt <- apply(PsiPsiOmegaOmega.SigmaInverse, 4:6, sum)
	PsiPsiOmegaOmega.SigmaInverse.tl <- apply(PsiPsiOmegaOmega.SigmaInverse.ltt, 1, FUN=function(x){-diag(x)})
	PsiPsiOmegaOmega.SigmaInverse.lt <- t(PsiPsiOmegaOmega.SigmaInverse.tl)

	PsiPsiOmegaOmega.SigmaInverse.data.mu <- MyOuterPsiPsiOmegaOmegaSigmaInverseDataMu(data, mu, sigma.trans, PsiPsiOmegaOmega)
	PsiPsiOmegaOmega.SigmaInverse.data.mu <- array(PsiPsiOmegaOmega.SigmaInverse.data.mu, dim = c(n, m, K, L, T))
	PsiPsiOmegaOmega.SigmaInverse.data.mu.lt <- apply(PsiPsiOmegaOmega.SigmaInverse.data.mu, 4:5, sum)

	dQ.dcl.1st <- -apply(PsiPsiOmegaOmega.SigmaInverse.data.mu[, , 2, , ], 3, sum)
	dQ.dcl.2nd <- -apply(PsiPsiOmegaOmega.SigmaInverse[, , 2, , , ], 3, sum)

	new.mu[1, , ] <- mu[1, , ] - PsiPsiOmegaOmega.SigmaInverse.data.mu.lt / PsiPsiOmegaOmega.SigmaInverse.lt
	new.mu[2, , ] <- new.mu[1, , ] + mu[2, , ] - mu[1, , ] + {dQ.dcl.1st / dQ.dcl.2nd} %o% rep(1,T)

	return(new.mu)
}

Newton.ClusterInteraction <- function(data, mu, sigma, PsiPsiOmegaOmega, K, L, n, m, T)
{
	c <- apply(mu[1, , ] - mu[2, , ], 2, mean)
	mu[2, , ] <- mu[1, , ] - rep(1,L) %o% c

	sigma.trans <- array(sigma, dim = c(K, L, T^2))
	new.mu <- array(NA, dim = c(K, L, T))

	PsiPsiOmegaOmega.SigmaInverse <- MyOuterPsiPsiOmegaOmegaSigmaInverse(sigma.trans, PsiPsiOmegaOmega, T)
	PsiPsiOmegaOmega.SigmaInverse <- array(PsiPsiOmegaOmega.SigmaInverse, dim = c(n, m, K, L, T, T))
	PsiPsiOmegaOmega.SigmaInverse.ltt <- apply(PsiPsiOmegaOmega.SigmaInverse, 4:6, sum)
	PsiPsiOmegaOmega.SigmaInverse.tl <- apply(PsiPsiOmegaOmega.SigmaInverse.ltt, 1, FUN=function(x){-diag(x)})
	PsiPsiOmegaOmega.SigmaInverse.lt <- t(PsiPsiOmegaOmega.SigmaInverse.tl)

	PsiPsiOmegaOmega.SigmaInverse.data.mu <- MyOuterPsiPsiOmegaOmegaSigmaInverseDataMu(data, mu, sigma.trans, PsiPsiOmegaOmega)
	PsiPsiOmegaOmega.SigmaInverse.data.mu <- array(PsiPsiOmegaOmega.SigmaInverse.data.mu, dim = c(n, m, K, L, T))
	PsiPsiOmegaOmega.SigmaInverse.data.mu.lt <- apply(PsiPsiOmegaOmega.SigmaInverse.data.mu, 4:5, sum)

	dQ.dc.1st <- -apply(PsiPsiOmegaOmega.SigmaInverse.data.mu[, , 2, , ], 4, sum)
	dQ.dc.2nd <- -apply(PsiPsiOmegaOmega.SigmaInverse[, , 2, , , ], 4:5, sum)
	dQ.dc.2nd.diag <- diag(dQ.dc.2nd)

	new.mu[1, , ] <- mu[1, , ] - PsiPsiOmegaOmega.SigmaInverse.data.mu.lt / PsiPsiOmegaOmega.SigmaInverse.lt
	new.mu[2, , ] <- new.mu[1, , ] + mu[2, , ] - mu[1, , ] + rep(1,L) %o% {dQ.dc.1st / dQ.dc.2nd.diag}

	return(new.mu)
}

eQTL.StepE <- function(data, mu, sigma, p, q, K, L, n, m)
{
	sigma.trans <- array(sigma, dim = c(K, L, T^2))
	Fijkl <- MyOuterMvnorm(data, mu, sigma.trans)
	Fijkl <- array(Fijkl, dim = c(n, m, K, L))
	Fijkl[is.nan(Fijkl)] <- 1e-300
	Fijkl[Fijkl < 1e-300] <- 1e-300

	Psi.ijkl <- apply(Fijkl,1:3,FUN=function(x){x*q/(x%*%q)})
	Psi.ijkl <- aperm(Psi.ijkl, c(2,3,4,1))
	F1.ijk <- apply(Fijkl,1:3,FUN=function(x){x%*%q})
	F1.log.ijk <- log(F1.ijk)
	F1.log.ik <- apply(F1.log.ijk,c(1,3),sum)
	F1.log.ik.adjtmp <- apply(F1.log.ik,1,
		FUN=function(x){
			tmpx <- max(x)-1400
			x[which(x<tmpx)]<-tmpx
			return(x)})	
	F1.log.ik.adj <- apply(F1.log.ik.adjtmp,2,FUN=function(x){x-max(x)/2-min(x)/2})
	F1.log.ik.adj <- t(F1.log.ik.adj)
	F1.ik <- exp(F1.log.ik.adj)*p
	F1.ik <- apply(F1.ik,1,FUN=function(x){x/sum(x)})
	Psi.ik <- t(F1.ik)

	Omega.ijkl <- apply(Fijkl,c(2,4),FUN=function(x){x*p/rowSums(x*p)})
	Omega.ijkl <- array(Omega.ijkl, dim=c(n,K,m,L))
	Omega.ijkl <- aperm(Omega.ijkl, c(1,3,2,4))
	F2.ijl <- apply(Fijkl,c(2,4),FUN=function(x){rowSums(x*p)})
	F2.log.ijl <- log(F2.ijl)
	F2.log.jl <- apply(F2.log.ijl,c(2,3),sum)
	F2.log.jl.adjtmp <- apply(F2.log.jl,1,
		FUN=function(x){
			tmpx <- max(x)-1400
			x[which(x<tmpx)]<-tmpx
			return(x)})	
	F2.log.jl.adj <- apply(F2.log.jl.adjtmp,2,FUN=function(x){x-max(x)/2-min(x)/2})
	F2.log.jl.adj <- t(F2.log.jl.adj)
	F2.jl <- apply(F2.log.jl.adj,1,FUN=function(x){exp(x)*q/(exp(x)%*%q)})
	Omega.jl <- t(F2.jl)

	F1.log.i.mean <- apply(F1.log.ik.adjtmp,2,FUN=function(x){max(x)/2+min(x)/2})
	F1.ik.tmp <- exp(F1.log.ik.adj)*p
	F1.i.tmp <- apply(F1.ik.tmp,1,sum)
	F1.i.log <- log(F1.i.tmp)
	F2.log.j.mean <- apply(F2.log.jl.adjtmp,2,FUN=function(x){max(x)/2+min(x)/2})
	F2.jl.tmp <- apply(F2.log.jl.adj,1,FUN=function(x){exp(x)*q})
	F2.j.tmp <- apply(F2.jl.tmp,2,sum)
	F2.j.log <- log(F2.j.tmp)
	LogLikelihood <- sum(F1.i.log)+sum(F2.j.log)+sum(F1.log.i.mean)+sum(F2.log.j.mean)

	return(list(Psi.ik = Psi.ik, Psi.ijkl = Psi.ijkl, 
			Omega.jl = Omega.jl, Omega.ijkl = Omega.ijkl, LogLikelihood = LogLikelihood))
}

eQTL.StepM <- function(data, mu, sigma, Psi.ik, Psi.ijkl, Omega.jl, Omega.ijkl, K, L, n, m, T, method)
{
	PsiPsi <- MyOuterPsiPsi(Psi.ik, Psi.ijkl)
	PsiPsi <- array(PsiPsi, dim = c(n, m, K, L))
	PsiPsi.jl <- apply(PsiPsi,c(2,4),sum)

	OmegaOmega <- MyOuterOmegaOmega(Omega.jl, Omega.ijkl)
	OmegaOmega <- array(OmegaOmega, dim = c(n, m, K, L))
	OmegaOmega.ik <- apply(OmegaOmega,c(1,3),sum)

	PsiOmega.kl <- apply(OmegaOmega+PsiPsi, 3:4, sum)

	PsiPsiOmegaOmega <- PsiPsi + OmegaOmega

	q <- apply(Omega.jl+PsiPsi.jl,2,sum)/m/(n+1)


	if(method == Method_H0)
	{
		new.mu <- Maximize.H0(data, PsiPsiOmegaOmega, PsiOmega.kl, K, L, n, m, T)
	}
	if(method == Method_H1Equal)
	{
		new.mu <- Newton.Equal(data, sigma, PsiPsiOmegaOmega, K, L, n, m, T)
	}
	if(method == Method_H1TissueInteraction)
	{
		new.mu <- Newton.TissueInteraction(data, mu, sigma, PsiPsiOmegaOmega, K, L, n, m, T)
	}
	if(method == Method_H1ClusterInteraction)
	{
		new.mu <- Newton.ClusterInteraction(data, mu, sigma, PsiPsiOmegaOmega, K, L, n, m, T)
	}

	PsiOmegaData <- MyOuterSigma(data, mu, PsiPsiOmegaOmega)
	PsiOmegaData <- array(PsiOmegaData,dim=c(n,m,K,L,T,T))
	PsiOmegaData.klts <- apply(PsiOmegaData, 3:6, sum)
	new.sigma <- apply(PsiOmegaData.klts, 3:4, FUN=function(x){x/PsiOmega.kl})
	new.sigma <- array(new.sigma, dim = c(K, L, T, T))

	return(list(q=q, mu=new.mu, sigma=new.sigma))
}

eQTL.InitTheta <- function(data, p, K, L)
{
	n <- dim(data)[1]
	m <- dim(data)[2]
	T <- dim(data)[3]

	tmp.data1 <- array(data, dim = c(n, m * T))
	res1 <- kmeans(tmp.data1, K)

	tmp.data2 <- aperm(data, c(2, 1, 3))
	tmp.data3 <- array(tmp.data2, dim = c(m, n * T))
	res2 <- kmeans(tmp.data3,L)
	
	q <- res2$size / sum(res2$size)	
	
	mu <- array(0, dim = c(K, L, T))
	sigma <- array(0, dim = c(K, L, T, T))
	for(k in 1:K)
	{
		for(l in 1:L)
		{
			tmp.ind1 <- which(res1$cluster == k)
			tmp.ind2 <- which(res2$cluster == l)
			tmp1 <- data[tmp.ind1, tmp.ind2, ]
			if(length(tmp.ind1) != 1 & length(tmp.ind2) != 1)
			{
				mu[k, l, ] <- apply(tmp1, 3, mean)
				tmp2 <- dim(tmp1)
				tmp3 <- array(tmp1, dim = c(tmp2[1] * tmp2[2], T))
				sigma[k, l, , ] <- cov(tmp3)
			}
			else
			{
				mu[k, l, ] <- apply(tmp1, 2, mean)
				sigma[k, l, , ] <- cov(tmp1)
			}
		}
	}

	return(list(mu = mu, sigma = sigma, q = q, row.class = res1$cluster, col.class = res2$cluster) )
}

eQTL.InitTheta.Marker <- function(data, p, K, L)
{
	n <- dim(data)[1]
	m <- dim(data)[2]
	T <- dim(data)[3]

	res1 <- p + 1

	tmp.data2 <- aperm(data, c(2, 1, 3))
	tmp.data3 <- array(tmp.data2, dim = c(m, n * T))
	res2 <- kmeans(tmp.data3,L)
	
	q <- res2$size / sum(res2$size)	
	
	mu <- array(0, dim = c(K, L, T))
	sigma <- array(0, dim = c(K, L, T, T))
	for(k in 1:K)
	{
		for(l in 1:L)
		{
			tmp.ind1 <- which(res1 == k)
			tmp.ind2 <- which(res2$cluster == l)
			tmp1 <- data[tmp.ind1, tmp.ind2, ]
			if(length(tmp.ind1) != 1 & length(tmp.ind2) != 1)
			{
				mu[k, l, ] <- apply(tmp1, 3, mean)
				tmp2 <- dim(tmp1)
				tmp3 <- array(tmp1, dim = c(tmp2[1] * tmp2[2], T))
				sigma[k, l, , ] <- cov(tmp3)
			}
			else
			{
				mu[k, l, ] <- apply(tmp1, 2, mean)
				sigma[k, l, , ] <- cov(tmp1)
			}
		}
	}

	return(list(mu = mu, sigma = sigma, q = q, col.class = res2$cluster) )
}

eQTL.RunEM <- function(data, p, K, L, method=Method_H0)
{
	n <- dim(data)[1]
	m <- dim(data)[2]
	T <- dim(data)[3]

	#if(method==Method_H0)
	#{
	#	#theta <- eQTL.InitTheta(data, p, K, L)	
	#	theta <- result0
	#}
	#if(method!=Method_H0)
	#{
	#	theta <- result0	
	#}
	theta <- result0
	mu <- theta$mu
	sigma <- theta$sigma
	q <- theta$q

	rpt <- 1
	LogLikelihood <- -Inf
	while(TRUE){
		OldLogLikelihood <- LogLikelihood

		EResult <- eQTL.StepE(data, mu, sigma, p, q, K, L, n, m)

		Psi.ik <- EResult$Psi.ik
		Psi.ijkl <- EResult$Psi.ijkl
		Omega.jl <- EResult$Omega.jl
		Omega.ijkl <- EResult$Omega.ijkl
		LogLikelihood <- EResult$LogLikelihood
	
		cat("K:",K,"  L:",L,"	rpt:", rpt, "\n")
		#cat("mu:")
		#print(mu)
		#cat("sigma:")
		#print(sigma)
		#cat("q:", q, "\n")
		cat("LogLikelihood:",LogLikelihood,"\n")
		#cat("\n")
		
		if (is.infinite(LogLikelihood))
			break
		if( (abs(1 - OldLogLikelihood/LogLikelihood) < RELATIVE_DIFF) ){
			cat("quit due to likelihood\n")
			cat("RELATIVE_DIFF:",RELATIVE_DIFF,"\n")
			break
		}
		if( rpt >= REPEAT_LIMIT ){
			cat("quit due to rpt\n")
			break
		}
		if( rpt > 2 & OldLogLikelihood > LogLikelihood ){
			cat("quit due to likelihood decreasing\n")
			break
		}

		theta <- eQTL.StepM(data, mu, sigma, Psi.ik, Psi.ijkl, Omega.jl, Omega.ijkl, K, L, n, m, T, method)

		mu <- theta$mu
		sigma <- theta$sigma
		q <- theta$q

		#identifiable
		#tmp.ind1 <- sort(colSums(mu),index.return=T)$ix
		#mu <- mu[,tmp.ind1]
		#sigma <- sigma[,tmp.ind1]
		#q <- q[tmp.ind1]
		#tmp.ind2 <- sort(rowSums(mu),index.return=T)$ix
		#mu <- mu[tmp.ind2,]
		#sigma <- sigma[tmp.ind2,]
		#p <- p[tmp.ind2]

		rpt <- rpt + 1
	}

	row.class <- max.col(Psi.ik)
	col.class <- max.col(Omega.jl)

	return(list(q=q, row.class = row.class, col.class = col.class, 
				LogLikelihood = LogLikelihood, mu = mu, sigma = sigma) )
}


simuGeno <- function(n.marker, length.chrom, loci, n)
{
	library(qtl)
	mapA <- sim.map(length.chrom, n.marker, include.x=F, eq.spacing=T)
	int <- mapA[[1]]
	set.seed(1111)
	simA <- sim.cross(mapA, n.ind=n, type="bc", map.function="morgan")
	set.seed(val)
	simB <- sim.geno(simA, n.draws=1, step=1, off.end=0, error.prob=0.0001, map.function="morgan")

	qtl <- NULL
	for(i in loci)
	{
		qtl <- cbind(qtl, simB$geno$"1"$draws[,,1][,i+1])
	}
	qtl <- qtl-1

	geno <- simB$geno$'1'$data
	geno <- geno-1
	
	qtl <- as.matrix(qtl)
	return(list(QTL=qtl, Geno=geno, Int=int))
}

simuPheno <- function(mu, sigma, n, m, qtl)
{
	K <- dim(mu)[1]
	L <- dim(mu)[2]
	T <- dim(mu)[3]
	
	pheno <- array(0,dim=c(n, m, T))
	for(l in 1:L)
	{
		for(k in 1:K)
		{
			k.tmp <- k-1
			ind <- which(qtl[,1]==k.tmp)
			#ind <- which(qtl[,l]==k.tmp)
			tmp <- rmvnorm(length(ind) * m, mean = mu[k, l, ], sigma = sigma[k, l, , ])
			pheno[ind, ((l - 1) * m + 1):(l * m), 1:T] <- array(tmp, dim = c(length(ind), m, T))
		}
	}	

	return(list(Pheno=pheno))
}

simuPheno <- function(mu, sigma, n, m, qtl, q)
{
	K <- dim(mu)[1]
	L <- dim(mu)[2]
	T <- dim(mu)[3]

	rep.col <- rmultinom(1, size=m, prob=q)
	
	pheno <- array(0, dim = c(n, m, T))
	for(k in 1:K)
	{
		ind <- which(qtl[, 1] == k - 1)
		Y <- NULL
		for (l in 1:L)
		{
			tmp1 <- rmvnorm(length(ind) * rep.col[l], mean = mu[k, l, ], sigma = sigma[k, l, , ])
			tmp2 <- array(tmp1, dim = c(length(ind), rep.col[l], T))
			Y <- abind(Y, tmp2, along = 2)
		}
		pheno[ind, , ] <- Y
	}

	return(list(Pheno=pheno))
}

simu.eQTL <- function(mu, sigma, q, length.chrom, n.marker, loci, n, m, method)
{
	K <- dim(mu)[1]
	L <- dim(mu)[2]
	T <- dim(mu)[3]

	res1 <- simuGeno(n.marker=n.marker, length.chrom=length.chrom, loci=loci, n)
	qtl <- res1$QTL
	geno <- res1$Geno
	Int <- res1$Int
	res2 <- simuPheno(mu=mu, sigma=sigma, n, m, qtl=qtl, q)
	pheno <- res2$Pheno

	LR.int <- NULL
	mu.int <- NULL
	sigma.int <- NULL
	q.int <- NULL
	for(i in 0:length.chrom)
	{
		print(i)
		tmp <- c(i,Int)
		tmp1 <- sort(tmp)
		ind <- which(tmp1==i)
		if(length(ind)==2)
		{
			p <- cbind(geno[,ind[1]],1-geno[,ind[1]])
		}
		if(length(ind)!=2)
		{
			dist <- i - tmp1[ind-1]
			r1 <- 0.5*(1-exp(-.02*dist))
			r <- 0.5*(1-exp(-.02*(tmp1[ind+1]-tmp1[ind-1])))
			theta <- r1/r
			tmp.p<-geno[,c(ind-1,ind)]
			p1 <- tmp.p[,1]-theta*(tmp.p[,1]-tmp.p[,2])
			p <- cbind(p1,1-p1)
		}
	
		
		if(method==Method_H1Equal)
		{
			if(i==0)
			{
				result0 <<- eQTL.InitTheta.Marker(pheno, p[,1], K, L)
			}
			if(i!=0)
			{
				result0 <<- eQTL.old1
			}
			eQTL.result1 <- eQTL.RunEM(data=pheno, p, K, L, method=Method_H1Equal)
			eQTL.old1 <- eQTL.result1

			if(i==0)
			{
				result0 <<- eQTL.InitTheta.Marker(pheno, p[,1], K, L)
			}
			if(i!=0)
			{
				result0 <<- eQTL.old2
			}
			eQTL.result2 <- eQTL.RunEM(data=pheno, p, K, L, method=Method_H0)
			eQTL.old2 <- eQTL.result2

			LR.int <- c(LR.int, 2*(eQTL.result2$LogLikelihood - eQTL.result1$LogLikelihood))

			mu.int <- cbind(mu.int, as.vector(eQTL.result2$mu))
			sigma.int <- cbind(sigma.int, as.vector(eQTL.result2$sigma))
			q.int <- cbind(q.int, as.vector(eQTL.result2$q))
		}
	}	


	est.QTL <- which.max(LR.int)
	est.mu <- mu.int[,est.QTL]
	est.sigma <- sigma.int[,est.QTL]
	est.q <- q.int[,est.QTL]

	return(list(Mu=est.mu, Sigma=est.sigma, Q=est.q, QTL=est.QTL-1, LR=LR.int))

}