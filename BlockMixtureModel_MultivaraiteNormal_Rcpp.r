REPEAT_LIMIT <<- 1001
RELATIVE_DIFF <<- 0.000001

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
				Rcpp::NumericVector PsiPsi,
				Rcpp::NumericVector OmegaOmega) {
	Rcpp::IntegerVector DimsDataR = DataR.attr("dim");
	arma::cube DataC(DataR.begin(), DimsDataR[0], DimsDataR[1], DimsDataR[2]);

	Rcpp::IntegerVector DimsPsiPsi = PsiPsi.attr("dim");
	
	int n_iter = DimsPsiPsi[0] * DimsPsiPsi[1] * DimsPsiPsi[2] * DimsPsiPsi[3], 
	n_iter_ijk = DimsPsiPsi[0] * DimsPsiPsi[1] * DimsPsiPsi[2], 
	n_iter_ij = DimsPsiPsi[0] * DimsPsiPsi[1], 
	n_iter_i = DimsPsiPsi[0],
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
			Output(h,t) = DataC(i,j,t) * (PsiPsi(h) + OmegaOmega(h));
		}
	}
	return Output;      
}', depends=c("RcppArmadillo"))

cppFunction('
arma::mat MyOuterSigma( Rcpp::NumericVector DataR,
				Rcpp::NumericVector MuR,
				Rcpp::NumericVector PsiPsi,
				Rcpp::NumericVector OmegaOmega) {
	Rcpp::IntegerVector DimsDataR = DataR.attr("dim");
	arma::cube DataC(DataR.begin(), DimsDataR[0], DimsDataR[1], DimsDataR[2]);

	Rcpp::IntegerVector DimsMuR = MuR.attr("dim");
	arma::cube MuC(MuR.begin(), DimsMuR[0], DimsMuR[1], DimsMuR[2]);

	Rcpp::IntegerVector DimsPsiPsi = PsiPsi.attr("dim");
	
	int n_iter = DimsPsiPsi[0] * DimsPsiPsi[1] * DimsPsiPsi[2] * DimsPsiPsi[3], 
	n_iter_ijk = DimsPsiPsi[0] * DimsPsiPsi[1] * DimsPsiPsi[2], 
	n_iter_ij = DimsPsiPsi[0] * DimsPsiPsi[1], 
	n_iter_i = DimsPsiPsi[0],
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
			Output(h,ts) = (PsiPsi(h) + OmegaOmega(h)) * (DataC(i, j, t) - MuC(k, l, t)) * (DataC(i, j, s) - MuC(k, l, s));
		}
	}
	return Output;      
}', depends=c("RcppArmadillo"))

cppFunction('
arma::mat MyOuterDMu( 	Rcpp::NumericVector DataR,
				Rcpp::NumericVector MuR,
				Rcpp::NumericVector SigmaR,
				Rcpp::NumericVector Fijkl) {
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
	n_iter_i = DimsDataR[0],
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
		arma::vec SigmaVec = SigmaC.subcube(k, l, 0, k, l, T2C);
		arma::mat SigmaMat(SigmaVec.begin(), DimsMuR[2], DimsMuR[2]);
		arma::vec DataMuVec = DataVec - MuVec;
		arma::vec SigmaInverseDataMuVec = arma::inv(SigmaMat) * DataMuVec;

		for(int t = 0; t < n_iter_t; t++){
			Output(h,t) = Fijkl(h) * SigmaInverseDataMuVec(t);
		}
	}
	return Output;      
}', depends=c("RcppArmadillo"))

cppFunction('
arma::mat MyOuterDSigma(Rcpp::NumericVector DataR,
				Rcpp::NumericVector MuR,
				Rcpp::NumericVector SigmaR,
				Rcpp::NumericVector Fijkl) {
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
	n_iter_i = DimsDataR[0],
	n_iter_ts = DimsDataR[2] * (DimsDataR[2] + 1) / 2;
	arma::mat Output(n_iter, n_iter_ts);
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
		arma::vec SigmaInverseDataMuVec = arma::inv(SigmaMat) * DataMuVec;
		arma::mat DSigma = SigmaInverseDataMuVec * arma::trans(SigmaInverseDataMuVec) - arma::inv(SigmaMat);

		for(int ts = 0; ts < n_iter_ts; ts++){
			double TMPSolv = sqrt(2 * ts + 2.25) - 1.5;
			int t = ceil(TMPSolv);
			int s = ts - t * (t + 1) / 2;
			Output(h, ts) = Fijkl(h) / 2 * DSigma(t, s);
		}
	}
	return Output;      
}', depends=c("RcppArmadillo"))

BM.StepE <- function(data, mu, sigma, p, q, K, L, n, m)
{
	sigma.trans <- array(sigma, dim = c(K, L, T^2))
	Fijkl <- MyOuterMvnorm(data, mu, sigma.trans)
	Fijkl <- array(Fijkl, dim = c(n, m, K, L))
	Fijkl[is.nan(Fijkl)] <- 1e-300
	Fijkl[Fijkl < 1e-300] <- 1e-300

	Psi.ijkl <- apply(Fijkl,1:3,FUN=function(x){x*q/(x%*%q)})
	Psi.ijkl <- aperm(Psi.ijkl, c(2,3,4,1))
	F1.log.ijk <- apply(Fijkl,1:3,FUN=function(x){log(x%*%q)})
	F1.log.ik <- apply(F1.log.ijk,c(1,3),sum)
	F1.log.ik.adjtmp <- apply(F1.log.ik,1,
		FUN=function(x){
			tmpx <- max(x)-1400
			x[which(x<tmpx)]<-tmpx
			return(x)})	
	F1.log.ik.adj <- apply(F1.log.ik.adjtmp,2,FUN=function(x){x-max(x)/2-min(x)/2})
	F1.log.ik.adj <- t(F1.log.ik.adj)
	F1.ik <- apply(F1.log.ik.adj,1,FUN=function(x){exp(x)*p/(exp(x)%*%p)})
	Psi.ik <- t(F1.ik)

	Omega.ijkl <- apply(Fijkl,c(1,2,4),FUN=function(x){x*p/(x%*%p)})
	Omega.ijkl <- aperm(Omega.ijkl, c(2,3,1,4))
	F2.log.ijl <- apply(Fijkl,c(1,2,4),FUN=function(x){log(x%*%p)})
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
	F1.ik.tmp <- apply(F1.log.ik.adj,1,FUN=function(x){exp(x)*p})
	F1.i.tmp <- apply(F1.ik.tmp,2,sum)
	F1.i.log <- log(F1.i.tmp)
	F2.log.j.mean <- apply(F2.log.jl.adjtmp,2,FUN=function(x){max(x)/2+min(x)/2})
	F2.jl.tmp <- apply(F2.log.jl.adj,1,FUN=function(x){exp(x)*q})
	F2.j.tmp <- apply(F2.jl.tmp,2,sum)
	F2.j.log <- log(F2.j.tmp)
	LogLikelihood <- sum(F1.i.log)+sum(F2.j.log)+sum(F1.log.i.mean)+sum(F2.log.j.mean)

	return(list(Psi.ik = Psi.ik, Psi.ijkl = Psi.ijkl, 
			Omega.jl = Omega.jl, Omega.ijkl = Omega.ijkl, LogLikelihood = LogLikelihood))
}

BM.StepM <- function(data, Psi.ik, Psi.ijkl, Omega.jl, Omega.ijkl, K, L, n, m, T)
{
	PsiPsi <- MyOuterPsiPsi(Psi.ik, Psi.ijkl)
	PsiPsi <- array(PsiPsi, dim = c(n, m, K, L))
	PsiPsi.jl <- apply(PsiPsi,c(2,4),sum)

	OmegaOmega <- MyOuterOmegaOmega(Omega.jl, Omega.ijkl)
	OmegaOmega <- array(OmegaOmega, dim = c(n, m, K, L))
	OmegaOmega.ik <- apply(OmegaOmega,c(1,3),sum)

	PsiOmega.kl <- apply(OmegaOmega+PsiPsi, 3:4, sum)

	p <- apply(Psi.ik+OmegaOmega.ik,2,sum)/n/(m+1)
	q <- apply(Omega.jl+PsiPsi.jl,2,sum)/m/(n+1)

	PsiOmegaData <- MyOuterMu(data, PsiPsi, OmegaOmega)
	PsiOmegaData <- array(PsiOmegaData, dim=c(n, m, K, L, T))
	PsiOmegaData.klt <- apply(PsiOmegaData, 3:5, sum)
	mu <- apply(PsiOmegaData.klt, 3, FUN = function(x){x / PsiOmega.kl})
	mu <- array(mu, dim = c(K, L, T))

	PsiOmegaData <- MyOuterSigma(data, mu, PsiPsi, OmegaOmega)
	PsiOmegaData <- array(PsiOmegaData,dim=c(n,m,K,L,T,T))
	PsiOmegaData.klts <- apply(PsiOmegaData, 3:6, sum)
	sigma <- apply(PsiOmegaData.klts, 3:4, FUN=function(x){x/PsiOmega.kl})
	sigma <- array(sigma, dim = c(K, L, T, T))

	return(list(p=p, q=q, mu=mu, sigma=sigma) )
}

BM.InitTheta <- function(data, K, L)
{
	n <- dim(data)[1]
	m <- dim(data)[2]
	T <- dim(data)[3]

	tmp.data1 <- array(data, dim = c(n, m * T))
	res1 <- kmeans(tmp.data1, K)

	tmp.data2 <- aperm(data, c(2, 1, 3))
	tmp.data3 <- array(tmp.data2, dim = c(m, n * T))
	res2 <- kmeans(tmp.data3,L)
	
	p <- res1$size / sum(res1$size)
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

	return(list(mu = mu, sigma = sigma, p = p, q = q, row.class = res1$cluster, col.class = res2$cluster) )
}

BM.RunEM <- function(data, K, L)
{
	n <- dim(data)[1]
	m <- dim(data)[2]
	T <- dim(data)[3]

	#theta <- BM.InitTheta(data, K, L)
	theta <- result0
	mu <- theta$mu
	sigma <- theta$sigma
	p <- theta$p
	q <- theta$q

	rpt <- 1
	LogLikelihood <- -Inf
	while(TRUE){
		OldLogLikelihood <- LogLikelihood

		EResult <- BM.StepE(data, mu, sigma, p, q, K, L, n, m)

		Psi.ik <- EResult$Psi.ik
		Psi.ijkl <- EResult$Psi.ijkl
		Omega.jl <- EResult$Omega.jl
		Omega.ijkl <- EResult$Omega.ijkl
		LogLikelihood <- EResult$LogLikelihood
	
		cat("K:",K,"  L:",L,"	rpt:", rpt, "\n")
		cat("mu:", mu, "\n")
		cat("sigma:", sigma, "\n")
		cat("p:", p, "\n")
		cat("q:", q, "\n")
		cat("LogLikelihood:",LogLikelihood,"\n")
		cat("\n")
		gc(verbose = T)
		
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
		if( OldLogLikelihood > LogLikelihood ){
			cat("quit due to likelihood decreasing\n")
			break
		}

		theta <- BM.StepM(data, Psi.ik, Psi.ijkl, Omega.jl, Omega.ijkl, K, L, n, m, T)

		mu <- theta$mu
		sigma <- theta$sigma
		p <- theta$p
		q <- theta$q

		#identifiable
		#tmp.ind1 <- sort(rowSums(colSums(mu)),index.return=T)$ix
		#mu <- mu[,tmp.ind1,]
		#sigma <- sigma[,tmp.ind1,,]
		#q <- q[tmp.ind1]
		#tmp.ind2 <- sort(rowSums(mu),index.return=T)$ix
		#mu <- mu[tmp.ind2,,]
		#sigma <- sigma[tmp.ind2,,,]
		#p <- p[tmp.ind2]

		rpt <- rpt + 1
	}

	row.class <- max.col(Psi.ik)
	col.class <- max.col(Omega.jl)

	return(list(p=p, q=q, row.class = row.class, col.class = col.class, 
				LogLikelihood = LogLikelihood, mu = mu, sigma = sigma) )
}

BM.BICsimu <- function(data, mu, sigma, p, q, K, L)
{
	n <- dim(data)[1]
	m <- dim(data)[2]
	T <- dim(data)[3]

	sigma.trans <- array(sigma, dim = c(K, L, T^2))
	Fijkl <- MyOuterMvnorm(data, mu, sigma.trans)
	Fijkl <- array(Fijkl, dim = c(n, m, K, L))
	Fijkl[is.nan(Fijkl)] <- 1e-300
	Fijkl[Fijkl < 1e-300] <- 1e-300

	Fijklt.dmu <- MyOuterDMu(data, mu, sigma.trans, Fijkl)
	Fijklt.dmu <- array(Fijklt.dmu, dim = c(n, m, K, L, T))
	Fijklt.dmu[is.nan(Fijklt.dmu)] <- 1e-300
	Fijklt.dmu[abs(Fijklt.dmu) < 1e-300] <- 1e-300

	TS <- T * {T + 1} / 2 
	Fijklt.dsigma <- MyOuterDSigma(data, mu, sigma.trans, Fijkl)
	Fijklt.dsigma <- array(Fijklt.dsigma, dim=c(n, m, K, L, TS))
	Fijklt.dsigma[is.nan(Fijklt.dsigma)] <- 1e-300
	Fijklt.dsigma[abs(Fijklt.dsigma) < 1e-300] <- 1e-300

	F1.ijk <- apply(Fijkl,1:3,FUN=function(x){x%*%q})
	F1.log.ijk <- log(F1.ijk)
	F1.log.ik <- apply(F1.log.ijk,c(1,3),sum)
	F1.log.ik.adjtmp <- apply(F1.log.ik,1,
		FUN=function(x){
			tmpx <- max(x)-1200
			x[which(x<tmpx)]<-tmpx
			return(x)})	
	F1.log.ik.adj <- apply(F1.log.ik.adjtmp,2,FUN=function(x){x-max(x)/2-min(x)/2})
	F1.log.ik.adj <- t(F1.log.ik.adj)
	F1.ik <- apply(F1.log.ik.adj,1,FUN=function(x){exp(x)/(exp(x)%*%p)})
	F1.ik <- t(F1.ik)
 	drow.dk <- apply(as.matrix(F1.ik[,-K]),2,FUN=function(x){x-F1.ik[,K]})


	F1.i <- apply(F1.log.ik.adj,1,FUN=function(x){exp(x)%*%p})

	F1.ijkl <- apply(Fijkl,1:3,FUN=function(x){x/(x%*%q)})
	F1.ijkl <- aperm(F1.ijkl, c(2,3,4,1))
	F1.ikl <- apply(F1.ijkl,c(1,3,4),sum)
	F1.ikl.2 <- apply(F1.ikl,3,FUN=function(x){x*exp(F1.log.ik.adj)})
	F1.ikl.2 <- array(F1.ikl.2,dim=c(n,K,L))
	F1.il <- apply(F1.ikl.2,c(1,3),FUN=function(x){x%*%p})
	F1.il.2 <- apply(F1.il,2,FUN=function(x){x/F1.i})
	drow.dl <- apply(as.matrix(F1.il.2[,-L]),2,FUN=function(x){x-F1.il.2[,L]})

	
	F1.lijkt.dmu <- apply(Fijklt.dmu,c(1,2,3,5),FUN=function(x){x*q})
	F1.ijklt.dmu <- apply(F1.lijkt.dmu,c(1,5),FUN=function(x){x/F1.ijk})
	F1.ijklt.dmu <- array(F1.ijklt.dmu,dim=c(n,m,K,L,T))
	F1.iklt.dmu <- apply(F1.ijklt.dmu,c(1,3,4,5),sum)
	F1.iklt.dmu.2 <- apply(F1.iklt.dmu,c(3,4),FUN=function(x){x*exp(F1.log.ik.adj)})
	F1.iklt.dmu.2 <- array(F1.iklt.dmu.2,dim=c(n,K,L,T))
	F1.kilt.dmu.3 <- apply(F1.iklt.dmu.2,c(1,3,4),FUN=function(x){x*p})
	F1.iklt.dmu.3 <- apply(F1.kilt.dmu.3,c(1,3,4),FUN=function(x){x/F1.i})
	drow.dmu <- apply(F1.iklt.dmu.3,1,FUN=function(x){as.vector(x)})
	drow.dmu <- t(drow.dmu)


	F1.lijkt.dsigma <- apply(Fijklt.dsigma,c(1,2,3,5),FUN=function(x){x*q})
	F1.ijklt.dsigma <- apply(F1.lijkt.dsigma,c(1,5),FUN=function(x){x/F1.ijk})
	F1.ijklt.dsigma <- array(F1.ijklt.dsigma,dim=c(n,m,K,L,TS))
	F1.iklt.dsigma <- apply(F1.ijklt.dsigma,c(1,3,4,5),sum)
	F1.iklt.dsigma.2 <- apply(F1.iklt.dsigma,c(3,4),FUN=function(x){x*exp(F1.log.ik.adj)})
	F1.iklt.dsigma.2 <- array(F1.iklt.dsigma.2,dim=c(n,K,L,TS))
	F1.kilt.dsigma.3 <- apply(F1.iklt.dsigma.2,c(1,3,4),FUN=function(x){x*p})
	F1.iklt.dsigma.3 <- apply(F1.kilt.dsigma.3,c(1,3,4),FUN=function(x){x/F1.i})
	drow.dsigma <- apply(F1.iklt.dsigma.3,1,FUN=function(x){as.vector(x)})
	drow.dsigma <- t(drow.dsigma)





	F2.ijl <- apply(Fijkl,c(1,2,4),FUN=function(x){x%*%p})
	F2.log.ijl <- log(F2.ijl)
	F2.log.jl <- apply(F2.log.ijl,c(2,3),sum)
	F2.log.jl.adjtmp <- apply(F2.log.jl,1,
		FUN=function(x){
			tmpx <- max(x)-1200
			x[which(x<tmpx)]<-tmpx
			return(x)})	
	F2.log.jl.adj <- apply(F2.log.jl.adjtmp,2,FUN=function(x){x-max(x)/2-min(x)/2})
	F2.log.jl.adj <- t(F2.log.jl.adj)
	F2.jl <- apply(F2.log.jl.adj,1,FUN=function(x){exp(x)/(exp(x)%*%q)})
	F2.jl <- t(F2.jl)
	dcol.dl <- apply(as.matrix(F2.jl[,-L]),2,FUN=function(x){x-F2.jl[,L]})


	F2.j <- apply(F2.log.jl.adj,1,FUN=function(x){exp(x)%*%q})

	F2.ijkl <- apply(Fijkl,c(1,2,4),FUN=function(x){x/(x%*%p)})
	F2.ijkl <- aperm(F2.ijkl, c(2,3,1,4))
	F2.jkl <- apply(F2.ijkl,c(2,3,4),sum)
	F2.jkl.2 <- apply(F2.jkl,2,FUN=function(x){x*exp(F2.log.jl.adj)})
	F2.jlk <- array(F2.jkl.2,dim=c(m,L,K))
	F2.jk <- apply(F2.jlk,c(1,3),FUN=function(x){x%*%q})
	F2.jk.2 <- apply(F2.jk,2,FUN=function(x){x/F2.j})
	dcol.dk <- apply(as.matrix(F2.jk.2[,-K]),2,FUN=function(x){x-F2.jk.2[,K]})


	F2.kijlt.dmu <- apply(Fijklt.dmu,c(1,2,4,5),FUN=function(x){x*p})
	F2.ijlkt.dmu <- apply(F2.kijlt.dmu,c(1,5),FUN=function(x){x/F2.ijl})
	F2.ijlkt.dmu <- array(F2.ijlkt.dmu,dim=c(n,m,L,K,T))
	F2.jlkt.dmu <- apply(F2.ijlkt.dmu,c(2,3,4,5),sum)
	F2.jlkt.dmu.2 <- apply(F2.jlkt.dmu,c(3,4),FUN=function(x){x*exp(F2.log.jl.adj)})
	F2.jlkt.dmu.2 <- array(F2.jlkt.dmu.2,dim=c(m,L,K,T))
	F2.ljkt.dmu.3 <- apply(F2.jlkt.dmu.2,c(1,3,4),FUN=function(x){x*q})
	F2.jlkt.dmu.3 <- apply(F2.ljkt.dmu.3,c(1,3,4),FUN=function(x){x/F2.j})
	F2.jklt.dmu.3 <- aperm(F2.jlkt.dmu.3, c(1,3,2,4))
	dcol.dmu <- apply(F2.jklt.dmu.3, 1, FUN = function(x){as.vector(x)})
	dcol.dmu <- t(dcol.dmu)


	F2.kijlt.dsigma <- apply(Fijklt.dsigma,c(1,2,4,5),FUN=function(x){x*p})
	F2.ijlkt.dsigma <- apply(F2.kijlt.dsigma,c(1,5),FUN=function(x){x/F2.ijl})
	F2.ijlkt.dsigma <- array(F2.ijlkt.dsigma,dim=c(n,m,L,K,TS))
	F2.jlkt.dsigma <- apply(F2.ijlkt.dsigma,c(2,3,4,5),sum)
	F2.jlkt.dsigma.2 <- apply(F2.jlkt.dsigma,c(3,4),FUN=function(x){x*exp(F2.log.jl.adj)})
	F2.jlkt.dsigma.2 <- array(F2.jlkt.dsigma.2,dim=c(m,L,K,TS))
	F2.ljkt.dsigma.3 <- apply(F2.jlkt.dsigma.2,c(1,3,4),FUN=function(x){x*q})
	F2.jlkt.dsigma.3 <- apply(F2.ljkt.dsigma.3,c(1,3,4),FUN=function(x){x/F2.j})
	F2.jklt.dsigma.3 <- aperm(F2.jlkt.dsigma.3, c(1,3,2,4))
	dcol.dsigma <- apply(F2.jklt.dsigma.3, 1, FUN = function(x){as.vector(x)})
	dcol.dsigma <- t(dcol.dsigma)


	drow <- cbind(drow.dk,drow.dl,drow.dmu,drow.dsigma)
	dcol <- cbind(dcol.dk,dcol.dl,dcol.dmu,dcol.dsigma)

	return(list(DROW=drow,DCOL=dcol))
}

BM.SimuData <- function(n, m, mu, sigma, p, q)
{
	K <- length(p)
	L <- length(q)
	T <- dim(mu)[3]

	rep.row <- round(n*p)
	rep.col <- round(m*q)

	Z <- NULL
	for (k in 1:K)
	{
		Y <- NULL
		for (l in 1:L)
		{
			if(rep.row[k] * rep.col[l] != 0)
			{
				tmp1 <- rmvnorm(rep.row[k] * rep.col[l], mean = mu[k, l, ], sigma = sigma[k, l, , ])
				tmp2 <- array(tmp1, dim = c(rep.row[k], rep.col[l], T))
				Y <- abind(Y, tmp2, along = 2)
			}
		}
		Z <- abind(Z, Y, along = 1)
	}

	Z <- as.array(Z)
	return(Z)
}

BM.BIC <- function(data, K, L, rep=30)
{
	n <- dim(data)[1]
	m <- dim(data)[2]
	T <- dim(data)[3]
	TS <- T * {T + 1} / 2

	emResult <- BM.RunEM(data, K, L)

	mu <- emResult$mu
	sigma <- emResult$sigma
	p <- emResult$p
	q <- emResult$q

	H <- NULL
	J <- NULL
	for(i in 1:rep)
	{
		tmpdata <- BM.SimuData.alt(n, m, mu, sigma, p, q)
		BICResult <- BM.BICsimu(tmpdata, mu, sigma, p, q, K, L)
		drow <- BICResult$DROW
		dcol <- BICResult$DCOL

		cat("Simulation:",i,"\n")
		#cat("dRow:",drow,"\n")
		#cat("dCol:",dcol,"\n")
		cat("\n")

		tmp1 <- t(drow)%*%drow+t(dcol)%*%dcol
		H <- cbind(H,as.vector(tmp1))
		tmp2 <- (colSums(drow)+colSums(dcol))%*%t(colSums(drow)+colSums(dcol))
		J <- cbind(J,as.vector(tmp2))
	}
	tmp <- K+L-2+K*L*{T+TS}
	H.simu <- apply(H,1,mean)
	J.simu <- apply(J,1,mean)
	H.simu <- array(H.simu,dim=c(tmp,tmp))
	J.simu <- array(J.simu,dim=c(tmp,tmp))

	tmpdiag <- diag(J.simu%*%solve(H.simu))

	BIC <- -2*emResult$LogLikelihood + (log(n)+log(m)+log(T))*sum(tmpdiag)

	cat("trace:",sum(tmpdiag),"\n")
	cat("relative DF:",sum(tmpdiag)/tmp,"\n")
	cat("BIC:",BIC,"\n")
	cat("LogLikelihood:",emResult$LogLikelihood,"\n")
	cat("\n")
	return(list(H = H, J = J, Hsimu=H.simu, Jsimu=J.simu, BIC = BIC , est = emResult))
}


BM.SimuData.alt <- function(n, m, mu, sigma, p, q)
{
	K <- length(p)
	L <- length(q)
	T <- dim(mu)[3]

	rep.row <- rmultinom(1, size=n, prob=p)
	rep.col <- rmultinom(1, size=m, prob=q)


	Z <- NULL
	for (k in 1:K)
	{
		Y <- NULL
		for (l in 1:L)
		{
			if(rep.row[k] * rep.col[l] != 0)
			{
				tmp1 <- rmvnorm(rep.row[k] * rep.col[l], mean = mu[k, l, ], sigma = sigma[k, l, , ])
				tmp2 <- array(tmp1, dim = c(rep.row[k], rep.col[l], T))
				Y <- abind(Y, tmp2, along = 2)
			}
		}
		Z <- abind(Z, Y, along = 1)
	}

	Z <- as.array(Z)
	return(Z)
}

