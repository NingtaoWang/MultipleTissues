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
Rcpp::NumericMatrix MyOuterLogPDF(	Rcpp::IntegerVector ind_a, 
					Rcpp::IntegerVector ind_b, 
					Rcpp::NumericVector DataR1, 
					Rcpp::NumericVector DataR2, 
					Rcpp::NumericVector MuR, 
					Rcpp::NumericVector SigmaR){
	Rcpp::IntegerVector DimsDataR1 = DataR1.attr("dim");
	arma::cube DataC1(DataR1.begin(), DimsDataR1[0], DimsDataR1[1], DimsDataR1[2]);
	Rcpp::IntegerVector DimsDataR2 = DataR2.attr("dim");
	arma::cube DataC2(DataR2.begin(), DimsDataR2[0], DimsDataR2[1], DimsDataR2[2]);
	int TC = DimsDataR1[2] - 1;
	int n_col_data1 = DimsDataR1[0];
	int n_col_data2 = DimsDataR2[0];

	Rcpp::IntegerVector DimsMuR = MuR.attr("dim");
	arma::cube MuC(MuR.begin(), DimsMuR[0], DimsMuR[1], DimsMuR[2]);

	Rcpp::IntegerVector DimsSigmaR = SigmaR.attr("dim");
	arma::cube SigmaC(SigmaR.begin(), DimsSigmaR[0], DimsSigmaR[1], DimsSigmaR[2]);
	int T2C = DimsSigmaR[2] - 1;

	int n_ind_a = ind_a.size(), n_ind_b = ind_b.size();
	Rcpp::NumericMatrix xab(n_ind_a, n_ind_b);
	for (int i = 0; i < n_ind_a; i++){
		for (int j = 0; j < n_ind_b; j++){
			double tmp1_sum = 0;
			for (int h_1 = 0; h_1 < n_col_data1; h_1++){
				arma::vec DataVec = DataC1.subcube(h_1, i, 0, h_1, i, TC);
				arma::vec MuVec = MuC.subcube(0, j, 0, 0, j, TC);
				arma::vec SigmaVec = SigmaC.subcube(0, j, 0, 0, j, T2C);
				arma::mat SigmaMat(SigmaVec.begin(), DimsMuR[2], DimsMuR[2]);
				arma::vec DataMuVec = DataVec - MuVec;
				double tmp1 = - 0.5 * DimsDataR1[2] * log(2 * PI);
				double tmp2 = - 0.5 * log(arma::det(SigmaMat));
				tmp1_sum += arma::as_scalar(- arma::trans(DataMuVec) * arma::inv(SigmaMat) * DataMuVec / 2) + tmp1 + tmp2;
			}

			double tmp2_sum = 0;
			for (int h_2 = 0; h_2 < n_col_data2; h_2++){
				arma::vec DataVec = DataC2.subcube(h_2, i, 0, h_2, i, TC);
				arma::vec MuVec = MuC.subcube(1, j, 0, 1, j, TC);
				arma::vec SigmaVec = SigmaC.subcube(1, j, 0, 1, j, T2C);
				arma::mat SigmaMat(SigmaVec.begin(), DimsMuR[2], DimsMuR[2]);
				arma::vec DataMuVec = DataVec - MuVec;
				double tmp1 = - 0.5 * DimsDataR1[2] * log(2 * PI);
				double tmp2 = - 0.5 * log(arma::det(SigmaMat));
				tmp2_sum += arma::as_scalar(- arma::trans(DataMuVec) * arma::inv(SigmaMat) * DataMuVec / 2) + tmp1 + tmp2;
			}
			xab(i,j) = tmp1_sum + tmp2_sum;
		}
	}
	return xab;
}', depends=c("RcppArmadillo"))

Newton.Equal <- function(data, sigma, Omega, p, K, L, n, m, T)
{
	new.mu <- array(NA, dim = c(K, L, T))

	sigma.inverse <- array(NA, dim = c(K, L, T, T))
	for(k in 1:K)
	{
		for(l in 1:L)
		{
			sigma.inverse[k, l, , ] <- solve(sigma[k, l, , ])
		}
	}

	colsum.Omega <- colSums(Omega)

	OmegaSigmaInverseData <- array(NA, dim = c(K, L, T))
	OmegaSigmaInverse <- array(NA, dim = c(K, L, T, T))
	for(k in 1:K)
	{
		ind <- which(p == k - 1)
		tmp.data <- apply(data[ind, , ], 2:3, sum)
		tmp <- t(Omega) %*% tmp.data

		for(l in 1:L)
		{
			OmegaSigmaInverseData[k, l, ] <- sigma.inverse[k, l, , ] %*% tmp[l, ]
			OmegaSigmaInverse[k, l, , ] <- colsum.Omega[l] * length(ind) * sigma.inverse[k, l, , ]
		}
	}

	OmegaSigmaInverseData.Ksums <- apply(OmegaSigmaInverseData, 2:3, sum)
	OmegaSigmaInverse.Ksums <- apply(OmegaSigmaInverse, 2:4, sum)

	for(l in 1:L)
	{
		new.mu[1, l, ] <- solve(OmegaSigmaInverse.Ksums[l, , ]) %*% OmegaSigmaInverseData.Ksums[l, ]
		new.mu[2, l, ] <- new.mu[1, l, ]
	}
	
	return(new.mu)
}

Newton.TissueInteraction <- function(data, mu, sigma, Omega, p, K, L, n, m, T)
{
	cl <- apply(mu[1, , ] - mu[2, , ], 1, mean)
	mu[2, , ] <- mu[1, , ] - cl %o% rep(1,T)

	new.mu <- array(NA, dim = c(K, L, T))

	sigma.inverse <- array(NA, dim = c(K, L, T, T))
	for(k in 1:K)
	{
		for(l in 1:L)
		{
			sigma.inverse[k, l, , ] <- solve(sigma[k, l, , ])
		}
	}

	colsum.Omega <- colSums(Omega)

	dQ.dmu.1st <- array(NA, dim = c(K, L, T))
	dQ.dmu.2nd <- array(NA, dim = c(K, L, T, T))
	for(k in 1:K)
	{
		ind <- which(p == k - 1)
		tmp.data <- apply(data[ind, , ], 2:3, sum)
		tmp <- t(Omega) %*% tmp.data

		for(l in 1:L)
		{
			dQ.dmu.1st[k, l, ] <- sigma.inverse[k, l, , ] %*% tmp[l, ] - colsum.Omega[l] * length(ind) * sigma.inverse[k, l, , ] %*% mu[k, l, ]
			dQ.dmu.2nd[k, l, , ] <- - colsum.Omega[l] * length(ind) * sigma.inverse[k, l, , ]
		}
	}
	dQ.dmu.1st.Ksums <- apply(dQ.dmu.1st, 2:3, sum)
	dQ.dmu.2nd.Ksums <- apply(dQ.dmu.2nd, 2:4, sum)
	dQ.dmu.2nd.Ksums.diag <- apply(dQ.dmu.2nd.Ksums, 1, FUN=function(x){diag(x)})
	dQ.dmu.2nd.Ksums.diag <- t(dQ.dmu.2nd.Ksums.diag)

	dQ.dcl.1st <- - apply(dQ.dmu.1st[2, , ], 1, sum)
	dQ.dcl.2nd <- apply(dQ.dmu.2nd[2, , , ], 1, sum)
	new.cl <- cl - dQ.dcl.1st / dQ.dcl.2nd

	new.mu[1, , ] <- mu[1, , ] - dQ.dmu.1st.Ksums / dQ.dmu.2nd.Ksums.diag
	new.mu[2, , ] <- new.mu[1, , ] - new.cl %o% rep(1,T)
	
	return(new.mu)
}


Newton.ClusterInteraction <- function(data, mu, sigma, Omega, p, K, L, n, m, T)
{
	c.vector <- apply(mu[1, , ] - mu[2, , ], 2, mean)
	mu[2, , ] <- mu[1, , ] - rep(1,L) %o% c.vector

	new.mu <- array(NA, dim = c(K, L, T))

	sigma.inverse <- array(NA, dim = c(K, L, T, T))
	for(k in 1:K)
	{
		for(l in 1:L)
		{
			sigma.inverse[k, l, , ] <- solve(sigma[k, l, , ])
		}
	}

	colsum.Omega <- colSums(Omega)

	dQ.dmu.1st <- array(NA, dim = c(K, L, T))
	dQ.dmu.2nd <- array(NA, dim = c(K, L, T))
	for(k in 1:K)
	{
		ind <- which(p == k - 1)
		tmp.data <- apply(data[ind, , ], 2:3, sum)
		tmp <- t(Omega) %*% tmp.data

		for(l in 1:L)
		{
			dQ.dmu.1st[k, l, ] <- sigma.inverse[k, l, , ] %*% tmp[l, ] - colsum.Omega[l] * length(ind) * sigma.inverse[k, l, , ] %*% mu[k, l, ]
			dQ.dmu.2nd[k, l, ] <- - colsum.Omega[l] * length(ind) * diag(sigma.inverse[k, l, , ])
		}
	}
	dQ.dmu.1st.Ksums <- apply(dQ.dmu.1st, 2:3, sum)
	dQ.dmu.2nd.Ksums <- apply(dQ.dmu.2nd, 2:3, sum)

	dQ.dc.1st <- - apply(dQ.dmu.1st[2, , ], 2, sum)
	dQ.dc.2nd <- apply(dQ.dmu.2nd[2, , ], 2, sum)
	new.c.vector <- c.vector - dQ.dc.1st / dQ.dc.2nd

	new.mu[1, , ] <- mu[1, , ] - dQ.dmu.1st.Ksums / dQ.dmu.2nd.Ksums
	new.mu[2, , ] <- new.mu[1, , ] - rep(1,L) %o% new.c.vector
	
	return(new.mu)	
}

eQTL.StepE <- function(data, mu, sigma, p, q, K, L, T, n, m)
{
	ind1 <- which(p==0)
	ind2 <- which(p==1)
	data1 <- data[ind1,,]
	data2 <- data[ind2,,]
	sigma.trans <- array(sigma, dim = c(K, L, T^2))
	Fjl <- MyOuterLogPDF(1:m, 1:L, data1, data2, mu, sigma.trans)
	Fjl[Fjl == Inf] <- max(Fjl[is.finite(Fjl)])
	Fjl[Fjl == -Inf] <- min(Fjl[is.finite(Fjl)])

	F.log.jl.adjtmp <- apply(Fjl,1,
			FUN=function(x){
				tmpx <- max(x)-1400
				x[which(x<tmpx)]<-tmpx
				return(x)})	
	F.log.jl.adj <- apply(F.log.jl.adjtmp,2,FUN=function(x){x-max(x)/2-min(x)/2})
	F.log.jl.adj <- t(F.log.jl.adj)
	

	F.j.log <- rep(0,m)
	Omega <- array(0,dim=c(m,L))
	for(j in 1:m){
		Fjtmp <- exp(F.log.jl.adj[j,])
		OmegaFtmp <- q %*% Fjtmp
		Omega[j,] <-(q * Fjtmp)/OmegaFtmp
		F.j.log[j] <- log(OmegaFtmp)		
	}
	Omega[is.nan(Omega)] <- 1e-300
	Omega[Omega < 1e-300] <- 1e-300

	F.log.j.mean <- apply(F.log.jl.adjtmp, 2, FUN=function(x){max(x)/2+min(x)/2})
	#F.jl.tmp <- apply(F.log.jl.adj,1,FUN=function(x){exp(x)*q})
	#F.j.tmp <- apply(F.jl.tmp,2,sum)
	#F.j.log <- log(F.j.tmp)
	LogLikelihood <- sum(F.j.log)+sum(F.log.j.mean)

	return(list(Omega = Omega, LogLikelihood = LogLikelihood))
}

eQTL.StepM <- function(data, mu, sigma, p, Omega, K, L, T, n, m, method)
{
	new.mu <- array(NA, dim = c(K, L, T))
	new.sigma <- array(NA, dim = c(K, L, T, T))

	colsum.Omega <- colSums(Omega)
	q <- colsum.Omega/m

	if(method == Method_H0)
	{
		for(k in 1:K)
		{
			ind <- which(p == k - 1)
			tmp.data <- apply(data[ind, , ], 2:3, sum)
			tmp <- t(Omega) %*% tmp.data
			new.mu[k, , ] <- tmp / length(ind) / colsum.Omega
		}
	}
	if(method == Method_H1Equal)
	{
		new.mu <- Newton.Equal(data, sigma, Omega, p, K, L, n, m, T)
	}
	if(method == Method_H1TissueInteraction)
	{
		new.mu <- Newton.TissueInteraction(data, mu, sigma, Omega, p, K, L, n, m, T)
	}
	if(method == Method_H1ClusterInteraction)
	{
		new.mu <- Newton.ClusterInteraction(data, mu, sigma, Omega, p, K, L, n, m, T)
	}

	#estimate sigma
	for(k in 1:K)
	{
		ind <- which(p == k - 1)
		Data.Data <- apply(data[ind, , ], 1:2, FUN = function(x){x %o% x})
		Omega.Data.Data.part1 <- apply(Data.Data, c(1, 3), sum) %*% Omega
		Omega.Data.Data.part1 <- array(t(Omega.Data.Data.part1), dim = c(L, T, T))
		Omega.Data.Data.part2.tmp <- t(Omega) %*% apply(data[ind, , ], 2:3, sum)
		Omega.Data.Data.part2 <- array(NA, dim = c(L, T, T))
		for(l in 1:L)
		{
			Omega.Data.Data.part2[l, , ] <- Omega.Data.Data.part2.tmp[l, ] %o% new.mu[k, l, ] + new.mu[k, l, ] %o% Omega.Data.Data.part2.tmp[l, ] 

		}
		Omega.Data.Data.part3 <- length(ind) * colsum.Omega * t(apply(new.mu[k, , ], 1, FUN = function(x){x %o% x}))
		Omega.Data.Data.part3 <- array(Omega.Data.Data.part3, dim = c(L, T, T))
		Omega.Data.Data <- Omega.Data.Data.part1 - Omega.Data.Data.part2 + Omega.Data.Data.part3
		new.sigma[k, , , ] <- Omega.Data.Data / length(ind) / colsum.Omega		
	}

	return(list(q = q, mu = new.mu, sigma = new.sigma))
}

eQTL.InitTheta.Marker <- function(data, p, K, L)
{
	set.seed(1)
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

eQTL.InitTheta.Marker.H1Equal <- function(data, p, K, L)
{
	set.seed(1)
	n <- dim(data)[1]
	m <- dim(data)[2]
	T <- dim(data)[3]

	res1 <- p + 1

	tmp.data2 <- aperm(data, c(2, 1, 3))
	tmp.data3 <- array(tmp.data2, dim = c(m, n * T))
	res2 <- kmeans(tmp.data3, L)
	
	q <- res2$size / sum(res2$size)	
	
	mu <- array(0, dim = c(K, L, T))
	sigma <- array(0, dim = c(K, L, T, T))
	for(k in 1:K)
	{
		for(l in 1:L)
		{
			tmp.ind1 <- which(res1 == k)
			tmp.ind2 <- which(res2$cluster == l)
			tmp1 <- data[, tmp.ind2,]
			tmp0 <- data[tmp.ind1, tmp.ind2,]
			if(length(tmp.ind2) != 1)
			{
				mu[k, l, ] <- apply(tmp1, 3, mean)
				tmp2 <- dim(tmp0)
				tmp3 <- array(tmp0, dim = c(tmp2[1] * tmp2[2], T))
				sigma[k, l, , ] <- cov(tmp3)
			}
			else
			{
				mu[k, l, ] <- apply(tmp1, 2, mean)
				sigma[k, l, , ] <- cov(tmp0)
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

	#theta <- eQTL.InitTheta.Marker(data, p, K, L)	
	theta <- result0
	mu <- theta$mu
	sigma <- theta$sigma
	q <- theta$q

	LogLikelihood <- -Inf

	rpt <- 1
	while(TRUE){
		OldLogLikelihood <- LogLikelihood

		EResult <- eQTL.StepE(data, mu, sigma, p, q, K, L, T, n, m)
		Omega <- EResult$Omega
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
			#cat("quit due to likelihood\n")
			#cat("RELATIVE_DIFF:",RELATIVE_DIFF,"\n")
			break
		}
		if( rpt >= REPEAT_LIMIT ){
			cat("quit due to rpt\n")
			break
		}
		if( rpt > 2 & OldLogLikelihood > LogLikelihood ){
			LogLikelihood <- OldLogLikelihood
			theta <- Oldtheta
			mu <- theta$mu
			sigma <- theta$sigma
			q <- theta$q
			cat("quit due to likelihood decreasing\n")
			break
		}
		Oldtheta <- theta		

		theta <- eQTL.StepM(data, mu, sigma, p, Omega, K, L, T, n, m, method)
		mu <- theta$mu
		sigma <- theta$sigma
		q <- theta$q

		#identifiable
		#tmp.ind <- sort(mu[1,],index.return=T)$ix
		#mu <- mu[,tmp.ind]
		#sigma <- sigma[,tmp.ind]
		#q <- q[tmp.ind]


		rpt <- rpt + 1
	}

	col.class <- max.col(Omega)

	return(list(q=q, col.class = col.class, 
				LogLikelihood = LogLikelihood, mu = mu, sigma=sigma) )
}


simuGeno <- function(n.marker, length.chrom, loci, n)
{
	library(qtl)
	mapA <- sim.map(length.chrom, n.marker, include.x=F, eq.spacing=T)
	#set.seed(1111)
	simA <- sim.cross(mapA, n.ind=n, type="bc", map.function="morgan")
	#set.seed(val)
	simB <- sim.geno(simA, n.draws=1, step=1, off.end=0, error.prob=0.0001, map.function="morgan")


	qtl <- NULL
	geno <- NULL
	int <- NULL
	for(j in 1:length(loci))
	{
		i <- loci[j]
		qtl <- cbind(qtl, simB$geno[[j]]$draws[,,1][,i+1])
		geno <- cbind(geno, simB$geno[[j]]$data)
		int <- c(int, mapA[[j]])
	}
	qtl <- qtl-1
	geno <- geno-1
	
	qtl <- as.matrix(qtl)
	return(list(QTL=qtl, Geno=geno, Int=int))
}

simuPheno <- function(mu, sigma, n, m, qtl)
{
	K <- dim(mu)[1]
	L <- dim(mu)[2]
	T <- dim(mu)[3]
	
	pheno <- array(0,dim=c(n, L*m, T))
	for(l in 1:L)
	{
		for(k in 1:K)
		{
			k.tmp <- k-1
			#ind <- which(qtl[,1]==k.tmp)
			ind <- which(qtl[,l]==k.tmp)
			tmp <- rmvnorm(length(ind) * m, mean = mu[k, l, ], sigma = sigma[k, l, , ])
			pheno[ind, ((l - 1) * m + 1):(l * m), 1:T] <- array(tmp, dim = c(length(ind), m, T))
		}
	}	

	return(list(Pheno=pheno))
}

simu.cutoff <- function(pheno, geno, K, L, rep, method)
{
	LDperm <- NULL
	n <- dim(geno)[1]
	p <- c(rep(0, round(n / 2)), rep(1, n - round(n / 2)))
	for(i in 1:rep)
	{
		print(i)
		set.seed(i)
		ind <- sample(1:dim(pheno)[1])
		tmp.pheno <- pheno[ind,,]
		if(method == Method_H1Equal)
		{
			result0 <<- eQTL.InitTheta.Marker(data=tmp.pheno, p, K, L=3)
			eQTL.result1 <-  eQTL.RunEM(data=tmp.pheno, p, K, L=3, method=Method_H1TissueInteraction)
			result0 <<- eQTL.result1
			eQTL.result0 <-  eQTL.RunEM(data=tmp.pheno, p, K, L=3, method=Method_H1Equal)

			LDperm[i] <- eQTL.result1$LogLikelihood - eQTL.result0$LogLikelihood
		}

		if(method == Method_H1TissueInteraction)
		{
			result0 <<- eQTL.InitTheta.Marker(data=tmp.pheno, p, K, L=3)
			eQTL.result1 <-  eQTL.RunEM(data=tmp.pheno, p, K, L=3, method=Method_H0)
			result0 <<- eQTL.result1
			eQTL.result0 <-  eQTL.RunEM(data=tmp.pheno, p, K, L=3, method=Method_H1TissueInteraction)
			LDperm[i] <- eQTL.result1$LogLikelihood - eQTL.result0$LogLikelihood
		}
	}

	return(LDperm)
}

