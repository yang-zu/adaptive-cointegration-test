/* adaptcoint.ox
   Adaptive test for cointegration with nonstationary volatility.
   Boswijk and Zu (2022), Journal of Business and Economic Statistics.

   Bivariate (p = 2) system, null of no cointegration (rank r = 0).
   Provides the Johansen trace test and the adaptive (spot-covariance weighted)
   LR test, each with wild-bootstrap and volatility-bootstrap p-values.

   Extracted from the paper's Monte Carlo code (lib2020.ox / lib2020k1.ox and
   the driver's pvals/pvalsk1) into a single callable library. */

#include <oxstd.oxh>

/* ================= k = 1 routines (no short-run dynamics) ================= */



tracetestk1(const DX, const LXC)
{
	//trace test
	//unrestricted OLS
	decl T = rows(DX);
	decl Z = LXC;
	decl temp = invertsym(Z'Z); 
	decl psi = temp*Z'DX;
	decl E = DX - Z*psi;

	//restricted OLS
	decl Er = DX;
	decl trace = T*log(determinant(Er'Er/T)/determinant(E'E/T));
	return {trace,Er,E};
}

mle00k1(const DX, const sig2invM)
{
// MLE estimation of restricted model when r=0
	decl T = rows(DX);
	decl p = columns(DX);

	decl i,sig2inv;

	decl E=0;
	decl resid = zeros(T,p);
	for (i=0;i<T;i++)
	{
		sig2inv = sig2invM[(i)*p : ((i+1)*p-1)][];
		resid[i][] = DX[i][];
		E = E + resid[i][] * sig2inv * resid[i][]';
	}
	
	return {E,resid};
}


mle1k1(const DX, const LXC, const sig2invM)
{
// MLE estimation of the unrestricted model, k=k version
	decl T = rows(DX);
	decl p = columns(DX);
	decl k = 1;	
	decl F = zeros((k*p+1)*p,(k*p+1)*p);
	decl G = zeros(p,k*p+1);

	decl z,sig2inv,i;
	for (i=0;i<T;i++)
		{
		z = (LXC[i][])';
		sig2inv = sig2invM[(i)*p : ((i+1)*p-1)][];
		F = F + (z * z') ** sig2inv;
		G = G + sig2inv* DX[i][]' * z';
		}
	decl temp = invertsym(F);
	decl vecpp = temp*vec(G);
	decl pp = shape(vecpp,p,k*p+1);
	decl pai = pp[][0: p];

	decl E = 0;
	decl resid = zeros(T,p);
	for (i=0;i<T;i++)
	{
		sig2inv = sig2invM[(i)*p : ((i+1)*p-1)][];
		resid[i][] = DX[i][] - LXC[i][] * pai';
		E = E + resid[i][] * sig2inv * resid[i][]';
	}
	return{pai,E,resid};
}


LRk1(const DX, const LXC, const sigma2)
{
// input:
//		  X0: T+k by p data matrix
//		  sigma2: error covariance matrices stacked up, Tp by p matrix

// output:
//		  LR statistic for r=0 specifically


	decl T = rows(DX);
	decl p = columns(DX);

	// inverse of covariance matrices
	decl i;
	decl sig2invM = zeros(T*p,p);
	for (i=0;i<T;i++)
	{
	sig2invM[i*p : ((i+1)*p-1)][] = invertsym(sigma2[i*p : ((i+1)*p-1)][]);
	}
	
	// mle estimation of the unrestricted model
	decl pai,E1,resid1;
	[pai,E1,resid1] = mle1k1(DX, LXC, sig2invM);
	
	// mle estimation of restricted models, for r=0 
	decl mu,alpha,beta,E0,resid,lr;

	[E0,resid] = mle00k1(DX, sig2invM);
	lr = E0-E1;
	//return the lr stat and the residuals of the restricted model
	return {lr,resid};			
}

ECM0(const psi, const eps, const DX0, const LXC0, const W0)
{
	// for k>=2 case
	decl p = columns(eps);
	decl T = rows(eps);
	decl k = round(columns(W0)/p)+1;
	
	decl i,t;
	
	decl DX = zeros(T,p);
	decl LXC = zeros(T,p)~ones(T,1);
	decl W = zeros(T,p*(k-1));
	decl WW;
	
	LXC[0][0:p-1] = LXC0[0:p-1]; W[0][] = W0;
	for (t=0; t<T; ++t)
	{
		DX[t][] = W[t][] * psi' + eps[t][];  // using model equation for DXb[t]

		if (t<T-1)
		{
		if (k==2)
		{W[t+1][] = DX[t][];}
		else
		{W[t+1][] = DX[t][]~W[t][0:p*(k-2)-1];}	
		
		LXC[t+1][0:p-1] = LXC[t][0:p-1] + DX[t][];
		}

				
	}
	return {DX,LXC,W};
}


pval_wbk1(const LRt, const LRa, const E,const k, const sigma2, const LXC0)
{
	decl T = rows(E);
	decl d = columns(E);


	decl seed, p, m, DX, LXC, lrt,lra,pt,pa;
	seed = ranseed(0);
	ranseed(12345678);

	pt = 0;
	pa = 0;
	m = 499;

	decl eta,eps,X,i;
	decl Ebt,Erbt;
	for (i=0; i<m; ++i)
	{
		eta = rann(T,1);
		eps = eta.*E;
		eps[0][]=LXC0[0:d-1];
		X = cumulate(eps, unit(d,d));
		DX = X - lag0(X,1);
		LXC = lag0(X,1)~ones(T,1);
		DX = DX[k:T-1][];
		LXC = LXC[k:T-1][];
		lrt = tracetestk1(DX,LXC);
		lra = LRk1(DX, LXC, sigma2);
		pt += (lrt[0].>LRt)/m;
		pa += (lra[0].>LRa)/m;
		//println(i);
	}	
	ranseed(seed);
	return double(pt)~double(pa);
}



pval_vbk1(const LRt, const LRa, const sigma,const k, const sigma2, const LXC0)
{

	// input sigma is n by 3 matrix	of the lower triangular version of the square root given by sig2tosig()
	decl T = rows(sigma);
	decl d = 2;

	decl seed, p, m, DX, LXC, W, lrt,lra,pt,pa;
	seed = ranseed(0);
	ranseed(12345678);

	pt = 0;
	pa = 0;
	m = 499;

	decl eta,eps,X,i;
	for (i=0; i<m; ++i)
	{
		eta = rann(T,2);
		eps = sigma[][0].*eta[][0]~(sigma[][1].*eta[][1]+sigma[][2].*eta[][0]);
		eps[0][]=LXC0[0:d-1];
		X = cumulate(eps, unit(d,d));
		DX = X - lag0(X,1);
		LXC = lag0(X,1)~ones(T,1);

		DX = DX[k:T-1][];
		LXC = LXC[k:T-1][];
		lrt = tracetestk1(DX,LXC);
		lra = LRk1(DX, LXC, sigma2);
		pt += (lrt[0].>LRt)/m;
		pa += (lra[0].>LRa)/m;		
		//println(i);
	}	
	ranseed(seed);
	return double(pt)~double(pa);
}

/* ============ general routines: spot cov, order, trace, adaptive LR, bootstrap ============ */
spcov_s(const tgrid, const t, const x2, const h)
{
// tgrid: 1 by m grid of time to be evaluate the covariance matrix
// t: 1 by n grid of observation time
// x: 2 by n of vector observations
// h: bandwidth
// output:3 by m matrix, each column is the estimated spot covariance matrix corresponding to tgrid
// this version only works for 2 dim case
	
	decl tt = t - tgrid'; 		//m by n matrix
	decl ktth = densn(tt/h)/h;	//m by n matrix, each row corrsponds to one time point to be estimated
	decl w = ktth./sumr(ktth);	 //m by n matrix weights
	//decl m = columns(tgrid);  // # points to be evaluated
	decl spv = x2*w'; // 3 by n times n by m, gives 3 by m matrix, each column is a estimated spot cov matrix correspondint to tgrid
	return spv;

}


cvbw_s(const t, const x2)
{
// t: 1 by n grid of observation time
// x: 3 by n of vector observations
// use the leave-one-out version of the spot covariance estimator
// this is based on the previous simulation, not for matrix to speed up
	
	decl n = columns(x2);
	decl dt = t[1]-t[0];
	decl hgrid = range(dt,dt*n/6,(dt*n/6-dt)/100);
	decl cv = zeros(columns(hgrid),1);
	decl i,tt,ktth,w,spv,D,Dnorm;
	for (i=0; i<columns(hgrid); i++)
	{
		tt = t - t'; 		
		ktth = densn(tt/hgrid[i])/hgrid[i];	//m by n matrix, each row corrsponds to one time point to be estimated
		w = ktth./sumr(ktth);	 //n by n matrix weights
		spv = x2*w';
		D = meanr(((x2-spv)./(1-diagonal(w))).^2);	  //3 by n difference matrix
		cv[i] = D[0]+D[1]+2*D[2]; 					  
	}
	decl h = hgrid[mincindex(cv)];
	return h;
}



sig2tosig(const pShat)
{
// taking square root of a 2 by 2 covariance matrix
// Shat is a 3 by n matrix with each column a stacked covariance matrix
decl Shat = pShat';
decl sig = sqrt(Shat[][0]~(Shat[][1]-(Shat[][2].^2)./Shat[][0]))~(Shat[][2]./sqrt(Shat[][0])); //sig is lower triangular, with the last row as the off diagonal element
return sig';
}

DGP(const alpha, const beta, const psi, const psigma)
{
	// sigma as a 3 by n matrix, with output to be x by n matrix
    decl n, eta, eps, An, X, DX, Z;
	decl sigma = psigma';
	n = rows(sigma);
    eta = rann(n,2);
	eps = sigma[][0].*eta[][0]~(sigma[][1].*eta[][1]+sigma[][2].*eta[][0]);
	//square root matrix organize in this way, making the variance matrix in right form
	An = unit(2) + alpha*beta'+psi;
	X = cumulate(eps, {An,-psi});
	return X;
}

icorder(const X, const kmax, const whichic)
{
	// order selection by BIC
	// minimum k is 1, as we work with ECM form, need k>=1
	decl n, p, T, Xtr, k,ct, consZ, preE, Sig, SC, sign,dimp;
	
	n = rows(X);
	p = columns(X);
	T = n-kmax;

	if (whichic == 1)//AIC
	{ct = 2;}
	else if (whichic == 2)//BIC
	{ct = log(T);}
	else if (whichic == 3)
	{ct = 2*log(log(T));}//HQIC

	
	Xtr = X[kmax:n-1][];
	consZ = ones(T,1);
	SC = zeros(kmax,1);
	decl temp;
	for (k=1; k<=kmax; ++k)
	{
		consZ = consZ~X[kmax-k:n-1-k][];//getting larger for each k. ensure n-kmax obs in every regression
		temp = 	invertsym(consZ'consZ);
		preE = Xtr - consZ*temp*consZ'Xtr;
		Sig = preE'preE/T;
		SC[k-1] = logdet(Sig,&sign) + p*columns(consZ)*ct/T;				
	}
	return mincindex(SC)+1;

}


tracetest(const DX, const LXC, const W)
{
	//trace test for r=0
	//unrestricted OLS
	decl T = rows(DX);
	decl Z = LXC~W;
	decl temp = invertsym(Z'Z);
	decl psi = temp*Z'DX;
	decl E = DX - Z*psi;

	//restricted OLS
	decl Z1 = W;	   // for restricted constant case, no ones(T,1) needed
	decl temp1 = invertsym(Z1'Z1);
	decl psir = temp1*Z1'DX;
	decl Er = DX - Z1*psir;
	decl trace = T*log(determinant(Er'Er/T)/determinant(E'E/T));
	return {trace,psir',Er,E};//return restricted coefficient, residuals for bootstrap, unrestricted residuals for vol estimation.
}







mle00(const DX, const W, const sig2invM)
{
// MLE estimation of restricted model when r=0
// mle0 and mle0 are actually fast to run, faster than the matrix version of wls implementation
	decl T = rows(DX);
	decl p = columns(DX);
	decl k = round(columns(W)/p)+1;


	decl A0 = zeros((p*(k-1))*p,(p*(k-1))*p);
	decl B0 = zeros(p,p*(k-1));
	decl Z = W;

	decl i,sig2inv;

	for (i=0;i<T;i++)
		{
		sig2inv = sig2invM[(i)*p : ((i+1)*p-1)][];
		A0 = A0 + (Z[i][]' * Z[i][]) ** sig2inv;
		B0 = B0 + sig2inv* DX[i][]' * Z[i][];
		}
	decl temp = invertsym(A0);
	decl vecpsi = temp * vec(B0);
	decl psimu = shape(vecpsi,p,p*(k-1));
	decl psi = psimu;
	

	decl E=0;
	decl resid = zeros(T,p);
	for (i=0;i<T;i++)
	{
		sig2inv = sig2invM[(i)*p : ((i+1)*p-1)][];
		resid[i][] = DX[i][]  - ( psimu * Z[i][]' )';
		E = E + resid[i][] * sig2inv * resid[i][]';
	}
//	decl mu = psimu[][p*(k-1)];	
	return {psi,E,resid}; //return coef estimate and residuals
}

mle1(const DX, const LXC, const W, const sig2invM)
{
// MLE estimation of the unrestricted model, k=k version
	decl T = rows(DX);
	decl p = columns(DX);
	decl k = round(columns(W)/p)+1;	
	decl F = zeros((k*p+1)*p,(k*p+1)*p);
	decl G = zeros(p,k*p+1);

	decl z,sig2inv,i;
	for (i=0;i<T;i++)
		{
		z = (LXC[i][] ~ W[i][])';
		sig2inv = sig2invM[(i)*p : ((i+1)*p-1)][];
		F = F + (z * z') ** sig2inv;
		G = G + sig2inv* DX[i][]' * z';
		}
	decl temp = invertsym(F);
	decl vecpp = temp*vec(G);
	decl pp = shape(vecpp,p,k*p+1);
	decl pai = pp[][0: p];
	decl psi = pp[][p+1:];

	decl E = 0;
	decl resid = zeros(T,p);
	for (i=0;i<T;i++)
	{
		sig2inv = sig2invM[(i)*p : ((i+1)*p-1)][];
		resid[i][] = DX[i][] - LXC[i][] * pai' - W[i][] * psi';
		E = E + resid[i][] * sig2inv * resid[i][]';
	}
	return{pai,psi,E,resid};//return coef estimate and residuals
}

LR(const DX, const LXC, const W, const sigma2)
{
// input:
//		  X0: T+k by p data matrix
//		  sigma2: error covariance matrices stacked up, Tp by p matrix

// output:
//		  LR statistic for a specific r


	decl T = rows(DX);
	decl p = columns(DX);


	// inverse of covariance matrices
	decl i;
	decl sig2invM = zeros(T*p,p);
	for (i=0;i<T;i++)
	{
	sig2invM[i*p : ((i+1)*p-1)][] = invertsym(sigma2[i*p : ((i+1)*p-1)][]);
	}

	
	// mle estimation of the unrestricted model
	decl pai,psi,E1,resid1;
	[pai,psi,E1,resid1] = mle1(DX, LXC, W, sig2invM);
	
	// mle estimation of restricted models, for r=0 case, specific to simulation 2019
	decl psir,E0,resid,lr;

	[psir,E0,resid] = mle00(DX, W, sig2invM);
	lr = E0-E1;

	//return the coefficient of restricted model for bootstrap, no residual needed for vol estimation
	return {lr,psir,resid};
			
}




pval_wb(const LRt, const LRa, const E,const psi, const sigma2, const DX0, const LXC0, const W0)
{
	// this implementation of wild boostrap trace test p value simulate the lag k dynamics
	decl T = rows(E);
//	decl d = columns(E);
	
	decl seed, p, m, DX, LXC, W, lrt,lra,pt,pa;
	seed = ranseed(0);
	ranseed(12345678);

	pt = 0;
	pa = 0;
	m = 499;

	decl eta,eps,X;
	decl psirbt,Ebt,Erbt,i;
	for (i=0; i<m; ++i)
	{
		eta = rann(T,1);
		eps = eta.*E;
		[DX,LXC,W] = ECM0(psi, eps, DX0, LXC0, W0); 

		lrt = tracetest(DX,LXC,W);
		lra = LR(DX, LXC, W, sigma2);
		pt += (lrt[0].>LRt)/m;
		pa += (lra[0].>LRa)/m;
		//println(i);
	}	
	ranseed(seed);
	return double(pt)~double(pa);
}



pval_vb(const LRt, const LRa, const sigma,const psi, const sigma2, const DX0, const LXC0, const W0)
{
	// this implementation of wild boostrap trace test p value simulate the lag k dynamics
	// input sigma is n by 3 matrix	of the lower triangular version of the square root given by sig2tosig()
	decl T = rows(sigma);


	decl seed, p, m, DX, LXC, W, lrt,lra,pt,pa;
	seed = ranseed(0);
	ranseed(12345678);

	pt = 0;
	pa = 0;
	m = 499;

	decl eta,eps,X,i;
	for (i=0; i<m; ++i)
	{
		eta = rann(T,2);
		eps = sigma[][0].*eta[][0]~(sigma[][1].*eta[][1]+sigma[][2].*eta[][0]);
		[DX,LXC,W] = ECM0(psi, eps, DX0, LXC0, W0); 

		lrt = tracetest(DX,LXC,W);
		lra = LR(DX, LXC, W, sigma2);
		pt += (lrt[0].>LRt)/m;
		pa += (lra[0].>LRa)/m;		
		//println(i);
	}	
	ranseed(seed);
	return double(pt)~double(pa);
}

/* ================= top-level pipelines (per lag case) ================= */

pvals(const X, const k)
{
	//k >=2
	//prepare data
	decl n = rows(X);
	decl DX = X - lag0(X,1);
	decl LXC = lag0(X,1)~ones(n,1);  	
	decl W = lag0(DX,1);
	decl i;
	for (i=2; i<k; ++i)
		{W = W~lag0(DX,i);}
	decl W0 = W[k][];
	decl DX0 = DX[k][];
	decl LXC0 = LXC[k][];
	
	W = W[k:n-1][];
	DX = DX[k:n-1][];
	LXC = LXC[k:n-1][];
	decl T = n-k;
	
	decl lrt, psir1,er1,e1;
	
	[lrt,psir1,er1,e1] = tracetest(DX,LXC,W);

	decl h = cvbw_s(range(1,T), (e1[][0].^2~e1[][1].^2~e1[][0].*e1[][1])');
	decl sigma2hat = spcov_s((range(1,T))', (range(1,T))',(e1[][0].^2~e1[][1].^2~e1[][0].*e1[][1])', h);//output is 3 by n matrix
	decl sigma2hat_su = sigma2hat[0][]'**<1,0;0,0>+sigma2hat[1][]'**<0,0;0,1>+sigma2hat[2][]'**<0,1;1,0>;	  // stacked up version of sigma2hat, 2n by 2 matrix

	decl lra,psir2,er2;
	[lra,psir2,er2] = LR(DX, LXC, W, sigma2hat_su);
	
	decl sigma = sig2tosig(sigma2hat); // square root of the 3 by n matrix version, still give a 3 by n matrix of square root
									// used in the variance bootstrap function	
	decl pval = zeros(1,4);	
	pval[0:1] = pval_wb(lrt, lra, e1,psir1, sigma2hat_su, DX0,LXC0,W0);
	pval[2:3] = pval_vb(lrt, lra, sigma',psir1, sigma2hat_su, DX0,LXC0,W0);
	return {lrt, lra, pval};
}


pvalsk1(const X, const k)
{
	//k=1 case, when k=1, there is no W needed
	//prepare data
	decl n = rows(X);
	decl DX = X - lag0(X,1);
	decl LXC = lag0(X,1)~ones(n,1);
	
	decl DX0 = DX[k][];
	decl LXC0 = LXC[k][];

	DX = DX[k:n-1][];
	LXC = LXC[k:n-1][];
	
	decl T = n-k;
	
	decl lrt, psir1,er1,e1;
	
	[lrt,er1,e1] = tracetestk1(DX,LXC);

	decl h = cvbw_s(range(1,T), (e1[][0].^2~e1[][1].^2~e1[][0].*e1[][1])');
	decl sigma2hat = spcov_s((range(1,T))', (range(1,T))',(e1[][0].^2~e1[][1].^2~e1[][0].*e1[][1])', h);//output is 3 by n matrix
	decl sigma2hat_su = sigma2hat[0][]'**<1,0;0,0>+sigma2hat[1][]'**<0,0;0,1>+sigma2hat[2][]'**<0,1;1,0>;	  // stacked up version of sigma2hat, 2n by 2 matrix

	decl lra,psir2,er2;
	[lra,er2] = LRk1(DX, LXC, sigma2hat_su);
	
	decl sigma = sig2tosig(sigma2hat); // square root of the 3 by n matrix version, still give a 3 by n matrix of square root
									// used in the variance bootstrap function	
	decl pval = zeros(1,4);	
	pval[0:1] = pval_wbk1(lrt, lra, e1,k, sigma2hat_su,  LXC0);
	pval[2:3] = pval_vbk1(lrt, lra, sigma',k, sigma2hat_su, LXC0);
	return {lrt, lra, pval};
}



/* ===================== user-facing entry point ===================== */
coint_test(const X, const kmax)
{
    /* X    : n x 2 matrix, the two I(1) series.
       kmax : maximum lag for BIC order selection (>= 1).
       Returns a 1 x 7 row:
         [ k, trace, adaptiveLR,
           p_wildboot_trace, p_wildboot_adaptive,
           p_volboot_trace,  p_volboot_adaptive ]
       Reject the no-cointegration null at level alpha when the p-value < alpha.
       The adaptive test with the volatility bootstrap is the paper's recommended
       procedure (column p_volboot_adaptive). */
    decl k = icorder(X, kmax, 2);   /* BIC */
    decl lrt, lra, pval;
    if (k == 1)
        [lrt, lra, pval] = pvalsk1(X, k);
    else
        [lrt, lra, pval] = pvals(X, k);
    return double(k) ~ double(lrt) ~ double(lra) ~ pval;
}
