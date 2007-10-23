# Rd
# description >> This function computes local fdr values by using a two-components mixture model with a semi-parametric density estimation. The code is freely inspired from the \link[stats]{density} function. For a simple use, we recommand the default setting (most parameters are optional). 
# argument
# item >> pv >> the vector of raw p-values.
# item >> x >> a transformation of \code{pv}. It can be given by the user or (if \code{NULL}) computed via the \code{trans} parameter    
# item >> trans >> the transformation to apply on \code{pv} to produce \code{x}: \code{"probit"} (by default) returns \code{qnorm(pv)} and \code{"log"} returns \code{log10(pv)}. 
# item >> f0 >> the sample density under the null hypothesis. Can be specified by the user. If \code{NULL} (by default) the density under H0 is determined according to \code{trans}: if \code{trans = "probit"} then \code{f0} is a standard Gaussian distribution; if \code{trans = "log"} then \code{f0} is a standard Exponential distribution; if \code{trans = "none"} then \code{f0} is a standard Uniform distribution
# item >> localfdr >> values to initiate the iterative algorithm. If \code{NULL} (by default) initial values are then sampled in a Uniform distribution [0,1]
# item >> pi1 >> a priori proportion of alternative hypothesis or a method (string) to compute it; by default it uses the method proposed by Storey and Tibshirani (2003).
# item >> lambda >> p-value threshold for the Storey's calculation of \code{pi1} (0.5 by default). See \link[qvalue]{qvalue} for more details. 
# item >> bw >> a bandwidth value or a method to determine it among \code{"nrd0", "nrd", "ucv", "bcv", "sj-ste", "sj-dpi"}. See \link[stats]{bandwidth} for more details. 
# item >> kernel >> the kernel used (string) among \code{"gaussian"} (by default), \code{"epanechnikov", "rectangular", "triangular", "biweight","cosine"}. For more details on kernels: \url{http://stat.genopole.cnrs.fr/sg/software/kerfdr/kernels}
# item >> truncat >> an interval on p-values to deal with truncated distributions such as those obtained with Monte-Carlo simulations.   
# item >> plot >> if \code{TRUE}, it returns graphics of local fdr estimations. Some plots are inspired from \link[qvalue]{qvalue}. 
# item >> cuts >> vector of significance values to use in \code{summary} (see below)
# value >> A list of parameters (\code{pv, x, pi1, bw, f0} ...) and the following results: 
# item >> f >> the observed mixture density
# item >> f1 >> the estimated density under H1
# item >> localfdr >> the local fdr values resulting from the algorithm 
# item >> summary >> a summary table comparing the number of significant calls for the raw p-values, Bonferroni and Benjamini-Hochberg corrections and for the calculated local fdr, using a set of cutoffs given by \code{cuts}
# author >> M Guedj, A Celisse, S Robin, G Nuel
# keyword >> nonparametric
# references >> \url{http://stat.genopole.cnrs.fr/sg/software/kerfdr}, Robin et al (2007), Strimmer (2008), Guedj et al (2009)
# examples >> # Example 1: kerfdr with different plots >> generating a vector of p-values with pi0 = 80%
# examples >> n = 10000
# examples >> pi0 = 0.8
# examples >> # plot in a probit scale (default)
# examples >> pv = 1-pnorm(c(rnorm(n*pi0), rnorm(n*(1-pi0), 4)))
# examples >> res = kerfdr(pv)
# examples >> res$pi0
# examples >> res$summary
# examples >> # plot in a log scale
# examples >> kerfdr(pv, trans = "log")
# examples >> # plot in the raw p-values scale
# examples >> kerfdr(pv, trans = "none")
# examples >> # Example 2: truncation on a vector of null p-values (resulting local fdr should be 1 for each point)
# examples >> n = 10000
# examples >> pv = runif(n)
# examples >> # truncation on [0.1;0.9] 
# examples >> pv[which(pv < 0.1)] = 0.1
# examples >> pv[which(pv > 0.9)] = 0.9
# examples >> # kerfdr WITHOUT taking the truncation into account (local fdr is hence badly estimated)
# examples >> kerfdr(pv, trans = "log")
# examples >> # kerfdr by taking the truncation into account (local fdr is then well estimated)
# examples >> kerfdr(pv, truncat = c(0.1, 0.9), trans = "log")
# end
kerfdr = function(
  pv, x = NULL, trans = c("probit", "log", "none"), f0 = NULL, localfdr = NULL, pi1 ="storey",
  lambda=0.5, bw = c("nrd0", "nrd", "ucv", "bcv", "sj-ste", "sj-dpi"),
  kernel = c("gaussian","epanechnikov", "rectangular", "triangular", "biweight","cosine"),
  truncat = c(0,1), plot = TRUE, cuts = c(1e-04, 0.001, 0.01, 0.025, 0.05, 0.1))  {
        # general parameters
  n = 512;
  adjust = 1.0
  give.Rkern = FALSE;
  cut = 3;
  tol = 1e-15;
  maxiter = 2000;
  bw = match.arg(bw)
	# parse options
	kernel = match.arg(kernel)
	trans = match.arg(trans)
	
	if (missing(pv))
        	stop("argument 'pv' is required")
	if (!is.numeric(pv))
        	stop("argument 'pv' must be numeric")
	if (is.null(x))
	   x = switch(tolower(trans), probit = qnorm(pv), log = log10(pv), none=pv);
	if (is.null(f0))
	   f0 = switch(tolower(trans), probit = dnorm(x), log = log(10)*dexp(-log(10)*x), none=dunif(x)); 
	if (give.Rkern)
        	return(switch(kernel, gaussian = 1/(2 * sqrt(pi)), rectangular = sqrt(3)/6, triangular = sqrt(6)/9, epanechnikov = 3/(5 * sqrt(5)), biweight = 5 * sqrt(7)/49, cosine = 3/4 * sqrt(1/3 - 2/pi^2), optcosine = sqrt(1 - 8/pi^2) * pi^2/16))
    	name = deparse(substitute(x))
    	x = as.vector(x)
    		
	# truncature if needed
  truncinf=(pv<=truncat[1]);
  truncsup=(pv>=truncat[2]);
  trunc=(!truncinf & !truncsup);
	nx=sum(trunc);
	if (nx==0)
		stop(paste("no p-values in the valid domains: ]",truncat[1], " ; ", truncat[2], " ["));
	xtrunc=x[trunc];
  
	# get working index
	if (is.null(localfdr)) {
	   localfdr=rep(NA,length(pv));
	}
	working=is.na(localfdr);
	
			
	# determination of bw according to the selected method
	if (is.character(bw)) {
	   if (nx < 2)
		stop("need at least 2 points to select a bandwidth automatically")
	   bw = switch(tolower(bw), nrd0 = bw.nrd0(xtrunc), nrd = bw.nrd(xtrunc),
		ucv = bw.ucv(xtrunc), bcv = bw.bcv(xtrunc), "sj-ste" = bw.SJ(xtrunc, method = "ste"),
		"sj-dpi" = bw.SJ(xtrunc, method = "dpi"), stop("unknown bandwidth rule"))
	}
	if (!is.finite(bw))
	   stop("non-finite 'bw'")
	bw = adjust * bw
	if (bw <= 0)
	   stop("'bw' is not positive.")
	
	# determination of pi1 according to the selected method
	if (is.character(pi1))
	   pi1= switch(tolower(pi1), storey=1-mean(pv[working]>lambda)/(1.0-lambda),stop("unknown pi1 rule"))	
	if (pi1<0.0) {
	   warning(paste("estimated pi1 =",round(pi1,digit=4),"set to 0.0"));
	   pi1=0.0;
	}
 	if (pi1>1.0) {
	   warning(paste("estimated pi1 =",round(pi1,digit=4),"set to 1.0"));
	   pi1=1.0;
	}
	if (sum(!working)!=0) {
	   pi1 = sum(working)/length(pv)*pi1 + sum(!working)/length(pv)*(1-mean(localfdr[!working]));	
	}
	# set pi0 from pi1
	pi0 = 1-pi1    	

	# random initialization of tau = 1.0 - localfdr
	tau=runif(x);
	tau[!working]=1.0-localfdr[!working];

	# compute working interval		
	from = min(xtrunc) - cut * bw
	to = max(xtrunc) + cut * bw
	if (!is.finite(from))
		stop("non-finite 'from'")
	if (!is.finite(to))
		stop("non-finite 'to'")
	lo = from - 4 * bw
	up = to + 4 * bw

	# compute truncature conditional probabilities
  p0inf=truncat[1];
  p0sup=1-truncat[2];
  if (pi1>0) {
    p1inf=(sum(truncinf)/length(pv)-pi0*p0inf)/pi1;
    p1sup=(sum(truncsup)/length(pv)-pi0*p0sup)/pi1;
  } else {
    warning(paste("pi1=0.0 hence both p1inf and p1sup are set to 0.0"));
    p1inf=0.0;
    p1sup=0.0;
  }
  if (p1inf<0.0) {
    warning(paste("estimated p1inf =",p1inf,"set to 0.0"));
    p1inf=0.0;
  }
  if (p1sup<0.0) {
    warning(paste("estimated p1sup =",p1sup,"set to 0.0"));
    p1sup=0.0;
  }
  if (p1inf>1.0) {
    warning(paste("estimated p1inf =",p1inf,"set to 1.0"));
    p1inf=1.0;
  }
  if (p1sup>1.0) {
    warning(paste("estimated p1sup =",p1sup,"set to 1.0"));
    p1sup=1.0;
  }  
  p0=1-p0inf-p0sup;
  p1=1-p1inf-p1sup;
 	# build kernel
	kords = seq(0, 2 * (up - lo), length = 2 * n)
	kords[(n + 2):(2 * n)] = -kords[n:2]
	kords = switch(kernel, gaussian = dnorm(kords, sd = bw),
		rectangular = {a = bw * sqrt(3); ifelse(abs(kords) < a, 0.5/a, 0)},
		triangular = { a = bw * sqrt(6); ax = abs(kords); ifelse(ax < a, (1 - ax/a)/a, 0)},
		epanechnikov = {a = bw*sqrt(5); ax = abs(kords); ifelse(ax < a, 3/4 * (1 - (ax/a)^2)/a, 0)},
		biweight = {a = bw * sqrt(7); ax = abs(kords); ifelse(ax < a, 15/16 * (1 - (ax/a)^2)^2/a, 0)},
		cosine = {a = bw/sqrt(1/3 - 2/pi^2); ifelse(abs(kords) < a, (1 + cos(pi*kords/a))/(2*a), 0)},
		optcosine = {a = bw/sqrt(1 - 8/pi^2); ifelse(abs(kords) < a, pi/4 * cos(pi*kords/(2*a))/a, 0)})
	xords = seq(lo, up, length = n)


	# main loop
	converged = FALSE
	iter = 0
	newtau=tau;
	f1=f0;
	while (!converged && iter < maxiter) {
		iter = iter+1;
		# mass distribution according to tau
		y = .C("massdist", x = as.double(xtrunc),
			xmass = as.double(tau[trunc]/sum(tau[trunc])),
			nx = nx, xlo = as.double(lo), xhi = as.double(up),
			y = double(2*n), ny = as.integer(n), PACKAGE = "stats")$y
		# fft convolution
		conv = fft(fft(y)*Conj(fft(kords)), inv = TRUE)
		conv = Re(conv)[1:n]/length(y)
		# f1 update
		f1[trunc] = approx(xords,conv,xtrunc)$y
		# localfdr update
		newtau[trunc & working] = pi1*p1*f1[trunc & working]/(pi1*p1*f1[trunc & working]+ pi0*p0*f0[trunc & working])
		# check for convergence
		test = max(abs(tau-newtau))
		if (test < tol) 
			converged = TRUE
		tau = newtau
		if (sum(tau[trunc])==0.0) {
		   converged= TRUE;
		}
	}
	    
	# gestion of the truncature if needed. 
	f0[!trunc]=NA;
	f1[!trunc]=NA;
	localfdr[trunc & working]=1.0-tau[trunc & working];
	localfdr[truncinf & working]=p0inf*pi0/(p0inf*pi0+p1inf*pi1);
	localfdr[truncsup & working]=p0sup*pi0/(p0sup*pi0+p1sup*pi1);
	x[truncinf]=min(xtrunc);
	x[truncsup]=max(xtrunc);

  bon = pv*length(pv)
  bon[which(bon > 1)] = 1
  bh = pv*length(pv)/rank(pv)
  bh[which(bh > 1)] = 1
  
  counts <- sapply(cuts, function(x) c("p-value" = sum(pv <= x), "FWER (Bonf.)"= sum(bon <= x), "FDR (BH)"= sum(bh <= x),"kerfdr" = sum(localfdr <= x)))
  colnames(counts) <- paste("<", cuts, sep = "")


	# results as a density object
    	results = list(pv = pv, x = x, f = pi1*f1+(1-pi1)*f0,
		f0 = f0, f1 = f1,pi0 = pi0, pi1 = pi1, p0 = p0, p1 = p1,
		bw = bw, iter = iter, call = match.call(), data.name = name,
		kords = kords, localfdr = localfdr, summary = counts)

	results = list(pv = pv, x = x, f = pi1*f1+(1-pi1)*f0,
		f0 = f0, f1 = f1,pi0 = pi0, pi1 = pi1, p0 = p0, p1 = p1,
		bw = bw, iter = iter, data.name = name,
		kords = kords, localfdr = localfdr, summary = counts)

	# graphic output
	if (plot == TRUE)
	{
	  layout(matrix(c(1,1, 2, 2,3,4), 3,2, byrow = TRUE))
		xlab=switch(trans,probit="probit(pv)",log="log10(pv)",none="pv");
		main=paste("kerfdr(): pi1 = ",round(1000*pi1)/1000," and bw = ",round(1000*bw)/1000);
		#par(mfrow = c(2,1))
		hist(x,freq = FALSE, col = "lightgrey", border= "darkgrey",
			xlab = xlab, ylab = "density", main = main, ylim=c(0,max(results$f[trunc])),
			nc=max(20,round(diff(range(xtrunc))/bw)))
		points(x,f0*pi0, pch = 18, cex = 0.5, col = "darkblue")
		points(x,f1*pi1, pch = 18, cex = 0.5, col = "red")
		I=sort(xtrunc,index.return=TRUE)$ix;
		lines(xtrunc[I],results$f[trunc][I], col = "black")
		plot(range(xtrunc),c(0,1.0),t="n",xlab=xlab, ylab="local fdr");
		points(x,localfdr, cex = 0.5, pch = 20)
		plot(pv, localfdr,xlab = "p-value", ylab = "local fdr", xlim = c(0,1), ylim = c(0,1), col = "darkgrey", cex = 0.7) 
		plot(sort(localfdr), 1:length(sort(localfdr)),xlab = "local fdr", ylab = "nb of tests", xlim = c(0,1), col = "darkgrey", cex = 0.7) 	
	}
	
	return(results)
}

