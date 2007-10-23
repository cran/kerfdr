kerfdr = function (
  pv, x = NULL, f0 = NULL, localfdr = NULL, pi1 ="storey",
  lambda=0.5, bw = "sj-dpi", adjust = 1.0,
  trans = c("probit", "log", "none"),
  kernel = c("gaussian","epanechnikov", "rectangular", "triangular", "biweight","cosine", "optcosine"),
  pvMin = 0.0, pvMax = 1.0, verbose = FALSE, plot = TRUE)  {
        # general parameters
  n = 512;
  give.Rkern = FALSE;
  cut = 3;
  tol = 1e-15;
  maxiter = 2000;

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
  truncinf=(pv<=pvMin);
  truncsup=(pv>=pvMax);
  trunc=(!truncinf & !truncsup);
	nx=sum(trunc);
	if (nx==0)
		stop(paste("no p-values in the valid domains: ]",pvMin, " ; ", pvMax, " ["));
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
	tau[!working]=1-localfdr[!working];

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
  p0inf=pvMin;
  p0sup=1-pvMax;
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
		if (verbose) 
			print(paste("iter =",round(iter),"; test =",test))
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

	# results as a density object
    	results = list(pv = pv, x = x, f = pi1*f1+(1-pi1)*f0,
		f0 = f0, f1 = f1,pi0 = pi0, pi1 = pi1, p0 = p0, p1 = p1, localfdr = localfdr,
		bw = bw, iter = iter, call = match.call(), data.name = name,
		kords = kords)

	# graphic output
	if (plot == TRUE)
	{
		xlab=switch(trans,probit="probit(pv)",log="log10(pv)",none="pv");
		main=paste("kerfdr(): pi1 = ",round(1000*pi1)/1000," and bw = ",round(1000*bw)/1000);
		par(mfrow = c(2,1))
		hist(x,freq = FALSE, col = "lightgrey", border= "darkgrey",
			xlab = xlab, ylab = "density", main = main, ylim=c(0,max(results$f[trunc])),
			nc=max(20,round(diff(range(xtrunc))/bw)))
		points(x,f0*pi0, pch = 18, cex = 0.5, col = "darkblue")
		points(x,f1*pi1, pch = 18, cex = 0.5, col = "red")
		I=sort(xtrunc,index.return=TRUE)$ix;
		lines(xtrunc[I],results$f[trunc][I], col = "black")
		plot(range(xtrunc),c(0,1.0),t="n",xlab=xlab, ylab="local fdr");
		points(x,localfdr, cex = 0.5, pch = 20)
		
	}
	
	return(results)
}

