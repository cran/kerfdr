\name{kerfdr}
\alias{kerfdr}
\title{kerfdr}
\description{This function computes local fdr values by using a two-components mixture model with a semi-parametric density estimation. The code is freely inspired from the \link[stats]{density} function. For a simple use, we recommand the default setting (most parameters are optional). }
\usage{kerfdr(pv, x=NULL, trans=c("probit", "log", "none"), f0=NULL,localfdr=NULL, pi1="storey", lambda=seq(0, 0.9, 0.05),bw=c("nrd0","nrd", "ucv", "bcv", "sj-ste", "sj-dpi"), kernel=c("gaussian","epanechnikov", "rectangular", "triangular", "biweight","cosine"),truncat=c(0, 1), plot=TRUE, cuts=c(1e-04, 0.001, 0.01, 0.025, 0.05, 0.1)) }
\arguments{
\item{pv}{the vector of raw p-values.}
\item{x}{a transformation of \code{pv}. It can be given by the user or (if \code{NULL}) computed via the \code{trans} parameter    }
\item{trans}{the transformation to apply on \code{pv} to produce \code{x}: \code{"probit"} (by default) returns \code{qnorm(pv)} and \code{"log"} returns \code{log10(pv)}. }
\item{f0}{the sample density under the null hypothesis. Can be specified by the user. If \code{NULL} (by default) the density under H0 is determined according to \code{trans}: if \code{trans = "probit"} then \code{f0} is a standard Gaussian distribution; if \code{trans = "log"} then \code{f0} is a standard Exponential distribution; if \code{trans = "none"} then \code{f0} is a standard Uniform distribution}
\item{localfdr}{values to initiate the iterative algorithm. If \code{NULL} (by default) initial values are then sampled in a Uniform distribution [0,1]}
\item{pi1}{a priori proportion of alternative hypothesis or a method (string) to compute it; by default it uses the method proposed by Storey and Tibshirani (2003).}
\item{lambda}{p-value threshold for the Storey's calculation of \code{pi1} (0.5 by default). See \link[qvalue]{qvalue} for more details. }
\item{bw}{a bandwidth value or a method to determine it among \code{"nrd0", "nrd", "ucv", "bcv", "sj-ste", "sj-dpi"}. See \link[stats]{bandwidth} for more details. }
\item{kernel}{the kernel used (string) among \code{"gaussian"} (by default), \code{"epanechnikov", "rectangular", "triangular", "biweight","cosine"}. For more details on kernels: \url{http://stat.genopole.cnrs.fr/sg/software/kerfdr/kernels}}
\item{truncat}{an interval on p-values to deal with truncated distributions such as those obtained with Monte-Carlo simulations.   }
\item{plot}{if \code{TRUE}, it returns graphics of local fdr estimations. Some plots are inspired from \link[qvalue]{qvalue}. }
\item{cuts}{vector of significance values to use in \code{summary} (see below)}
}
\value{A list of parameters (\code{pv, x, pi1, bw, f0} ...) and the following results: 
\item{f}{the observed mixture density}
\item{f1}{the estimated density under H1}
\item{localfdr}{the local fdr values resulting from the algorithm }
\item{summary}{a summary table comparing the number of significant calls for the raw p-values, Bonferroni and Benjamini-Hochberg corrections and for the calculated local fdr, using a set of cutoffs given by \code{cuts}}
}
\author{M Guedj, A Celisse, ML Martin-Migniette, S Robin, G Nuel}
\keyword{nonparametric}
\references{\url{http://stat.genopole.cnrs.fr/sg/software/kerfdr}, Robin et al (2007), Strimmer (2008), Guedj et al (2009)}
\examples{# Example 1: kerfdr with different plots
n = 10000
pi0 = 0.8
# plot in a probit scale (default)
pv = 1-pnorm(c(rnorm(n*pi0), rnorm(n*(1-pi0), 4)))
res = kerfdr(pv)
res$pi0
res$summary
# plot in a log scale
kerfdr(pv, trans = "log")
# plot in the raw p-values scale
kerfdr(pv, trans = "none")
# Example 2: truncation on a vector of null p-values (resulting local fdr should be 1 for each point)
n = 10000
pv = runif(n)
# truncation on [0.1;0.9] 
pv[which(pv < 0.1)] = 0.1
pv[which(pv > 0.9)] = 0.9
# kerfdr WITHOUT taking the truncation into account (local fdr is hence badly estimated)
kerfdr(pv, trans = "log")
# kerfdr by taking the truncation into account (local fdr is then well estimated)
kerfdr(pv, truncat = c(0.1, 0.9), trans = "log") }
