\name{kerfdr}
\alias{kerfdr}
\title{kerfdr: a semi-parametric kernel-based approach to local fdr estimations}
\description{this function computes local fdr values by using a classical two-components mixture model with a semi-parametric density estimation. The code is freely inspired from the \code{density} function.} 
\usage{
kerfdr(pv, x = NULL, f0 = NULL, localfdr = NULL, pi1 = "storey", lambda = 0.5, bw = "sj-dpi", adjust = 1.0, trans = c("probit", "log", "none"), kernel = c("gaussian","epanechnikov", "rectangular", "triangular", "biweight","cosine", "optcosine"), pvMin = 0.0, pvMax = 1.0, verbose = FALSE, plot = TRUE)
}
\arguments{
  \item{pv}{a vector of p-values}
  \item{trans}{the transformation to apply on \code{pv} to produce \code{x}}
  \item{x}{transformation of pv, \code{qnorm(pv)} with "probit" (by default), \code{log10(pv)} with "log"}
  \item{f0}{the sample density under the null hypothesis}
  \item{localfdr}{initial values for the local fdr}
  \item{pi1}{proportion of alternative hypothesis (float) or a method (string) to compute it; by default it uses the method introduced by Storey}
  \item{bw}{a bandwidth value (float) or a method (string) to determine it}
  \item{lambda}{p-value threshold for the Storey's calculation of \code{pi1}; 0.5 by default}
  \item{adjust}{the bandwidth used is \code{adjust*bw}; 1.0 by default}
  \item{kernel}{the kernel to use (string); "gaussian" by default}
  \item{plot}{if TRUE, it produces a graphical display of the local fdr estimations}
  \item{pvMin}{the p-value under which to take a truncature in consideration; 0.0 by default}
  \item{pvMax}{the p-value upper which to take a truncature in consideration; 1.0 by default}
  \item{verbose}{if TRUE, it activates the verbose mode; FALSE by default}
}
\value{kerfdr returns a list including \code{pv}, \code{x}, \code{pi1}, \code{pi0}, \code{bw} and the following components:
  \item{localfdr}{local fdr estimations}
  \item{f0}{the estimated sample density under the null hypothesis}
  \item{f1}{the estimated sample density under the alternative hypothesis}
  \item{f}{the resulting mixed density}
}
\references{\url{http://stat.genopole.cnrs.fr/sg/software/kerfdr}}
\author{Mickael Guedj, Gregory Nuel}
\keyword{nonparametric}
\keyword{models}
\examples{
  # p-values under the null hypothesis:
  pv0 = runif(10000)
  # p-values under the alternative hypothesis:
  pv1 = runif(500, 0,0.05)
  # local fdr estimations:
  myPv = c(pv0, pv1)
  kerfdr(pv = myPv)
}







