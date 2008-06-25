\name{msBIC}
\alias{msBIC}
\alias{print.msBIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{PHM model selection with BIC}
\description{This function fits models for different proportionality restrictions.
}
\usage{
msBIC(x, K, method = "all", Sdist = "weibull", cutpoint = NULL, EMoption = "classification", EMstop = 0.01, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data frame or matrix of dimension n*p with survival times (\code{NA}'s allowed).}
  \item{K}{A vector with number of mixture components.}
  \item{method}{A vector with the methods provided in \code{phmclust}: 
  With \code{"separate"} no restrictions are imposed, \code{"main.g"} relates to a group main effect, 
  \code{"main.p"} to the variables main effects. \code{"main.gp"} reflects the proportionality assumption over groups
  and variables. \code{"int.gp"} allows for interactions between groups and variables. If method is \code{"all"}, each model is fitted.}
  \item{Sdist}{Various survival distrubtions such as \code{"weibull"}, \code{"exponential"}, and \code{"rayleigh"}.} 
  \item{cutpoint}{Cutpoint for censoring}
  \item{EMoption}{\code{"classification"} is based on deterministic cluster assignment, 
  \code{"maximization"} on deterministic assignment, and \code{"randomization"} 
  provides a posterior-based randomized cluster assignement.}  
  \item{EMstop}{Stopping criterion for EM-iteration.}
  \item{maxiter}{Maximum number of iterations.}
}
\details{Based on the output BIC matrix, model selection can be performed in terms of the
number of mixture components and imposed proportionality restrictions. 
}
\value{
  Returns an object of class \code{BICmat} with the following values:
  \item{BICmat}{Matrix with BIC values}
  \item{K}{Vector with different components}
  \item{method}{Vector with proportional hazard methods}
  \item{Sdist}{Survival distribution}
}

\seealso{\code{\link{screeBIC}}}

\examples{

##Fitting 3 Weibull proportional hazard models (over groups, pages) for K=2,3 components
data(webshop)
res <- msBIC(webshop, K = c(2,3), method = c("main.p","main.g"), maxiter = 10)
res 
}

\keyword{ models }
