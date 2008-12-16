\name{phmclust}
\alias{phmclust}
\alias{summary.mws}
\alias{print.mws}
\alias{logLik.mws}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fits mixtures of proportional hazard models}
\description{This function allows for the computation of proportional hazards models with different distribution assumptions 
on the underlying baseline hazard. Several options for imposing proportionality restrictions on the hazards are provided.
This function offers several variations of the EM-algorithm regarding the posterior computation in the M-step.
}
\usage{
phmclust(x, K, method = "separate", Sdist = "weibull", cutpoint = NULL, EMstart = NA, 
EMoption = "classification", EMstop = 0.01, maxiter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data frame or matrix of dimension n*p with survival times (\code{NA}'s allowed).}
  \item{K}{Number of mixture components.} 
  \item{method}{Imposing proportionality restrictions on the hazards: 
  With \code{"separate"} no restrictions are imposed, \code{"main.g"} relates to a group main effect, 
  \code{"main.p"} to variable main effects. \code{"main.gp"} reflects the proportionality assumption over groups
  and variables. \code{"int.gp"} allows for interactions between groups and variables.} 
  \item{Sdist}{Various survival distrubtions such as \code{"weibull"}, \code{"exponential"}, and \code{"rayleigh"}.} 
  \item{cutpoint}{Integer value with upper bound for observed dwell times. Above this cutpoint, values are regarded as censored. If NULL, no censoring is performed}
  \item{EMstart}{Vector of length n with starting values for group membership, 
  \code{NA} indicates random starting values.}
  \item{EMoption}{\code{"classification"} is based on deterministic cluster assignment, 
  \code{"maximization"} on deterministic assignment, and \code{"randomization"} 
  provides a posterior-based randomized cluster assignement.}  
  \item{EMstop}{Stopping criterion for EM-iteration.}
  \item{maxiter}{Maximum number of iterations.}
}
\details{The method \code{"separate"} corresponds to an ordinary mixture model. \code{"main.g"} imposes proportionality 
restrictions over variables (i.e., the group main effect allows for free-varying variable hazards). \code{"main.p"} 
imposes proportionality restrictions over groups (i.e., the variable main effect allows for free-varying group hazards).
If clusters with only one observation are generated, the algorithm stops. 
}
\value{
  Returns an object of class \code{mws} with the following values:
  \item{K}{Number of components}
  \item{iter}{Number of EM iterations}
  \item{method}{Proportionality restrictions used for estimation}
  \item{Sdist}{Assumed survival distribution}
  \item{likelihood}{Log-likelihood value for each iteration}
  \item{pvisit}{Matrix of prior probabilities due to \code{NA structure}}
  \item{se.pvisit}{Standard errors for priors}
  \item{shape}{Matrix with shape parameters}
  \item{scale}{Matrix with scale parameters}
  \item{group}{Final deterministic cluster assignment}
  \item{posteriors}{Final probabilistic cluster assignment}
  \item{npar}{Number of estimated parameters}
  \item{aic}{Akaike information criterion}
  \item{bic}{Bayes information criterion}
  \item{clmean}{Matrix with cluster means}
  \item{se.clmean}{Standard errors for cluster means}
  \item{clmed}{Matrix with cluster medians}
}
\references{
Mair, P., and Hudec, M. (2008). Analysis of dwell times in Web Usage Mining. 
Proceedings of the 31st Annual Conference of the German Classification Society 
on Data Analysis, Machine Learning, and Applications.

Collett, D. (2003). Modelling Survival Data in Medical Research. Boca Raton, FL: Chapman & Hall.

Celaux, G., and Govaert, G. (1992). A classification EM algorithm for clustering and
two stochastic versions. Computational Statistics and Data Analysis, 14, 315-332.
}

\seealso{\code{\link{stableEM}}, \code{\link{msBIC}} }
\examples{

data(webshop)

## Fitting a Weibll mixture model (3 components) is fitted with classification EM 
## Observations above 600sec are regarded as censored

res1 <- phmclust(webshop, K = 3, cutpoint = 600)
res1
summary(res1)

## Fitting a Rayleigh Weibull proportional hazard model (2 components, proportional over groups)
res2 <- phmclust(webshop, K = 2, method = "main.p", Sdist = "rayleigh") 
res2
summary(res2)

}

\keyword{ models }

