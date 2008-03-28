\name{stableEM}
\alias{stableEM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stable EM solution}
\description{
This function performs the clustering for different EM starting values in order
to find a stable solution.
}
\usage{
stableEM(x, K, numEMstart = 5, method = "separate", Sdist = "weibull", 
EMoption = "classification", EMstop = 0.0001, maxiter = 1000, print.likvec = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data frame or matrix of dimension n*p with survival times (\code{NA}'s allowed).}
  \item{K}{Number of mixture components.}
  \item{numEMstart}{Number of different starting solutions}
  \item{method}{Imposing proportionality restrictions on the hazards: 
  With \code{separate} no restrictions are imposed, \code{main.g} relates to a group main effect, 
  \code{main.p} to the variables main effects. \code{main.gp} reflects the proportionality assumption over groups
  and variables. \code{int.gp} allows for interactions between groups and variables.}
  \item{Sdist}{Various survival distrubtions such as \code{weibull}, \code{exponential}, and \code{rayleigh}.} 
  \item{EMoption}{\code{classification} is based on deterministic cluster assignment, 
  \code{maximization} on deterministic assignment, and \code{randomization} 
  provides a posterior-based randomized cluster assignement.}  
  \item{EMstop}{Stopping criterion for EM-iteration.}
  \item{maxiter}{Maximum number of iterations.}
  \item{print.likvec}{If \code{TRUE} the likelihood values for different starting solutions are printed.}
}

\details{
After the computation of the models for different starting solutions using the
function \code{phmclust} the best model is chosen, i.e., the model with the largest likelihood value. 
The output values refer to this final model.
}

\value{
  Returns an object of class \code{mws} with the following values:
  \item{K}{Number of components}
  \item{iter}{Number of EM iterations}
  \item{method}{Method with propotionality restrictions used for estimation}
  \item{Sdist}{Assumed survival distribution}
  \item{likelihood}{Log-likelihood value for each iteration}
  \item{pvisit}{Matrix of prior probabilities due to \code{NA structure}}
  \item{shape}{Matrix with shape parameters}
  \item{scale}{Matrix with scale parameters}
  \item{group}{Final deterministic cluster assignment}
  \item{posteriors}{Final probabilistic cluster assignment}
  \item{npar}{Number of estimated parameters}
  \item{aic}{Akaike information criterion}
  \item{bic}{Bayes information criterion}
  \item{clmean}{Matrix with cluster means}
  \item{clmed}{Matrix with cluster medians}
}


\seealso{\code{\link{phmclust}},\code{\link{msBIC}}}
\examples{

## Exponental mixture model with 2 components for 4 different starting solutions
data(webshop)
res <- stableEM(webshop, K = 2, numEMstart = 4, Sdist = "exponential")
res
summary(res)
}

\keyword{ models }
