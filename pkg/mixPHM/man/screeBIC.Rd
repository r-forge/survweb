\name{screeBIC}
\alias{screeBIC}
\title{Scree plot of BIC's}
\description{
This function produces a scree plot on the basis of the BIC values in \code{msBIC}.
}

\usage{
screeBIC(x, lty = 1, col = NA, pch = 19, type = "b", main = "BIC Screeplot", 
xlab = "Number of Components", ylab = "BIC", legpos = "topright", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{mws} from \code{msBIC}}
  \item{lty}{Line type}
  \item{col}{Line colors; if \code{NA}, colors are determined automatically}
  \item{pch}{Value for plotting points}
  \item{type}{Type of plot}
  \item{main}{Plot title}
  \item{xlab}{Label for x-axis}
  \item{ylab}{Label for y-axis} 
  \item{legpos}{position of the legend} 
  \item{\dots}{Additional plot parameters}
}

\seealso{\code{\link{msBIC}}}
\examples{

##Fitting all Weibull proportional hazard models for K=2,3,4 components
data(webshop)
res <- msBIC(webshop, K = c(2,3,4), method = "all", maxiter = 5)
screeBIC(res)
}
\keyword{ hplot }
