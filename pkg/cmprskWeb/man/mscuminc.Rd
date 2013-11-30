\name{mscuminc}
\alias{mscuminc}
\alias{plot.mscif}
\alias{print.mscif}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tests on the cumulative incidence function}
\description{Computes the cumulative incidence function (CIF) and performs tests of subgroup equality. 
It is a wrapper function for \code{cuminc()} in the package \code{cmprsk}. 
}
\usage{
mscuminc(msdata, cutpoint = NULL, censoring = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{msdata}{Data prepared with \code{msprep()}}
  \item{cutpoint}{Integer value which denotes the lower dwell time bound}
  \item{censoring}{If \code{TRUE} dwell times larger than cutpoint will be regarded as censored, if \code{FALSE} eliminated.}
}

\value{
  \item{msdata}{Input data after cutpoint treatment}
  \item{msres}{List with testing results for cumulative incidence function (from \code{cuminc}) for each page transition}
}

\seealso{\code{\link{msprep}}}
\examples{
data(webshopms)
msres <- mscuminc(webshopms, cutpoint = 1000, censoring = TRUE)
msres
plot(msres, from = 5.2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
