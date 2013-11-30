\name{linkfreq}
\alias{linkfreq}
\alias{print.linkfreq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transition frequencies}
\description{Utility function for determining the transition frequencies.
}
\usage{
linkfreq(from, to, linkcut = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{from}{vector with starting pages}
  \item{to}{vector with landing pages}
  \item{linkcut}{Cutpoint for excluding links (see details)}
}
\details{A value of \code{linkcut = 0.001} means that links used in more than 0.1\% of the transitions 
are taken into account for frequency computation.
}
\value{Returns an object of class \code{"linkfreq"}
  \item{freqlink}{Transition frequencies as data frame.}
  \item{linktable}{Transition frequencies as table.}
  \item{fromSum}{Vector with starting page margins.}
  \item{toSum}{Vector with landing page margins.}
  ...
}

\seealso{\code{\link{msprep}}}
\examples{
data(webshop)
attach(webshop)
res <- linkfreq(from, to, linkcut = 0.01)
res
}

\keyword{manip}
