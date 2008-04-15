\name{msprep}
\alias{msprep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multistage data preparation}
\description{Data preparation utility function for taking into account the number of visits (within
the clickstream) of the starting pages.
}
\usage{
msprep(linkdata, linkfreq.res, stage.depth = 2, elim = TRUE, sess = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{linkdata}{Data frame with session ID, user ID, from, to, dwelltime, covariate}
  \item{linkfreq.res}{Result from \code{linkfreq()} analysis}
  \item{stage.depth}{Integer denoting the stage depth}
  \item{elim}{It \code{TRUE} it eliminates PI's larger than stage depth, if \code{FALSE}, 
  successive PI's go into last stage}
  \item{sess}{If \code{TRUE}, analysis is performed at session level, otherwise at user level}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{Returns a data frame suited for multistage competing risks modeling.
}

\seealso{\code{\link{mscuminc}}}
\examples{
data(webshop)
webshopsub <- webshop[1:1000,]
library("gtools")
freqsub <- linkfreq(webshopsub[,3], webshopsub[,4])
msdata <- msprep(webshopsub, freqsub, stage.depth = 4)
}
\keyword{manip}
