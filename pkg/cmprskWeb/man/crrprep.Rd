\name{crrprep}
\alias{crrprep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Competing risks data preparation}
\description{Adds covariates for clickstream position and number of page visits to the dataset and prepars censoring.
}
\usage{
crrprep(linkdata, from, cutpoint = NULL, censoring = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{linkdata}{Data frame with session ID, user ID, from, to, dwelltime, covariate}
  \item{from}{Starting page for competing risks regression}
  \item{cutpoint}{Upper bound for page dwell time}
  \item{censoring}{If \code{TRUE} PI's with dwell times larger than \code{cutpoint} are regarded as censored.If \code{FALSE}, 
  these observations are deleted}
}

\value{Returns a data frame suited for competing risks regression.
}

\seealso{\code{\link{msprep}}}
\examples{
data(webshop)
crdata <- crrprep(webshop, from = 2, cutpoint = 300)
head(crdata)
}
\keyword{manip}
