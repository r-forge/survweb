\name{webshop}
\alias{webshop}
\docType{data}
\title{Webshop Data}
\description{Artificial webshop dataset with page dwell times.}
\usage{data(webshop)}
\format{
  A data frame with 20000 observations on the following 6 variables.
  \describe{
    \item{\code{session}}{Session ID}
    \item{\code{user}}{User ID}
    \item{\code{from}}{Starting page}
    \item{\code{to}}{Landing page}
    \item{\code{time}}{Dwell time on starting page}
    \item{\code{buy}}{buyer vs. non-buyer}
  }
}

\examples{
data(webshop)
}
\keyword{datasets}
