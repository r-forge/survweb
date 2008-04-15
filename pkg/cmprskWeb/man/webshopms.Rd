\name{webshopms}
\alias{webshopms}
\docType{data}
\title{Webshop data for multistage competing risks modeling}
\description{This data set for multistage computations is prepared with the function \code{msprep()}  
}
\usage{data(webshopms)}
\format{
  A data frame with 19836 observations on the following 7 variables.
  \describe{
    \item{\code{session}}{Session ID}
    \item{\code{user}}{User ID}
    \item{\code{from}}{Starting page (see details)}
    \item{\code{to}}{Landing page}
    \item{\code{time}}{Dwell time on starting page}
    \item{\code{nvisit}}{visit number within clickstream}
    \item{\code{buy}}{buyer vs. non-buyer}
  }
}
\details{The starting pages are prepared in terms of the number of visits within a 
clickstream, e.g. \code{"1.4"} denotes the transitions after the 4th visit of page 1.
}

\examples{
data(webshopms)
## maybe str(webshopms) ; plot(webshopms) ...
}
\keyword{datasets}
