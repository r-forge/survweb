\name{webshop}
\alias{webshop}
\docType{data}
\title{Webshop dataset for mixPHM package}
\description{This artificial data set represents dwell times in seconds of 333 sessions on 7 webpage 
categories of a webshop. Missing values indicate that the corresponding session did not visit a particular page.
}
\usage{data(webshop)}
\format{Numeric matrices of data frames with subjects as rows and variables as columns. 
Missing values are coded as \code{NA} (which corresponds to 0 survival time).
}

\examples{
data(webshop)
str(webshop)
}
\keyword{datasets}
