\name{cmprskWeb-package}
\alias{cmprskWeb-package}
\alias{cmprskWeb}
\docType{package}
\title{Competing Risks Models for Web Data
}
\description{Provides various functions for competing risks analyses on Web data.
}
\details{
\tabular{ll}{
Package: \tab cmprskWeb\cr
Type: \tab Package\cr
Version: \tab 0.2-0\cr
Date: \tab 2013-11-01\cr
License: \tab GPL\cr
}

}
\author{Patrick Mair
}

\keyword{ package }

\examples{
data(webshopms)
msres <- mscuminc(webshopms, cutpoint = 1000, censoring = TRUE)
msres
plot(msres, from = 5.2)
}
