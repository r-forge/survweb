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
Version: \tab 0.1.2\cr
Date: \tab 2008-04-08\cr
License: \tab GPL\cr
}

}
\author{Patrick Mair
}
\references{Mair, P. and Hudec, M. (2008). Competing Risks Approaches for Modeling Dwell Time Based Page
Transitions in Web Usage Mining. Under editorial review in Journal of Data Mining and Knowledge Discovery. 
}
\keyword{ package }
\seealso{\code{\link[cmprsk:cmprsk-package]{cmprsk}} ~~
}
\examples{
data(webshopms)
library("cmprsk")
msres <- mscuminc(webshopms, cutpoint = 1000, censoring = TRUE)
msres
plot(msres, from = 5.2)
}
