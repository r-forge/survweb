\name{webshop}
\alias{webshop}
\docType{data}
\title{Webshop dataset for mixPHM package}
\description{This data set represents sojourn times in seconds of 333 sessions on 7 webpage 
categories of a webshop. Missing values indicate that the corresponding session did not visit a particular page.
}
\usage{data(webshop)}
\format{Numeric matrices of data frames with subjects as rows and variables as columns. 
Missing values are coded as \code{NA} or \code{0} (which corresponds to 0 survival time).
}

\source{
Mair, P., and Hudec, M. (2007). Analysis of dwell times in Web Usage Mining. 
Proceedings of the 31st Annual Conference of the German Classification Society 
on Data Analysis, Machine Learning, and Applications. 
}
\examples{
data(webshop)
str(webshop)
}
\keyword{datasets}
