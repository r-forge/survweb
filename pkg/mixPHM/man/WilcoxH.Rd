\name{WilcoxH}
\alias{WilcoxH}
\alias{print.wilcoxh}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tests of Zero Correlations Among P Variables}
\description{This function computes Wilcox H-test and the Steiger-Hakstian-Test for testing H0: R = I. 
}
\usage{WilcoxH(x, use = "pairwise.complete.obs")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Data frame or matrix of dimension n*p with survival times (\code{NA}'s allowed).}
  \item{use}{Treatment of \code{NA}'s for the computation of the correlation matrix (see \code{cor()}). Either \code{"all.obs"}, \code{"complete.obs"}, or \code{"pairwise.complete.obs"}} 
}
\details{This test is robust against violations of normality. Since \code{phmclust()} assumes independence across pages, this test can be used to explore the appropriateness of the data.
}
\value{
  Returns an object of class \code{"wilcoxh"} with the following values:
  \item{Rmat}{Correlation matrix}
  \item{SH.res}{Results for Steiger-Hakstian-Test}
  \item{WH.res}{Results for Wilcox H-test}
}
\references{
Wilcox, R. (1997). Tests of independence and zero correlations among P variables. Biometrical Journal, 2, 183-193.
}

\seealso{\code{\link{phmclust}}}
\examples{

data(webshop)
res <- WilcoxH(webshop)
res

}

\keyword{ models }

