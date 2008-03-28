\name{plot.hazard}
\alias{plot.hazard}
\alias{plot.profile}
\alias{plot.survival}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot functions}
\description{
Plotting functions for hazard rates, survival times and cluster profiles. 
}
\usage{
plot.hazard (x, gr.subset, var.subset, group = TRUE, xlim = NA, ylim = NA, xlab = "Survival Time", ylab = "Hazard Function", 
main = "Hazard Functions", type = "l", lty = 1, lwd = 1, col = NA, legpos = "right", ...)

plot.survival (x, gr.subset, var.subset, group = TRUE, xlim = NA, ylim = NA, xlab = "Survival Time", ylab = "Survival Function", 
main = "Survival Functions", type = "l", lty = 1, lwd = 1, col = NA, legpos = "right", ...)

plot.profile (x, method = "mean", type = "b", pch = 19, lty = 1, lwd = 1, col = NA, 
xlab = "Variables", leglab = NA, ylab = NA, main = NA, legpos = "topright", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class \code{mws} from \code{phmclust}}
  \item{gr.subset}{Optional vector for plotting subset of clusters}
  \item{var.subset}{Optional vector for plotting subset of variables}
  \item{group}{if \code{TRUE} hazard/survival plots are produced for each group, if \code{FALSe} for each variable}
  \item{method}{\code{"mean"} for cluster mean profile plot and \code{"median"} for cluster median profile plot}
  \item{xlim}{limits for x-axis}
  \item{ylim}{limits for y-axis}
  \item{xlab}{label for x-axis}
  \item{ylab}{label for y-axis}
  \item{main}{title of the plot}
  \item{leglab}{label for the legend}
  \item{type}{type of plot}
  \item{lty}{line type}
  \item{lwd}{line width}
  \item{pch}{type of plotting points}
  \item{col}{colors; if \code{NA} it is determined in the function}
  \item{legpos}{position of the legend; \code{"topright"},\code{"topleft"},\code{"bottomright"},
  \code{"bottomleft"},\code{"left"},\code{"right"},\code{"top"}, or \code{"center"}}
  \item{\dots}{Additional plot options}
}


\seealso{\code{\link{phmclust}}}
\examples{

##Plots for mixture Weibull model with 3 components
data(webshop)
res <- phmclust(webshop, 3)

##Hazard plot for first and third group, all pages
plot.hazard(res, gr.subset = c(1,3), group = TRUE, xlab = "Dwell-Time")

##Survival plot for each group, first 6 pages
plot.survival(res, var.subset= 1:6, group = FALSE, xlab = "Dwell-Time")

##Cluster profile plot
plot.profile(res, xlab = "Pages", ylab = "Mean Dwell-Time", main = "Cluster Profile")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ hplot }
