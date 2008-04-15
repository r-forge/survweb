`plot.mscif` <-
function(x, from, main = paste("CIF from ",from), curvlab, ylim = c(0, 1), xlim, wh = 2, xlab= "Seconds",
         ylab = "Probability", lty, color = 1, lwd = 1, ...)
{
#msres...result from the previous ms-computation
#from...starting page
  msres1 <- x$msres
  msdata1 <- x$msdata
  fstatvec <- sort(unique(msdata1[(msdata1[,3]==from),4]))

  spagevec <- NULL
  for (i in 1:length(msres1)) {
    spagevec <- c(spagevec,msres1[[i]][[1]][1])
  }
  pp <- which(spagevec == from)                                 #plotting position
  plotres <- msres1[[pp]]

  if (missing(lty)) lty <- 1:length(fstatvec)
  if (missing(curvlab)) curvlab <- names(plotres[[2]])[-(length(names(plotres[[2]])))]
  if (missing(xlim)) xlim <- c(0,100)
   
  plot.cuminc(plotres[[2]], xlab = xlab, ylab = ylab, main = main, color = color,
              curvlab = curvlab, wh = wh, lwd = lwd, ylim = ylim, xlim = xlim,...)

}

