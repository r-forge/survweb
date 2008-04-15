print.linkfreq <- function(x, ...)
{
#x ... object of class linkfreq
  cat("Transition frequencies: \n")
  print(x$linktable)
  cat("\n")
}