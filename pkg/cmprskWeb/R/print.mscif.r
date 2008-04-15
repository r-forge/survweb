print.mscif <- function(x, ...)
{
 #print method for mscif
 
 nres <- length(x$msres)
 testlist <- rep(as.list(NA),nres)
 listnames <- NULL
 for (i in 1:nres) {
   if (is.null(x$msres[[i]][[2]]$Tests)) x$msres[[i]][[2]]$Tests <- NA
   list.el <- round(x$msres[[i]][[2]]$Tests,4)
   if(!is.na(list.el[1])) colnames(list.el) <- c("X^2", "p-value", "df")
   testlist[[i]] <- list.el
   listnames <- c(listnames,x$msres[[i]][[1]][,1])
 }
 names(testlist) <- listnames
 listnames <- paste("Transition from", listnames)
 
 cat("CIF tests: \n\n")
 for (i in 1:nres) {
   cat("\n")
   print(listnames[i])
   print(testlist[[i]])
 }
 
}