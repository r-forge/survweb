`mscuminc` <-
function(msdata, cutpoint = NULL, censoring = TRUE)                                                
{
#msdata...output from the previous data (list)
#clean... eliminate PI's with dwell times larger

msdata1 <- msdata
msdata1[,3] <- as.factor(msdata1[,3])

if (!is.null(cutpoint)) {
  if (censoring) {
    censpos <- which(msdata1[,5] > cutpoint)
    if (!is.numeric(msdata1[,4])) {
      ccode <- "cens"
      levels(msdata1[,4]) <- c(levels(msdata1[,4]),ccode)
      msdata1[censpos,4] <- ccode
    } else {                           #0 as censcode
      ccode <- 0
      msdata1[censpos,4] <- ccode
    }      
  } else {
    msdata1 <- msdata1[(msdata1[,5] < cutpoint),]            #keep only dwell times < clean
  }
}

fromvec <- unique(msdata[,3])           #starting pages
ind <- 1:(length(fromvec))
msL <- tapply(fromvec,ind,function(from) {
                            fromdata <- msdata1[msdata1[,3]==from,]
                            if (length(fromdata) > 0) {
                              fstatus <- fromdata[,4]              #failure types
                              ftimes <- fromdata[,5]
                              outms <- cuminc(ftimes,fstatus,group=fromdata[,7], cencode = ccode) #test for the cumulative incidence function
                              frommat <- as.matrix(from)
                              rownames(frommat) <- "starting page"
                              list(frommat,outms)
                            }
                            })
res <- list(msdata = msdata1, msres = msL)                                       #output list (starting page, estimates)
class(res) <- "mscif"
res
}

