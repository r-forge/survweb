#data preparation function for competing risks regression (pos and nvisit covariate)

`crrprep` <- function(linkdata, from, cutpoint = NULL, censoring = TRUE)
{
  #prepares the data set for competing risks by adding the covariates pos and nvisit from a particular spage
  #spage ... starting page from which the outgoing link will be analyzed
  #cutpoint ... specify dwell time cutpoint: either a numeric value (seconds), or NULL for censoring

  spage <- from
  linkdata1 <- linkdata
  if (!is.null(cutpoint)) {
  if (censoring) {
    censpos <- which(linkdata1[,5] > cutpoint)
    if (!is.numeric(linkdata1[,4])) {
      ccode <- "cens"
      levels(linkdata1[,4]) <- c(levels(linkdata1[,4]),ccode)
      linkdata1[censpos,4] <- ccode
    } else {                           #0 as censcode
      ccode <- 0
      linkdata1[censpos,4] <- ccode
    }
  } else {
    linkdata1 <- linkdata1[(linkdata1[,5] < cutpoint),]            #keep only dwell times < clean
  }
  }
  
  #cslist <- split(linkdata1[,1],linkdata1[,1])                            #S_ID as list
  #pos <- as.vector(unlist(lapply(cslist,function(x) (1:length(x)))))    #position in clickstream covariate
  #sessfreq <- sequence(table(xx))[order(xx)]


  indpos <- rank(unique(linkdata1[,1]))               #generate clickstream position covariate
  sessfreq1 <- table(linkdata1[,1])
  sessfreq <- table(linkdata1[,1])[indpos]
  pos <- sequence(sessfreq)
  linkdata1 <- cbind(linkdata1,pos)

  ldfilt <- linkdata1[linkdata1[,3]==spage,]            #filter transistions from spage

  indpos <- rank(unique(ldfilt[,1]))               #generate nvisit covariate
  sessfreq1 <- table(ldfilt[,1])
  sessfreq <- table(ldfilt[,1])[indpos]
  nvisit <- sequence(sessfreq)
  ldfilt1 <- cbind(ldfilt,nvisit)
  
  ldfilt1[,6] <- as.numeric(ldfilt[,6])
  print(table(ldfilt[,4]))
  ldfilt1[,4] <- as.numeric(ldfilt[,4])
  print(table(ldfilt[,4]))

  return(crdata=ldfilt1)
}

#crdata can be used to compute competing risks models such as in "crrweb_xxx_"

#competing risks analysis for the crdata-data generated in "datprepcr"
#covariates are buyer, clickstream position and nvisit.
#function call:
#res <- comprisk(crdata,fcode=7)
