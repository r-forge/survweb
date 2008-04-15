`linkfreq` <-
function(from, to, linkcut = 0.001)
{
#from ... vector of starting pages
#to ... vector of landing pages (transitions)  
#performs link analysis for page super types                         
#linkcut...cutting value for determining frequent transitions,
#e.g. 0.005 means that links used in more than 0.5% of the transitions

#returns list containing:
# $freqlink: transition frequencies (driected; as matrix)
# $linktable: table of transition frequencies (directed)
# $fromSum: outgoing transitions (from)
# $toSum: incoming transitions (to)

  fsup <- from                                     #from
  tsup <- to                                      #to
  linkcut1 <- length(from)*linkcut
  suptable <- table(from,to)
  pvek <- rownames(suptable)
  fromSum <- rowSums(suptable)                              #marginals
  toSum <- colSums(suptable)

#---determining frequent transitions (above linkcut)
  N <- dim(suptable)[1]
  ind1 <- permutations(N,2,v=pvek,repeats.allowed=TRUE)     #linked pages
  TFmat <- suptable > linkcut1
  TFvec <- as.vector(t(TFmat))
  freqtrans <- ind1[TFvec,]
  freq1 <- t(suptable)[t(TFmat)]
  freqlink <- as.data.frame(cbind(freqtrans,freq1))
  colnames(freqlink) <- c("from","to","freq")

  res <- list(freqlink = freqlink,linktable = suptable, fromSum = fromSum, toSum = toSum)
  class(res) <- "linkfreq"
  res  
}

