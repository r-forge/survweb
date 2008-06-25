`Mclass` <-
function(x,shape,scale,prior,K)
{
lik <- NULL
for (k in 1:K) {
  shape.vec <- shape[k,]
  scale.vec <- scale[k,]
  prior.vec <- prior[k,]

  likmat <- t(apply(x,1,function(y){
                    la.vec <- prior.vec*(mapply(dweibull,y,shape.vec,scale.vec))        #pages visited by session  (for weibull, exponential, and rayleigh
                    
                    la.vec[is.na(la.vec)] <- 0
                    la.vec[la.vec==0|la.vec==Inf] <- 1-prior.vec[la.vec==0|la.vec==Inf] #dweibull either 0 or Inf for 0 dwell time (pages not visited by session)
                    return(la.vec)
                    }))
  if ((dim(likmat)[1]) == 1) likmat <- t(likmat) #sanity check if x is vector (only 1 page)
  likmat <- log(likmat)
  lik <- cbind(lik,apply(likmat,1,sum))        #multiplying prob over pages
}

lik.n <- apply(lik,1,max)                      #maximum log-likelihood value for each session               
lik.tot <- sum(lik.n)

newgr <- apply(lik, 1, function(y){              #M-step by assigning session to group with max likelihood value
                         h<-(1:K)[y==max(y)]
                         if (length(h)==1) h else sample(h)[1]})
list(newgr=newgr,lik.tot=lik.tot)
}

