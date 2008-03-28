`phmclust` <-
function(x, K, method="separate", Sdist= "weibull", EMstart=NA, EMoption="classification", 
EMstop=0.0001, maxiter=1000)
{

# Input:
#        x ........ n x p data-matrix each column representing dwell-time of an individual
#                   session in site-area 1 to p. pages not visited should have a dwell time of 0 or NA
#        K ........ scalar with number of mixture components
#        method...  denotes the model to be fitted: possible values are "separate", "main.g", "main.p", "int.gp","main.gp", 
#                   while in method "separate" distributions of individual groups are 
#                   estimated independently, method "main.g" assumes that there is a common
#                   base-line hazard which is common to all groups, "main.p" the same for the sites. method "main.gp" 
#                   fits a main effects model whereas in model "int.gp" allows interaction effects.
#        EMstart .. Vector of length n with starting values for group membership
#        EMoption.. "classification" is based on the probabilistic cluster assignment, "maximization" on deterministic assignment,
#                   "randomization" provides a posterior-based randomized cluster assignement. 
#        Sdist .... Survival distribution for WPHM model. These include "weibull", "exponential", "rayleigh". 
#        
# Output:
#        list containig the following elements
#
#
# Packages Needed: MASS, Survival

if (is.data.frame(x)) x <- as.matrix(x)
if (is.vector(x)) stop("x must be a data frame or a matrix with more than 1 columns!")

pvisit.est <- (any(is.na(x)) | any(x==0))                          # TRUE if visiting prob estimated
n <- nrow(x)                                                       # n ... number of sessions
p <- ncol(x)                                                       # p ... number of site-areas                            

d0 <- s.check(x=x,K=K,n=n,EMstart=EMstart,EMoption=EMoption,method=method,Sdist=Sdist)    #sanity checks

x <- d0$x
EMstart <- d0$EMstart[1:(dim(x)[1])]
method <- d0$method
                                                                                                    
likelihood <- numeric(maxiter)
iter <- 0                                
ConvergEM <- FALSE  


#=============================EM-estimation================================

if (EMoption == "classification") {                                  #maximization EM
    while (ConvergEM == FALSE)
    {
       iter <- iter + 1
       #print(iter)
       d1 <- Eclass(x, EMstart, K=K, method=method, Sdist=Sdist,p)	     #E-Step maximization
       d2 <- Mclass(x, d1$shape,d1$scale,d1$prior,K=K)                    #M-Step maximization
       
       likelihood[iter+1] <- d2$lik.tot                            #likeihood in the current iteration  
       
       if ((iter >= maxiter) || (abs(likelihood[iter+1]-likelihood[iter]) < EMstop)) {   
          ConvergEM <- TRUE 
       } else {
          EMstart <- d2$newgr
          }
    }
postmat <- NULL
newgr <- d2$newgr
}

#========================================================================

if (EMoption == "maximization") {                                #classification EM
    while (ConvergEM == FALSE)
    {
       iter <- iter + 1
       d1 <- Emax(x, EMstart, K=K, method=method, Sdist=Sdist,p)                          #E-Step maximization
       d2 <- Mmax(x, d1$shape,d1$scale,d1$prior,K=K)                  #M-Step maximization
       
       likelihood[iter+1] <- d2$lik.tot                            #likeihood in the current iteration  
       
       if ((iter >= maxiter) || (abs(likelihood[iter+1]-likelihood[iter]) < EMstop)) {   
          ConvergEM <- TRUE 
       } else {
          EMstart <- d2$postmat
       }
    }
postmat <- d2$postmat
newgr <- apply(postmat, 1, function(y){ind <-(1:K)[y==max(y)]})    #final group assignement                         
}

#========================================================================

if (EMoption == "randomization") {                                  #maximization EM
    while (ConvergEM == FALSE)
    {
       iter <- iter + 1
       d1 <- Eclass(x, EMstart, K=K, method=method, Sdist=Sdist,p)	                   #E-Step randomization (=classification)
       d2 <- Mrandom(x, d1$shape,d1$scale,d1$prior,K=K)                    #M-Step maximization
       
       likelihood[iter+1] <- d2$lik.tot                            #likeihood in the current iteration  
       
       if ((iter >= maxiter) || (abs(likelihood[iter+1]-likelihood[iter]) < EMstop)) {   
          ConvergEM <- TRUE 
       } else {
          EMstart <- d2$newgr
       }
    }
postmat <- d2$postmat
newgr <- d2$newgr
}
#========================================================================

if ((iter == maxiter) && (maxiter == 1000)) warning("EM did not converge! Maximum iteration limit reached!")

if (pvisit.est) { 
   anzpar <- d1$anzpar + K*p       #if NA's in x --> number of estimated visiting probabilities added
} else {
   anzpar <- d1$anzpar             #no NA's in x
}


likconv <- likelihood[2:(iter+1)]
aic <- -2*(likelihood[iter+1]-anzpar)
bic <- (-2*likelihood[iter+1])+anzpar*log(n)

clmean.l <- by(x,newgr,function(y) {apply(y,2,mean)})
clmean <- matrix(unlist(clmean.l),nrow=K,byrow=TRUE)

clmed.l <- by(x,newgr,function(y) {apply(y,2,median)})
clmed <- matrix(unlist(clmed.l),nrow=K,byrow=TRUE)

if (is.null(colnames(x))) {
  colnames(d1$scale) <- paste("V",1:dim(d1$scale)[2],sep="")
  colnames(d1$shape) <- paste("V",1:dim(d1$shape)[2],sep="")
} else {
  colnames(d1$scale) <- colnames(d1$shape) <- colnames(x)
}

rownames(d1$shape) <- rownames(d1$scale) <- paste("Cluster",1:K,sep="") 

result <- list(K=K, iter=iter, method=method, Sdist=Sdist, likelihood=likconv,
     pvisit=d1$prior, shape=d1$shape, scale=d1$scale, group=newgr,posteriors=postmat,npar=anzpar,aic=aic,bic=bic,clmean=clmean,clmed=clmed)
class(result) <- "mws"
result
}

