`plot.survival` <-
function(x, gr.subset, var.subset, group = TRUE, xlim = NA, ylim = NA, xlab = "Survival Time", ylab = "Survival Function", 
         main = "Survival Functions", type = "l", lty = 1, lwd = 1, col = NA, legpos = "right", ...)


{

#x...object of class "mws"

scale <- 1/x$scale						       #different parameterization
K     <- dim(scale)[1]
pages <- dim(scale)[2]

if (missing(gr.subset)) gr.subset <- 1:K                 #no subset defined for group-wise plot (group title)
if (missing(var.subset)) var.subset <- 1:pages            #no subset defined for page-wise plot (page title)

if ((any(is.na(col))) && (group)) col <- var.subset
if ((any(is.na(col))) && (!group)) col <- gr.subset
if (any(is.na(xlim))) xlim <- c(0,mean(x$clmean))

maxtime <- round(xlim[2],0)
Ti <- maxtime                        #upper limit for sojourn time on x-axis
survgroup <- array(NA,c(Ti,pages,K))  #init plot array

if (group) {                           #1 plot for each group
  for (k in 1:K) {                     #Compute hazards for each page in group k
    for (j in 1:pages) {
      survgroup[,j,k] <- exp(-((scale[k,j])*(1:Ti))^(x$shape[k,j]))		
   }}
} else {                               #1 plot for each page
for (j in 1:pages) {
  for (k in 1:K) {
     survgroup[,j,k] <- exp(-((scale[k,j])*(1:Ti))^(x$shape[k,j]))	
   }}
}
dimnames(survgroup) <- list(NULL,colnames(x$shape),rownames(x$shape))

#----------------------------produce plots-------------------------- 
if (any(is.na(ylim))) ylim <- c(0,max(survgroup))

K.lab <- dimnames(survgroup[,,gr.subset])[3][[1]]            #labels
p.lab <- dimnames(survgroup[,var.subset,])[2][[1]]
K.ss <- length(gr.subset)                                   #maximum index
p.ss <- length(var.subset)

survvec <- as.vector(survgroup[,var.subset,gr.subset])                       #survival array as vector
page.g <- factor(rep(rep(1:p.ss, each = Ti),K.ss),labels = p.lab)
group.g <- factor(rep(1:K.ss, each = Ti*p.ss),labels = K.lab)
x1 <- rep(1:Ti,p.ss*K.ss)


if (group) {                                    #Group-wise plot
  xyplot(survvec~x1|group.g, groups=page.g, type=type,xlab=xlab,ylab=ylab,main=main,
      lty=lty,lwd=lwd,ylim=ylim,col=col,
      key=list(text=list(lab=p.lab,col=1),lines=TRUE,lty=lty,col=col,space=legpos),...)      
} else {
  xyplot(survvec~x1|page.g, groups=group.g, type=type,xlab=xlab,ylab=ylab,main=main,
      lty=lty,lwd=lwd,ylim=ylim,col=col,
      key=list(text=list(lab=K.lab,col=1),lines=TRUE,lty=lty,col=col,space=legpos),...)  
}    
}

