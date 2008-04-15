`msprep` <-
function(linkdata, linkfreq.res, stage.depth = 2, elim = TRUE, sess=TRUE)
{
#linkdata must be of the following structure: "sessionID", "userID", "from", "to", "dwelltime' 
#freqlink...output from linkanalysis
#cycle.depth...split up cyclic surfing behavior into e.g. 2 different stages,
#i.e. surfing 1-2-1-2-1-2-... leads to 1.1 and 1.2 (second time on 1)
#elim...eliminates PI's > cycle depth, if FALSE, successive PI's go into cycle.depth
#sess...if TRUE, session-wise analysis, if FALSE, user-wise analysis

#--------------------------- begin link preparation----------------------
  freqlink <- linkfreq.res$freqlink
  ftfreq <- freqlink[,1:2]
  TFfreq <- freqlink[,1]==freqlink[,2]
  selfcon <- freqlink[freqlink[,1]==freqlink[,2],][,1]        #extract all selfconnected nodes
  othercon <- freqlink[freqlink[,1]!=freqlink[,2],][,1:2]     #extract all transistions to other pages
  #othervec <- as.list(as.vector(othercon))

  #-----------FIXME---------
  #determine selfonly as a vector with only self connected nodes
  selfonly <- NULL
  #------------end FIXME------
  
  if (length(selfonly)==0) {
    rellink1 <- freqlink[,1:2] 
  } else {
    selfonlyL <- as.list(selfonly)
    TFmats <- sapply(selfonlyL, function(x) {TFvecs <- x==freqlink[,1]})
    TFvecs <- apply(TFmats,1,sum)
    rellink1 <- freqlink[TFvecs==0,][,1:2]                      #matrix with relevant links (1)
  }

  tabrel1 <- as.data.frame(table(rellink1[,1],deparse.level=2))
  elimpage <- tabrel1[(tabrel1[,2]==1),]                      #eliminate pages with only 1 outgoing link 
  if (dim(elimpage)[1]==0) {
    rellink2 <- rellink1 
  } else {
    elvecL <- as.list(elimpage[,1])
    TFmate <- sapply(elvecL, function(x) {TFvece <- x==rellink1[,1]})
    TFvece <- apply(TFmate,1,sum)
    rellink2 <- rellink1[TFvece==0,]                             #final matrix with relevant links (2)
  }

#---------------------------------end link preparation-------------------------------------

  
  TFmatl <- apply(rellink2,1, function(x) {
                               TFvec2 <- t(apply(linkdata[,3:4],1,function(y) x==y)) #TF matrix with 2 col
                               TFvecl <- TFvec2[,1] & TFvec2[,2]       #only TRUE&TRUE=TRUE
                               })
  TFvecl <- apply(TFmatl,1,sum)                                          #sum = 1 --> relevant link
  ldnew <-  linkdata[TFvecl==1,]                                         #data with relevant links only

  if (sess==TRUE) {
    ldnewL <- split(ldnew,ldnew[,1])             #session split as list
  } else  {
    ldnewL <- split(ldnew,ldnew[,2])             #user split as list
  }    

  sufinalL <- lapply(ldnewL,function(s_uL) {                       #cycle analysis
                            s_uLs <- s_uL[order(s_uL[,3]),]
                            nvisit <- sequence(table(s_uLs[,3]))
    
                            if (!elim) {                   #eliminate sessions with higher cycle depths
                              s_uLn <- cbind(s_uLs,nvisit)
                              s.u <- s_uLn[nvisit <= stage.depth,]
                            } else {                           #higher cycle depths as max cycle depth
                                nvisit[nvisit > stage.depth] <- stage.depth
                                s.u <- cbind(s_uLs,nvisit)
                            } 

                            s.u[,3] <- paste(s.u[,3],".", s.u[,7], sep = "")
                            #s.u[,3] <- s.u[,3]+s.u[,7]/10     #only the cycle of "from" site is relevant for transition!!!
                            return(s.u)
                      })
  #sumat <- NULL
  #for (i in 1:length(sufinalL)) {sumat <- rbind(sumat,sufinalL[[i]])}
  sudf <- do.call(rbind, sufinalL)

  dummy <- sudf[,6]                                              #put buy in the last column
  sudf <- (sudf[,-6])
  sudf <- cbind(sudf,dummy)
#colnames(sumat) <- c("sess_ID","user_ID","spage","lpage","dwell","nvisit","buy")
return(msdata = sudf)                                        
}

