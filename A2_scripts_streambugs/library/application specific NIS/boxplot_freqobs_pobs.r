boxplot.freqobs.pobs <- function(freq.obs,p.1.mat,cex.lab=1.0,mfrow=c(1,1),...)
{
  par(mfrow=mfrow,mgp=c(2.5, 1, 0),
      mar=c(5, 4, 2, 1) + 0.1)  #c(bottom, left, top, right)
  
  y <- p.1.mat
  
  ind.p1 <- list()
  
  for(i in 1:nrow(freq.obs))
  {
    site <- rownames(freq.obs)[i]
    for(j in 1:ncol(freq.obs))
    {
      invert <- colnames(freq.obs)[j]
      ind.p1[[paste(site,invert,sep="_")]] <- c(intersect(grep(paste(site,"_",sep=""),colnames(p.1.mat)),
                                                          grep(paste("_",invert,"_",sep=""),colnames(p.1.mat))))
    }
  }
  
  x <- NULL
  
  for(i in 1:nrow(freq.obs))
  {
    x.i <- freq.obs[i,]
    names(x.i) <- paste(rownames(freq.obs)[i],colnames(freq.obs),sep="_")
    
    x <- c(x,x.i)
  }
  
  ind.0    <- which(x<0.25)
  ind.0.25 <- which(x>=0.25 & x<0.5)
  ind.0.5  <- which(x>=0.5 & x<0.75)
  ind.0.75 <- which(x>=0.75)
    
  boxplot(list("0-0.25"  = unlist(y[,unlist(ind.p1[names(ind.0)])]),
               "0.25-0.5"= unlist(y[,unlist(ind.p1[names(ind.0.25)])]),
               "0.5-0.75"= unlist(y[,unlist(ind.p1[names(ind.0.5)])]),
               "0.75-1"  = unlist(y[,unlist(ind.p1[names(ind.0.75)])])),
          cex.lab=cex.lab,
          ylim=c(0,1),
          ylab=expression(paste("p.obs=p(y=1)",sep="")),
          xlab="frequency of observation",
          #varwidth=TRUE,
          ...)  
  
#   ind.ord <- order(x)
#   names(ind.ord) <- names(x)[ind.ord]
#   f.obs <- x[ind.ord]
#   
#   list.bx <- list()
#   list.mean <- list()
#   
#   for(i in 1:length(ind.ord))
#   {
#     list.bx[[paste(names(ind.ord)[i])]] <- unlist( y[,unlist(ind.p1[names(ind.ord)[i]])]) 
#     list.mean[[paste(names(ind.ord)[i])]] <- mean(unlist( y[,unlist(ind.p1[names(ind.ord)[i]])]) )
#   }
#   
#   mean.p1 <- unlist(list.mean)
#   
#   plot(f.obs,mean.p1)
#   
#   data <- cbind(f.obs,mean.p1)
#   
#   boxplot(mean.p1 ~ f.obs,data=data,xlab="f.obs",ylab="mean(p(y=1))",varwidth=TRUE)
#   
#   boxplot(list.bx)
 
}

#plot(x,y)