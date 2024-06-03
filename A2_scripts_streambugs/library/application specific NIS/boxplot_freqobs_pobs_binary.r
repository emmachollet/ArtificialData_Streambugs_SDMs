boxplot.freqobs.pobs.binary <- function(freq.obs,p.1.mat,cex.lab=1.0,mfrow=c(1,1),...)
{
  par(mfrow=mfrow,mgp=c(2.5, 1, 0),
      mar=c(5, 4, 2, 1) + 0.1)  #c(bottom, left, top, right)
  
  y <- p.1.mat
  
  ind.p1 <- list()
  
  for(i in 1:nrow(freq.obs))
  {
    site_hab <- rownames(freq.obs)[i]
    for(j in 1:ncol(freq.obs))
    {
      invert <- colnames(freq.obs)[j]
      ind.p1[[paste(site_hab,invert,sep="_")]] <- c(intersect(grep(paste(site_hab,"_",sep=""),colnames(p.1.mat)),
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
  
  ind.0    <- which(x<0.5)
  ind.1    <- which(x>=0.5)
  
  boxplot(list("0-<0.5" = unlist(y[,unlist(ind.p1[names(ind.0)])]),
               "0.5-1"  = unlist(y[,unlist(ind.p1[names(ind.1)])])),
          cex.lab=cex.lab,
          ylim=c(0,1),
          ylab=expression(paste("p.obs=p(y=1)",sep="")),
          xlab="frequency of observation",
          #varwidth=TRUE,
          ...) 
  
  plot(NA, ylim=c(0,1), xlim=c(0,1), ylab=expression(paste("p.obs=p(y=1)",sep="")),xlab="frequency of observation",
       cex.lab=cex.lab,...)#  xaxt="n",
  
  vioplot(x = unlist(y[,unlist(ind.p1[names(ind.0)])]),
          x2 = unlist(y[,unlist(ind.p1[names(ind.1)])]),
          names=c("f.obs:0-<0.5","f.obs:0.5-1"),
          col=gray(0.7),at=c(0.25,0.75),
          drawRect=FALSE,add=TRUE,wex=0.62)
  
  
  #   hist(unlist(y[,unlist(ind.p1[names(ind.0)])]),breaks=100)
  #   hist(unlist(y[,unlist(ind.p1[names(ind.1)])]),breaks=100)
  #   
}