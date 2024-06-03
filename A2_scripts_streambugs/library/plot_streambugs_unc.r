# 27.2.2023 adapted to run based on streambugs package, see lines 23,30,45,47,48

plot.streambugs.unc <- function(res.det,
                                res.quant.l=NA,
                                res.quant.u=NA,
                                res.median=NA,
                                par,
                                inp=NA)
{
  
  # res.det <- res.array[,,1]
  # res.quant.u <- res.quant.95
  # res.quant.l <- res.quant.05
  # res.median  <- res.quant.50
  # par<-par.samp[,1]
  # inp=NA
  
  res <- res.det

  sys.def <- streambugs.get.sys.def(y.names=colnames(res)[-1],par=par,inp=inp)
  y.names <- sys.def$y.names
  
  par.envcond.w <- streambugs:::get.inpind.parval.envcond.reach(
    par.names = c("w"),
    y.names   = y.names,
    par       = par,
    inp       = inp,
    required  = c("w"))
  
  par.envcond.fA <- streambugs:::get.inpind.parval.envcond.habitat(
    par.names = c("fA"),
    y.names   = y.names,
    par       = par,
    inp       = inp,
    defaults  = c(fA=1))

  # evaluate time-dependent inputs and merge them with parameters
  # (overwriting parameters if one exists with the same name):
  # -------------------------------------------------------------
  
  for ( j in 1:nrow(res) )
  {
    # update (time-dependent) parameters:
    
    inpvals <- streambugs:::interpolate.inputs(inp,res[j,1])
    
    w  <- streambugs:::streambugs.update.envcond.reach(par.envcond.w,inpvals)[,"w"]
    fA <- streambugs:::streambugs.update.envcond.habitat(par.envcond.fA,inpvals,y.names$ind.fA)[,"fA"]
    
    # convert results in mass per unit length into mass per unit 
    # surface area:
    
    res[j,-1] <- res[j,-1]/(w*fA)

    if(sum(!is.na(res.quant.l))>0) { res.quant.l[j,-1] <- res.quant.l[j,-1]/(w*fA) }
    if(sum(!is.na(res.quant.u))>0) { res.quant.u[j,-1] <- res.quant.u[j,-1]/(w*fA) }
    if(sum(!is.na(res.median)) >0) { res.median[j,-1]  <- res.median[j,-1]/(w*fA)  }
  }
    
  # determine maxima of converted results:
  
  y.max <- rep(NA,length(y.names$y.names))
  for ( i in 1:length(y.names$y.names) )
  {
    if(sum(is.na(res.quant.u))==0)
    {
      y.max[i] <- 1.1*max(res.quant.u[,1+i])
    } else
    {
      y.max[i] <- 1.1*max(res[,1+i])
    }
    
  }
  names(y.max) <- y.names$y.names
  
  # get time:
  
  t    <- res[,1]
  
  # plot output by y:
  
  par.def <- par(no.readonly=TRUE)
  #par(mfrow=c(length(y.names$habitats),length(y.names$groups)))
  par(mfrow=c(5,3))   #c(bottom, left, top, right) ,oma=c(2,1,2,1)
  
  for ( i in 1:length(y.names$y.names) )
  {
    plot(numeric(0),numeric(0),type="n",
         xlim=c(min(t),max(t)),ylim=c(0,y.max[i]),
         xlab="time [yr]",ylab="biomass [gDM/m2]",
         main=paste(y.names$y.reaches[i],y.names$y.habitats[i],y.names$y.taxa[i]))
    
    if(sum(!is.na(res.quant.u))>0)
    {
      polygon(c(res[1,1],res[,1],res[,1][length(res[,1])]),
              c(0,res.quant.u[,1+i],0),
              col=grey(0.75),lty="blank")
    }
    
    if(sum(!is.na(res.quant.l))>0)
    {     
      polygon(c(res[1,1],res[,1],res[,1][length(res[,1])]),
              c(0,res.quant.l[,1+i],0),
              col="white",lty="blank")
    }
    
    lines(t,res[,1+i],lty=1)
    
    if(sum(!is.na(res.median))>0)
    {
      lines(t,res.median[,1+i],lty=2,col=1)
    }
  }
  
  par(par.def)
}