calc.prior.dens <- function(par,par.prior.delta,prior.par.defs,log=FALSE)
    
  # calc prior dens
  
  # arguments:
  # par: named vector with parameter values of all inference parameters
  # prior.par.defs (globally defined)
  # log: if TRUE the log of density is returned
  # returns:
  # the joint prior density constructed from independent marginals 
  # (for discrete distributions the probabilities and not densities are taken)
  
{
  
#   source("libraries/convert_CSusPOM.r")
#   source("libraries/sysanal.r")
  
  if(sum(grepl("_DSusPOM",names(prior.par.defs)))>0) 
  {
    par.conv <- revert.CSusPOM(par=c(par,par.prior.delta)) 
    par <- par.conv(names(par))
  }
  
  d <- NULL
  
  for(i in 1:length(prior.par.defs))
  {
       
    indpar <- which(names(par)==names(prior.par.defs[i]))
    
    if(prior.par.defs[[i]][1]=="Discrete" | prior.par.defs[[i]][1]=="discrete")
    {
      #cat(i,names(prior.par.defs[i]))
      nclasses <- (length(prior.par.defs[[i]])-1)/2
      classes <- as.matrix(prior.par.defs[[i]][-1][1:nclasses])
      probs <- prior.par.defs[[i]][(nclasses+2):length(prior.par.defs[[i]])]
      
      rownames(classes) <- paste(rep(names(prior.par.defs[i]),times=nclasses))
      
      classes.numeric <- encode.envpars(classes)
      
      indpar <- which(names(par)==unique(rownames(classes.numeric)))
      
      indclass <- which(classes.numeric==par[indpar])
      
      d.i <- as.numeric(probs[indclass])
      if(log==TRUE) d.i <- log(d.i)
      
    } else { 
      
      d.i <- sysanal.calcpdf(x=par[indpar],distpar=prior.par.defs[[i]],log=log)
      
    }
    
    d <- c(d,as.numeric(d.i))
  }
  
  if(log==FALSE)
  {
    d.mult <- prod(d,na.rm=TRUE)
  } else {d.mult <- sum(d,na.rm=FALSE)}
  
  
  return(d.mult)
}