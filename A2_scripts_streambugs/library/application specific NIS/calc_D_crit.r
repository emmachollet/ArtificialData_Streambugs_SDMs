
calc.D.crit <- function(y.names,par,par.prior.delta=NULL)
{
  if(!is.list(y.names)) y.names <- decode.statevarnames(y.names)
  
  if(length(par.prior.delta)>0) par.fix <- c(par,par.prior.delta)  
  
  par.fix <- back.transform.parameters(par.fix)
  
  M.taxa <- streambugs.get.sys.def(y.names=y.names,par=par.fix)$par.taxaprop.direct$parvals[,"M"]
  
  Abund.crit <- 0.5  #iWaQa 1/(3*3*0.048)m2 AWEL 0.5# No of individuals/m2
  D.crit <- NA       # vector Dcrit der inverts in gDM/m2
  
  for (i in 1:length(y.names$y.names))
  {
    y.name <-y.names$y.names[i]
    ind <- grep(paste(y.name),names(M.taxa)) #rownames(M.taxa))
    M.i <- unique(M.taxa[ind])
    if( length(M.i)==1 )
    { 
      D.crit[i] <- M.i*Abund.crit
      names(D.crit)[i] <- y.name
    } else 
    { warning("Mass of ",y.name," is not unique: ",signif(M.i,2)) }
  }
  
  return(D.crit)
}