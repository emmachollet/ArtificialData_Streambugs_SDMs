# update for package


calc.res.steadystate <- function(res,times,par,y.names)                    #                                                               
  
  # for testing: times=tout;par=par.fix
  
{
  
  if(length(times)>nrow(res)) stop("length(times) larger than nrow(res), simulations not successful")
  
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
  
  t <- c(ceiling(length(times)*19/20):length(times))
 
  # Construct matrix D.mod with results converted from mass per unit length into mass per unit 
  # surface area:
  
  D.mod <- matrix(NA,nrow=length(t),ncol=ncol(res)-1)
  colnames(D.mod)<-colnames(res)[-1]
  rownames(D.mod)<-times[t]
  
  for ( j in 1:length(t))
  {
    # update w and fA:
    
    inpvals <- streambugs:::interpolate.inputs(inp,res[j,1])
    
    w  <- streambugs:::streambugs.update.envcond.reach(par.envcond.w,inpvals)[,"w"]
    fA <- streambugs:::streambugs.update.envcond.habitat(par.envcond.fA,inpvals,y.names$ind.fA)[,"fA"]
    
    # convert results in mass per unit length into mass per unit 
    # surface area:
    
    D.mod[j,] <- res[t[j],-1]/(w*fA)
  }
  
  D.mean.mod <- colMeans(D.mod)  #gDM/m2
  
  # check steady state:
  D.sd.mod <- apply(D.mod,2,sd)
  if(any(D.sd.mod>0.001)) {warning("steady state probably not reached for :\n",
                                  paste(colnames(D.mod)[which(D.sd.mod>0.001)],
                                  "\n")) }
  
  return(D.mean.mod)
}