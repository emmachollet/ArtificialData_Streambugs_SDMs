# update for package


calc.res.steadystate.update <- function(res, times, par, y.names, ratio.tail = 1/20, threshold = 0.001, warn.pom = T){
  
  # res=temp.res
  # times=tout
  # par=temp.parfix
  # y.names=temp.y.names
  # ratio.tail <- 1/20
  # threshold <- 0.001
  # warn.pom <- F
  
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
  
  # t <- c(ceiling(length(times)*19/20):length(times)) # old code
  t <- c(ceiling(length(times)*(1-ratio.tail)):length(times))
  
  # Construct matrix D.mod with results converted from mass per unit length into mass per unit 
  # surface area:
  
  # Could be simplified for our case with constant width and fA (only one habitat)
  D.mod <- matrix(NA,nrow=length(t),ncol=ncol(res)-1)
  colnames(D.mod)<-colnames(res)[-1]
  rownames(D.mod)<-times[t]
  
  for ( j in 1:length(t)){
    # update w and fA:
    
    inpvals <- streambugs:::interpolate.inputs(inp,res[j,1])
    
    w  <- streambugs:::streambugs.update.envcond.reach(par.envcond.w,inpvals)[,"w"]
    fA <- streambugs:::streambugs.update.envcond.habitat(par.envcond.fA,inpvals,y.names$ind.fA)[,"fA"]
    
    # convert results in mass per unit length into mass per unit 
    # surface area:
    
    D.mod[j,] <- res[t[j],-1]/(w*fA)
  }
  
  # calculate the mean of the selected tail of results for each state variable
  D.mean.mod <- colMeans(D.mod)  #gDM/m2
  
  # Check for steady state by calculating the standard deviation of each state variable.
  D.sd.mod <- apply(D.mod,2,sd) # calculate standard deviation of each column (state variable) of D.mod
  if(warn.pom){ # if we include warnings for POM
    c.ind <- which(D.sd.mod>threshold)
  } else { # if we exclude warnings from POM
    c.ind <- which(D.sd.mod>threshold & !grepl("POM", names(D.sd.mod)))
  }
  # if(any(D.sd.mod>0.001)) {warning("steady state probably not reached for :\n", # old code
  #                                 paste(colnames(D.mod)[which(D.sd.mod>0.001)],
  #                                 "\n")) }
  if(any(c.ind)) {
    warning("steady state probably not reached for :\n",
                                       paste(colnames(D.mod)[c.ind],
                                             "\n")) 
    vect.warnings <- colnames(D.mod)[c.ind]
    }
  
  return(list("D.mean.mod" = D.mean.mod, "Warnings" = vect.warnings))
}
