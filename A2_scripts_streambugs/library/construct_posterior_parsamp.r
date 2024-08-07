# Bayesian inference:
# ===================
# 
# generate posterior parsamp by
# saving parameter combinations which lead to a foodweb that is in 
# accordance with the observed foodwebs
# -----------------------------------------------------------------

construct.posterior.parsamp <- function(y.names,res.array,par.samp.nona,observations,y.exclude=NA)
  
{
  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)
  
  nres <- dim(res.array)[3] 
  npar <- dim(par.samp.nona)[2]
  if(nres != npar) warning("number of res and par samples not equal, nres: ",nres,", npar: ",npar)
  
  # always observed taxa:
  ind.always <- which(observations=="always") 
  ind.tax <-  grep(paste(names(ind.always),collapse="|"),y.names$y.names)
  always.obs  <- y.names$y.names[ind.tax]  
  
  # never observed taxa:
  ind.never <- which(observations=="never") 
  ind.tax <-  grep(paste(names(ind.never),collapse="|"),y.names$y.names)
  never.obs  <- y.names$y.names[ind.tax]  
  
  if (length(y.exclude)>1 & sum(is.na(y.exclude))==0)
  {
    ind.excl <- grep(paste(y.exclude$taxa.underest,collapse="|"),always.obs)
    always.obs <- always.obs[-ind.excl]
    
    ind.excl2 <- grep(paste(y.exclude$taxa.overest,collapse="|"),never.obs)
    never.obs <- never.obs[-ind.excl2]    
    
    cat("no of taxa excluded: ",length(c(ind.excl,ind.excl2)),"\n")
  }
  
  ind.match <- numeric(0)
  
  for (i in 1:nres)
  {
    survivals <- check.survival(res=res.array[,,i],times=tout,
                                par=par.samp.nona[,i],Dens.crit=D.crit)                                                              
    
      criterion1 <- sum(survivals[always.obs] == "survived") == length(always.obs)
      criterion2 <- sum(survivals[never.obs]  == "extinct")  == length(never.obs)
    
    if (criterion1 & criterion2) ind.match <- c(ind.match,n)       
  }
  
  par.samp.post <- par.samp.nona[ind.match]
  
  cat("no of matches:",length(ind.match),"of",nres)
  
  return(par.samp.post)
  
}
