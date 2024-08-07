generate.par.samp.vector <- function(prior.par.defs)
{
  
  # prior.par.defs
  
  par <- rep(NA,length(prior.par.defs))
  names(par) <- names(prior.par.defs)
  for ( p in 1:length(prior.par.defs))
  {
    
    if(length(prior.par.defs[[p]]) > 1)
    {
      par[p] <- unlist(sysanal.randsamp(1,dist="Indep",distdef=prior.par.defs[p])$sample)
      
    } else {
      
      par[p]  <- prior.par.defs[[p]]
    } 
    
  }
  
# # convert strings to numerics 
#   par <- encode.envpars(as.matrix(par))  #(not longer needed)
#  class(par) <- "numeric"
  
# # convert T in degC to K and CSusPOM to DSusPOM by multiplying with Filt_scope:
#  par.m <- convert.TdegC.CSusPOM(as.matrix(par))   #(not longer needed here)
#  par <- as.vector(par)
#  names(par) <- rownames(par.m)

  return(par)
}

