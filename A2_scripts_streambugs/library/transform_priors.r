
# transform logarithmic priors to normal priors

transform.priors<-function(prior.par.defs)
{
  for(i in 1:length(prior.par.defs))
  {
    parname <- names(prior.par.defs[i])
    
    # if (parname=="Algae_M"  | parname=="Algae_Ea" | grepl("_CSusPOM",parname) | 
    #       grepl("_CSusPOM",parname) | grepl("_notfit",parname) | 
    #       grepl("_nonfit",parname) )
    # {
    prior.i <- prior.par.defs[i]
    
    if (prior.i[[1]][1]=="Lognormal" |
          prior.i[[1]][1]=="lognormal" ) 
    {
      prior.par.defs[[i]][1] <- "Normal"
      
      mu    <- as.numeric(prior.i[[1]][2])
      sigma <- as.numeric(prior.i[[1]][3])
      
      prior.par.defs[[i]][2] <- log(mu)-(1/2)*log(1+(sigma^2/mu^2))
      prior.par.defs[[i]][3] <- sqrt(log(1+(sigma^2/mu^2)))
      
      names(prior.par.defs)[i] <- paste(names(prior.par.defs)[i],"transformed",sep="_")
      
    } else {
      
      if(prior.i[[1]][1]=="LognormalTrunc" |
           prior.i[[1]][1]=="lognormaltrunc" )
      {
        prior.par.defs[[i]][1] <- "NormalTrunc"
        
        mu    <- as.numeric(prior.i[[1]][2])
        sigma <- as.numeric(prior.i[[1]][3])
        
        prior.par.defs[[i]][2] <- log(mu)-(1/2)*log(1+(sigma^2/mu^2))
        prior.par.defs[[i]][3] <- sqrt(log(1+(sigma^2/mu^2)))
        prior.par.defs[[i]][4] <- log(as.numeric(prior.par.defs[[i]][4]))
        prior.par.defs[[i]][5] <- log(as.numeric(prior.par.defs[[i]][5]))
        
        names(prior.par.defs)[i] <- paste(names(prior.par.defs)[i],"transformed",sep="_")
      }
    }
  }
  return(prior.par.defs)
}


back.transform.parameters <- function(par)
{
  for(i in 1:length(par))
  {
    if(grepl("_transformed",names(par)[i])) 
    {   
      par[i] <- exp(par[i])
      names(par)[i] <- substr(names(par)[i], start=1, stop=nchar(names(par)[i])-12)
    }
  }  
  return(par)  
}
