divide.parameters <- function(par.unc=par.unc)
  # 1.
  # divides parameter list "par.unc" that contains of values as well as distribution definitions
  # into a vector with fixed parameters and a list with parameter distribution definitions
  
{
  prior.delta <- NULL 
  ind.delta   <- NULL
  
  prior.par.defs  <- list()
  
  for ( p in 1:length(par.unc))
  {
    
    if(length(par.unc[[p]]) == 1) 
    {
      if (!is.na(par.unc[[p]])) 
      {
        prior.delta <- c(prior.delta,par.unc[[p]])
        names(prior.delta)[length(prior.delta)] <- names(par.unc[p])
        ind.delta <- c(ind.delta,p)
      } 
    } else {
      
      if(length(par.unc[[p]]) > 1)
      {
        if(par.unc[[p]][1]=="delta"|par.unc[[p]][1]=="Delta")
        {
          prior.delta <- c(prior.delta,as.numeric(par.unc[[p]][2]))
          names(prior.delta)[length(prior.delta)] <- names(par.unc[p])
          ind.delta <- c(ind.delta,p)
        }       
      }
    }
  }
  
  ind.defs <- setdiff(c(1:length(par.unc)),ind.delta)
  ind <- 0
  
  for( p in ind.defs )
  {   
    
    if(length(par.unc[[p]]) > 1)
    {
      if(( par.unc[[p]][1]=="Normal"         | par.unc[[p]][1]=="normal"|
             par.unc[[p]][1]=="NormalTrunc"    | par.unc[[p]][1]=="normaltrunc"|
             par.unc[[p]][1]=="Lognormal"      | par.unc[[p]][1]=="lognormal"|
             par.unc[[p]][1]=="LognormalTrunc" | par.unc[[p]][1]=="lognormaltrunc") &
           par.unc[[p]][3]==0)
      {
        prior.delta <- c(prior.delta,as.numeric(par.unc[[p]][2]))
        names(prior.delta)[length(prior.delta)] <- names(par.unc[p])
        
      } else {
        
        if(  par.unc[[p]][1]=="uniform"          | par.unc[[p]][1]=="Uniform"|
               par.unc[[p]][1]=="Normal"         | par.unc[[p]][1]=="normal"|
               par.unc[[p]][1]=="NormalTrunc"    | par.unc[[p]][1]=="normaltrunc"|
               par.unc[[p]][1]=="Lognormal"      | par.unc[[p]][1]=="lognormal"|
               par.unc[[p]][1]=="LognormalTrunc" | par.unc[[p]][1]=="lognormaltrunc"|
               par.unc[[p]][1]=="Inv"            | par.unc[[p]][1]=="inv"|
               par.unc[[p]][1]=="Exponential"    | par.unc[[p]][1]=="exponential"|
               par.unc[[p]][1]=="Discrete"       | par.unc[[p]][1]=="discrete")
        {
          ind=ind+1
          prior.par.defs[ind] <- par.unc[p]
          names(prior.par.defs)[ind] <- names(par.unc[p])
          
        } else {
          
          cat(names(par.unc[p]),": dist",par.unc[[p]][1] , "not yet implemented\n")
        }  
      }
    } 
  }
  return(list("prior.delta"=prior.delta,"prior.par.defs"=prior.par.defs) ) 
}
