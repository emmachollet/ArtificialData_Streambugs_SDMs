generate.par.samp.matrix <- function(par.unc,n.samp=10,verbose=FALSE)
{
  
  par.samp           <- matrix(NA,nrow=length(par.unc),ncol=n.samp+1) 
  rownames(par.samp) <- rep(NA,length(par.unc))
  
  if(verbose == TRUE)
  {
    print("Generating parameters sampling")
    pb <- txtProgressBar(min=0, max=length(par.unc), char = "=", label = "", style = 3)
  }
  
  for ( p in 1:length(par.unc))
  {
    # print(p)
    if(length(par.unc[[p]]) > 1) # parameter defined as distribution
    {
      
      # choose mean value as first value in the parameter sample
      
      if(par.unc[[p]][1]=="discrete"|par.unc[[p]][1]=="Discrete") 
        # search for the class with highest probality
      {
        nprobs <- (length(par.unc[[p]])-1)/2 
        probs <- as.numeric(par.unc[[p]][(length(par.unc[[p]])-nprobs+1):length(par.unc[[p]])]) 
        classmax <-par.unc[[p]][which(probs==max(probs))+1]
        if(length(classmax)>1) 
        {
          classmax <- classmax[1]
          warning("Discrete parameter ",names(par.unc[p]),"has no unique max, first class is used" )
        }
        par.samp[p,1]  <- classmax
      } else {
        
        if(par.unc[[p]][1]=="uniform"|par.unc[[p]][1]=="Uniform") 
          #take mean between min and max
        {
          minv <- as.numeric(par.unc[[p]][2])
          maxv <- as.numeric(par.unc[[p]][3])
          par.samp[p,1] <- mean(c(minv,maxv))
        } else {
          
          par.samp[p,1]  <- par.unc[[p]][2]
        }
      }
      
      # sampling 
      
      par.samp[p,-1] <- t(sysanal.randsamp(sampsize=n.samp,dist="Indep",distdef=par.unc[p])$sample)
      
      rownames(par.samp)[p] <- names(par.unc)[p]
      
    } else {
      
      par.samp[p,1]  <- par.unc[[p]]
      par.samp[p,-1] <- rep(par.unc[[p]],times=n.samp)
      rownames(par.samp)[p] <- names(par.unc)[p]
    } 
    
    if(verbose == TRUE)
    {
      setTxtProgressBar(pb, p)
    }
#      if(verbose == TRUE) cat(p,"  ")
    
  }
  
  return(par.samp)
  
}