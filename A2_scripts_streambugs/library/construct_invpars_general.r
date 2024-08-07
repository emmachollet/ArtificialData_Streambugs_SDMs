construct.invpars.general <- function(inverts)
{
    
    invpars <- list()
    ind <- 0
    
    for (i in 1:length(inverts) )
    {
      invert <- inverts[i]
      
      ind <- ind+1
      invpars[[ind]] <- list( dist="NormalTrunc", dist.par= c(1,0.25,0,1e10) )   
      names(invpars)[ind] <- paste(invert,"_fbasali",sep="")    
      
      ind <- ind+1
      invpars[[ind]] <- list( dist="NormalTrunc", dist.par= c(1,0.25,0,1e10) )   
      names(invpars)[ind] <- paste(invert,"_fconsi",sep="")
      
      ind <- ind+1
      invpars[ind] <- 1    
      names(invpars)[ind] <- paste("unc_",invert,"_hdens",sep="")  
      
      ind <- ind+1
      invpars[ind] <- 1    
      names(invpars)[ind] <- paste("unc_",invert,"_Kfood",sep="")
      
      
    }
   
    return(invpars)
  
}