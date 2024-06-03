calc.loglikelihood <- function(res.ss,y.names,M.taxa,observed.abund,p.obs,p.abs,K.abs,D.drift)
  
{
  loglikeli <- 0
  
  p.0 <- NULL
  p.1 <- NULL
  p.y <- NULL
  
  for(i in 1:length(res.ss))
  {
    ind.y <- which(y.names$y.names==names(res.ss)[i])
    
    group <- y.names$y.groups[ind.y]
    
    if(group=="Invertebrates")
    {
      site <- y.names$y.reaches[ind.y]
      hab  <- y.names$y.habitats[ind.y]
      tax  <- y.names$y.taxa[ind.y]
       
      ind.j <- which(observed.abund[,"ReachID"] == site & observed.abund[,"Habitat"]==hab)
      
      dates <- observed.abund[ind.j,"Date"]
      
      M.i <- M.taxa[names(res.ss)[i]] # gDM
      
      if(length(dates)>0)
      {
        for(k in ind.j)
        {
          A.jk <- observed.abund[k,"area_m2"]
          
          if(sum(A.jk)==0) warning("area not given, A.jk:",A.jk)
          
          if(res.ss[i]<0) 
          {
            if(res.ss[i] < -1e-7) warning("res.ss of ",names(res.ss)[i]," <0: ",res.ss[i]," set to zero")
            res.ss[i] <- 0
          }
          
          p.absent <-  p.abs*exp(-log(2)*(res.ss[i]/M.i)/K.abs)
          
          p.ijk.0 <- p.absent + (1-p.absent) * (1-p.obs)^(A.jk*((res.ss[i]/M.i)+D.drift))
          p.ijk.1 <- 1 - p.ijk.0
          
          p.0 <- c(p.0, p.ijk.0)
          names(p.0)[length(p.0)] <- paste(y.names$y.names[ind.y],observed.abund[k,"Date"],sep="_")
          
          p.1 <- c(p.1, p.ijk.1)
          names(p.1)[length(p.1)] <- paste(y.names$y.names[ind.y],observed.abund[k,"Date"],sep="_")
                    
          if(observed.abund[k,tax]>0) p.ijk <- p.ijk.1 else {
            if(observed.abund[k,tax]==0) p.ijk <- p.ijk.0 else {warning(y.names$y.names,"not observed at",observed.abund[k,"Date"] )}}
          
          p.y <- c(p.y,p.ijk)
          names(p.y)[length(p.y)] <- paste(y.names$y.names[ind.y],observed.abund[k,"Date"],sep="_")
          
#           cat(log(p.ijk),names(M.i),i,k,"\n")
          
          loglikeli <- loglikeli + log(p.ijk)
        
        }        
      }
    }
  }
  
return(list("loglikeli"=loglikeli,"p.0"=p.0,"p.1"=p.1,"p.y"=p.y))
}