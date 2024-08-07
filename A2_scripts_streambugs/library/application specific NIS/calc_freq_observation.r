calc.freq.observation <- function(observed.abund,Invertebrates)
{
  sites <- sort(unique(observed.abund[,"ReachID"]))
  habitats <- sort(unique(observed.abund[,"Habitat"]))
  
  #invertebrates <- colnames(observed.abund)[c(-(1:4),-ncol(observed.abund))]
  
  freq.obs <- matrix(NA,nrow=length(sites)*length(habitats),ncol=length(Invertebrates))
  colnames(freq.obs) <- Invertebrates
  rownames(freq.obs) <- apply(expand.grid(sites, habitats), 1, paste, collapse="_")
  
  for(i in 1:length(sites))
  {
    for(k in 1:length(habitats))
    {
      
      rind <- which(observed.abund[,"ReachID"]==sites[i] & observed.abund[,"Habitat"]==habitats[k] )
      
      for(j in 1:length(Invertebrates))
      {
        cind <- which(colnames(observed.abund)==Invertebrates[j])
        
        rind.freq.obs <- which(rownames(freq.obs)==paste(sites[i],habitats[k],sep="_"))
        
        freq.obs[rind.freq.obs,j] <- sum(observed.abund[rind,cind]>0)/length(!is.na(observed.abund[rind,cind]))
        
      }
    }
  }
  
  return(freq.obs)
}