calc.max.correct.modres <- function(observed.abund,Invertebrates)
{
  sites <- sort(unique(observed.abund[,"ReachID"]))
  habitats <- sort(unique(observed.abund[,"Habitat"]))
  #invertebrates <- colnames(observed.abund)[c(-(1:4),-ncol(observed.abund))]
  
  n.opt.sum=0
  n.tot.sum=0
  
  for(i in 1:length(sites))
  {
    for(k in 1:length(habitats))
    {
      rind <- which(observed.abund[,"ReachID"]==sites[i] & observed.abund[,"Habitat"]==habitats[k] )
      
      for(j in 1:length(Invertebrates))
      {
        cind <- which(colnames(observed.abund)==Invertebrates[j])
        
        n.obs    <- sum(observed.abund[rind,cind]>0)
        n.notobs <- sum(observed.abund[rind,cind]==0)
        n.tot    <- n.obs+n.notobs
        
        if(n.tot != length(!is.na(observed.abund[rind,cind]))) warning("total number of obs not correct")
        
        n.opt <- max(n.obs,n.notobs)
        
        n.opt.sum <- n.opt.sum + n.opt
        n.tot.sum <- n.tot.sum + n.tot
      }
    }
  }
  
  cat("max no of correctly modelled observations:", n.opt.sum,"\n")
  cat("total no of observations:", n.tot.sum,"\n")
  cat("max fraction of correctly modelled observations:",round(n.opt.sum/n.tot.sum,digits=2),"\n\n")
  
  return(optim.mod.res <- list("n.opt"=n.opt.sum,
                               "n.tot"=n.tot.sum))
}