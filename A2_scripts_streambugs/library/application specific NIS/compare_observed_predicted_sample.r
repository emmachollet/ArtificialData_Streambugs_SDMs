#compare survivals and observed

compare.observed.predicted.sample <- function(p.y.mat,observed.abund,Invertebrates)
{  
  
  n.opt <- calc.max.correct.modres(observed.abund,Invertebrates)$n.opt
  
  
  # 1. calculating the fraction of parameter samples that lead to correct results
  
  f.correct <- rep(NA,times=ncol(p.y.mat))
  names(f.correct) <- colnames(p.y.mat)
  
  for(j in 1:ncol(p.y.mat))
  {
    ind.correct <- which(p.y.mat[,j]>0.5)
     
    f.correct[j] <- length(ind.correct)/nrow(p.y.mat)
    
  }
 
  # 2. calculating the number of observations where more than 50% of the model realisations 
  #    lead to correct results 
  
  n.cor <- length(which(f.correct>0.5))
  n.incor <- length(which(f.correct<=0.5))
  
  cat("analyzing fraction of correct results per observation: \n")
  
  cat("> 50% n.obs correct:",n.cor,"\n" )
  
  cat("> 50% n.obs incorrect:",n.incor,"\n")
  
  cat("fraction of max number of taxa at sites and times correct:",signif(n.cor/n.opt,2),"\n\n" )
  
  # 3. averaging the observation probability and then check if the mean p.y is in accordance with data
  
  compare.mod.obs <- observed.abund[,-ncol(observed.abund)]
  compare.mod.obs[,-(1:4)] <- NA 
  
  mean.p.1 <- colMeans(p.1.mat)
  
  for(i in 1:length(mean.p.1))
  {
    
    site <- unlist(strsplit(names(mean.p.1)[i],split="_"))[1]
    hab  <- unlist(strsplit(names(mean.p.1)[i],split="_"))[2]
    date <- unlist(strsplit(names(mean.p.1)[i],split="_"))[5]
    
    taxon <- unlist(strsplit(names(mean.p.1)[i],split="_"))[3]
    
    rind.obs <- which(observed.abund[,"ReachID"]==site & 
                        observed.abund[,"Habitat"]==hab &
                        observed.abund[,"Date"]==date)
    cind.obs <- which(colnames(observed.abund)==taxon)
    
    if(!is.na(rind.obs))
    {
      
      if(observed.abund[rind.obs,cind.obs] > 0  & mean.p.1[i] > 0.5) compare.mod.obs[rind.obs,cind.obs] <- 0 # correct (observed and high surv prob)
      if(observed.abund[rind.obs,cind.obs] == 0 & mean.p.1[i] <= 0.5) compare.mod.obs[rind.obs,cind.obs] <- 0 # correct (not observed and low surv prob)
      if(observed.abund[rind.obs,cind.obs] > 0  & mean.p.1[i] <= 0.5) compare.mod.obs[rind.obs,cind.obs] <--1 # underestimated by the model: observed but low surv prob
      if(observed.abund[rind.obs,cind.obs] == 0 & mean.p.1[i] > 0.5) compare.mod.obs[rind.obs,cind.obs] <- 1 # overestimated by the model: not observed but high surv. prob
    } 
  }
  
  n.over   <- length(which(compare.mod.obs==1))
  n.und  <- length(which(compare.mod.obs==-1))
  n.corr <- length(which(compare.mod.obs==0))
  
  cat("analyzis of mean p.obs: \n")
  cat("no of taxa at sites and times over est based on mean p.obs:",n.over,"\n")
  cat("no of taxa at sites and times under est based on mean p.obs:",n.und,"\n")
  cat("no of taxa at sites and times correct based on mean p.obs:",n.corr,"\n\n")
  
  p.corr <- n.corr/(n.corr+n.over+n.und)
  p.over <- n.over/(n.corr+n.over+n.und)
  p.und  <- n.und/(n.corr+n.over+n.und)
  
  cat("fraction of taxa at sites and times over est based on mean p.obs:",round(p.over,digits=2),"\n")
  cat("fraction of taxa at sites and times under est based on mean p.obs:",round(p.und,digits=2),"\n")
  cat("fraction of taxa at sites and times correct based on mean p.obs:",round(p.corr,digits=2),"\n\n")
  
  cat("fraction of max number of taxa at sites and times correct based on mean p.obs:",round(n.corr/n.opt,digits=2),"\n\n")
  
  return(list("f.correct"=f.correct,
              "compare.mod.obs"=compare.mod.obs))     
}
