#compare survivals and observed

compare.observed.predicted <- function(p.y,observed.abund,Invertebrates)
{  
  
  #   source("libraries/calc_max_correct_modres.r")
  n.opt <- calc.max.correct.modres(observed.abund,Invertebrates)$n.opt
  
  n.corr <- length(which(p.y>0.5))
  
  compare.mod.obs <- observed.abund[,-ncol(observed.abund)]
  compare.mod.obs[,-(1:4)] <- NA 
  
  for (i in 1:length(p.y) )
  {
    site <- unlist(strsplit(names(p.y)[i],split="_"))[1]
    hab  <- unlist(strsplit(names(p.y)[i],split="_"))[2]
    date <- unlist(strsplit(names(p.y)[i],split="_"))[5]
    
    taxon <- unlist(strsplit(names(p.y)[i],split="_"))[3]
    
    rind.obs <- which(observed.abund[,"ReachID"]==site & 
                        observed.abund[,"Habitat"]==hab &
                        observed.abund[,"Date"]==date)
    cind.obs <- which(colnames(observed.abund)==taxon)
    
    if(!is.na(rind.obs))
    {
      
      if(observed.abund[rind.obs,cind.obs] > 0  & p.y[i] > 0.5) compare.mod.obs[rind.obs,cind.obs] <- 0 # correct (observed and high surv prob)
      if(observed.abund[rind.obs,cind.obs] == 0 & p.y[i] > 0.5) compare.mod.obs[rind.obs,cind.obs] <- 0 # correct (not observed and low surv prob)
      if(observed.abund[rind.obs,cind.obs] > 0  & p.y[i] <= 0.5) compare.mod.obs[rind.obs,cind.obs] <--1 # underestimated by the model: observed but low surv prob
      if(observed.abund[rind.obs,cind.obs] == 0 & p.y[i] <= 0.5) compare.mod.obs[rind.obs,cind.obs] <- 1 # overestimated by the model: not observed but high surv. prob
    }
  }
  
  n.over   <- length(which(compare.mod.obs==1))
  n.und  <- length(which(compare.mod.obs==-1))
  n.corr <- length(which(compare.mod.obs==0))
  
  cat("no of taxa at sites and times over est:",n.over,"\n")
  cat("no of taxa at sites and times under est:",n.und,"\n")
  cat("no of taxa at sites and times correct:",n.corr,"\n\n")
  
  p.corr <- n.corr/(n.corr+n.over+n.und)
  p.over <- n.over/(n.corr+n.over+n.und)
  p.und  <- n.und/(n.corr+n.over+n.und)
  
  cat("fraction of taxa at sites and times over est:",round(p.over,digits=2),"\n")
  cat("fraction of taxa at sites and times under est:",round(p.und,digits=2),"\n")
  cat("fraction of taxa at sites and times correct:",round(p.corr,digits=2),"\n\n")
  
  cat("fraction of max number of taxa at sites and times correct:",round(n.corr/n.opt,digits=2),"\n\n")
  
  y.over <- NULL
  y.under <- NULL
  y.corr <- NULL
  
  for(i in 1:nrow(compare.mod.obs))
  {
    indo <- which(compare.mod.obs[i,]==1)
    if(length(indo)>0) y.over <- c(y.over,paste(compare.mod.obs[i,"ReachID"],
                                                compare.mod.obs[i,"Habitat"],
                                                colnames(compare.mod.obs)[indo],
                                                "Invertebrates",sep="_"))
    
    indu <- which(compare.mod.obs[i,]==-1)
    if(length(indu)>0) y.under <- c(y.under,paste(compare.mod.obs[i,"ReachID"],
                                                  compare.mod.obs[i,"Habitat"],
                                                  colnames(compare.mod.obs)[indu],
                                                  "Invertebrates",sep="_"))
    
    indc <- which(compare.mod.obs[i,]==0)
    if(length(indc)>0) y.corr <- c(y.corr,paste(compare.mod.obs[i,"ReachID"],
                                                compare.mod.obs[i,"Habitat"],
                                                colnames(compare.mod.obs)[indc],
                                                "Invertebrates",sep="_"))
  }
  
  y.res <- list("y.over"=y.over,"y.under"=y.under,"y.corr"=y.corr,"n.opt"=n.opt)
  
  return(list("compare.mod.obs"=compare.mod.obs,"y.res"=y.res))     
}
