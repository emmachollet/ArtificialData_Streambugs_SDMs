# 01.03.23 K.abs added to arguments and calc likelihood

simulate.streambugs.return.likeli <- function(par,y.names=y.names,C=TRUE,verbose=TRUE,
                                              par.prior.delta=par.prior.delta,
                                              M.taxa=M.taxa,
                                              return.res=FALSE,p.obs,p.abs,D.drift,K.abs)
{
  
  if(!is.list(y.names)) y.names <- decode.statevarnames(y.names)
  
  par.fix <- c(par,par.prior.delta) 
  
  par.fix <- back.transform.parameters(par.fix)
  
  par.fix <- convert.CSusPOM(par.fix)
  
  # run streambugs:
  # ===============
  
  # calc stoichiometric parameters
  
  par.stoich.out <- calc.stoich(par=as.list(par.fix),returns="parout")
  
  # assign calculated group stoichiometric parameters to the taxa:
  
  par.stoich.taxa <- assign.par.stoich(par.invtraits,par.stoich.out,y.names,CPOM.dif=CPOM.dif)
  
  # combine parameters
  
  par.fix <- c(par.fix,par.stoich.taxa)
  
  reslikeli <- list()
  
  # calculate results if stoich is ok, else write NAs
  
  if(CPOM.dif) 
  {
    Death_Invertebrates_O2 <- par.stoich.out["Death2_Invertebrates_O2"]
  } else {
    Death_Invertebrates_O2 <- par.stoich.out["Death1_Invertebrates_O2"]
  }
  
  if( round(par.stoich.out["Death_Algae_O2"],digits=10) < 0 |  
        round(Death_Invertebrates_O2,digits=10) < 0 )
  {
    reslikeli$loglikeli <- -Inf
    reslikeli$p.0 <- NA
    reslikeli$p.1 <- NA
    reslikeli$p.y <- NA
    reslikeli$res.ss <- rep(NA,times=length(y.names$y.names))
    names(reslikeli$res.ss) <- y.names$y.names
    
  } else {
    
    loglikeli <- NA
    
    # calculate results for each reach
    
    for ( i in 1:length(y.names$reaches))
    {
      # i <- 1
      # verbose <- T
      # C <- T
      reach <- y.names$reaches[i]
      ind.reach <- y.names$y.reaches==reach
      
      res.i <- run.streambugs(y.names = y.names$y.names[ind.reach], 
                              times   = tout,
                              par     = par.fix,
                              inp     = inp,
                              C       = C,
                              verbose=verbose)$res
      
      
      if (length(tout) > nrow(res.i)) 
      {
        stop("!!! simulations for reach ",reach," not completed, nrow(res) < length(tout)")
      }
      if(i==1) res <- res.i else  res <- cbind(res,res.i)
    }
    
    # remove doubled time-columns if necessary:
    ind.time <- which(colnames(res)=="time")
    if(length(ind.time)>1)
    {
      ind.time <- ind.time[-1]
      res      <- res[,-ind.time]
    }
    
    # calc mean res steadystate
    
    res.ss <- calc.res.steadystate(res=res,times=tout,par=par.fix,y.names=y.names)
    
    # calc loglikelihood
    
    reslikeli <- calc.loglikelihood(res.ss=res.ss,y.names=y.names,M.taxa=M.taxa,observed.abund=observed.abund,
                                    p.obs=p.obs,p.abs=p.abs,D.drift=D.drift,K.abs=K.abs) 
    
    reslikeli$res.ss <- res.ss
    
  }
  
  if(return.res) 
  {  
    return(reslikeli)
  } else return(reslikeli$loglikeli)  
}
