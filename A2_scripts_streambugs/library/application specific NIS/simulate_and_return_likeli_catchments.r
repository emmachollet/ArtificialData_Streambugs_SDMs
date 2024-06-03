simulate.streambugs.return.likeli.catchments <- function(par,y.names=y.names,C=TRUE,verbose=TRUE,
                                                         par.prior.delta=par.prior.delta,
                                                         M.taxa=M.taxa,
                                                         return.res=FALSE,p.obs,p.abs,K.abs,D.drift,
                                                         catchments,Invertebrates.list)
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
      reach <- y.names$reaches[i]
      
      catchment <- NA
      for(k in 1:length(catchments))
      {
        if (grepl(catchments[k],reach)) catchment <- catchments[k] 
        stop
      }
      
      ind.reach   <- y.names$y.reaches==reach
      ind.inverts <- y.names$y.taxa %in% Invertebrates.list[[catchment]]
      ind.tax <- y.names$y.groups!="Invertebrates"
      
      y.names.i <- y.names$y.names[ind.reach & (ind.inverts|ind.tax)]
      
      y.names.i <- decode.statevarnames(y.names.i)
      
      par.stoich.taxa <- assign.par.stoich(par.invtraits,par.stoich.out,y.names.i,
                                           CPOM.dif=CPOM.dif,ratio_pred_prey=par.fix["ratio_pred_prey"])
      
      # combine parameters
      
      par.fix.i <- c(par.fix,par.stoich.taxa)
      
      res.i <- run.streambugs(y.names = y.names.i, 
                              times   = tout,
                              par     = par.fix.i,
                              inp     = inp,
                              C       = C,
                              verbose = verbose)$res
      
      # calc mean res steadystate
      res.ss.i <- calc.res.steadystate(res=res.i,times=tout,par=par.fix.i,y.names=y.names.i)
      
      # calc loglikelihood
      reslikeli.i <- calc.loglikelihood(res.ss=res.ss.i,y.names=y.names.i,M.taxa=M.taxa,observed.abund=observed.abund,
                                        p.obs=p.obs,p.abs=p.abs,K.abs=K.abs,D.drift=D.drift) 
      
      if (length(tout) > nrow(res.i)) 
      {
        stop("!!! simulations for reach ",reach," not completed, nrow(res) < length(tout)")
      }
      if(i==1) 
      {
        res.ss    <- res.ss.i
        reslikeli <- reslikeli.i 
      } else { 
        res.ss              <- c(res.ss,res.ss.i) 
        
        reslikeli$loglikeli <- reslikeli$loglikeli+reslikeli.i$loglikeli
        reslikeli$p.0       <- c(reslikeli$p.0,reslikeli.i$p.0) 
        reslikeli$p.1       <- c(reslikeli$p.1,reslikeli.i$p.1) 
        reslikeli$p.y       <- c(reslikeli$p.y,reslikeli.i$p.y) 
      }
    }
    
    reslikeli$res.ss <- res.ss
    
  }
  
  if(return.res) 
  {  
    return(reslikeli)
  } else return(sum(reslikeli$loglikeli))  
}
