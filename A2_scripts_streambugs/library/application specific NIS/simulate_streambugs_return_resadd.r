simulate.streambugs.return.resadd <- function(par,y.names=y.names,C=TRUE,verbose=TRUE,par.prior.delta=par.prior.delta,
                                              M.taxa=M.taxa,
                                              p.obs,p.abs,D.drift,
                                              tout.add=tout[length(tout)],
                                              outputfolder="output",
                                              observed.abund)
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
    stop("death consumes oxygen, res=NA")
    
  } else {
    
    loglikeli <- NA
    
    res.both <- run.streambugs(y.names  = y.names,
                               times    = tout,
                               par      = par.fix,
                               inp      = inp,
                               C        = TRUE,
                               #file.add = paste(outputfolder,"/streambugs_res_add_",name.run,".dat",sep=""),
                               return.res.add=TRUE,
                               tout.add=tout[length(tout)])
    
    res     <- res.both$res
    res.add <- res.both$res.add
    
    
    if (length(tout) > nrow(res)) 
    {
      stop("!!! simulations not completed, nrow(res) < length(tout)")
    }
    
    res.add.t <- construct.matrix.res.add.t(res.add,t=tout.add,
                                            y.names,
                                            file=paste(outputfolder,"/res_add_t_",name.run,".dat",sep=""))
    
    
  }
  
  # calc mean res steadystate
  
  res.ss <- calc.res.steadystate(res=res,times=tout,par=par.fix,y.names=y.names)
  
  # calc loglikelihood
  
  reslikeli <- calc.loglikelihood(res.ss=res.ss,y.names=y.names,M.taxa=M.taxa,observed.abund=observed.abund,
                                  p.obs=p.obs,p.abs=p.abs,D.drift=D.drift) 
  
  reslikeli$res.ss <- res.ss
  
  comp.obs.pred <- compare.observed.predicted(reslikeli$p.y,observed.abund,Invertebrates)
  
  y.res <- comp.obs.pred$y.res
  
  res.add.all <- get.factors.inverts.wrong.predicted(y.res,
                                                     res.add.t,
                                                     file.fact=paste(outputfolder,"/res_add_sorted_",name.run,".dat",sep=""),
                                                     sort="taxa")
  
  return(res.add.all)
  
}
