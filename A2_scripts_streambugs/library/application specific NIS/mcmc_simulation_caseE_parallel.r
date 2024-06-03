# MCMC Metropolis-Hastings, Marjoram et al. 2003 PNAS vol 100, no.26,p. 15324-15328, Case E
# In contrary to Cedrics function this function does always start a new file for parameters instead of updating the old parameter file, 
# Furthermore, it writes model results to a file.
# Note the different naming of par.fix=par.prior.delta, as well as different arguments for the likelihoodfunction
# Note that the parfile in my case has less extra columns (4 instead of 7)

mcmc.simulation.caseE.parallel <- function(par.unc,
                                           par.prior.delta,
                                           prior.par.defs,
                                           f.likeli, 
                                           d.prior, 
                                           n.sample = 10,
                                           name.run = "test1",
                                           Covar,
                                           beta,
                                           p.obs,
                                           p.abs,
                                           D.drift,
                                           C=C,
                                           verbose = verbose,
                                           y.names = y.names,
                                           M.taxa = M.taxa, 
                                           update = FALSE,
                                           ...)
  {
  
  if(update == FALSE)
  {
    name.run2 = paste(name.run, Sys.getpid(), 1, sep="_")
    
    fileout <- paste("output/post_par_sample_",name.run2,".dat",sep="") 
    fileout.res <- paste("output/post_res_sample_",name.run2,".dat",sep="")   
    
    while(file.exists(fileout) | file.exists(fileout.res))
    {
      name.run2 = paste(name.run2,"new",Sys.Date(),sep="_")
      
      fileout <- paste("output/post_par_sample_",name.run2,".dat",sep="") 
      fileout.res <- paste("output/post_res_sample_",name.run2,".dat",sep="")   
    }
    
  } else {
    
    fileout     <- paste("output/post_par_sample_",name.run,"_update.dat",sep="") 
    fileout.res <- paste("output/post_res_sample_",name.run,"_update.dat",sep="")  
    
    while(file.exists(fileout) | file.exists(fileout.res))
    {
      name.run2 = paste(name.run,"new",Sys.Date(),sep="_")
      
      fileout <- paste("output/post_par_sample_",name.run,"_update_new_",Sys.Date(),".dat",sep="") 
      fileout.res <- paste("output/post_res_sample_",name.run,"_update_new_",Sys.Date(),".dat",sep="") 
    }
  }
  
  par.samp <- matrix(par.unc,nrow=1)
  colnames(par.samp) <- names(par.unc)
  
  # calculate initial results and write them to files
  
  log.prior.dens <- d.prior(par.unc,
                            par.prior.delta=par.prior.delta,
                            prior.par.defs=prior.par.defs,
                            log=TRUE)
  
  if(log.prior.dens==-Inf) stop ("initial prior density is -Inf, find better starting values" )
  
  res.ini <- f.likeli(par= par.unc,y.names=y.names,C=C,verbose=verbose,par.prior.delta=par.prior.delta, 
                      M.taxa=M.taxa,p.obs=p.obs,p.abs=p.abs,D.drift=D.drift,
                      return.res=TRUE,...)
  
  log.likeli <- res.ini$loglikeli
  
  res.ss <- res.ini$res.ss
  
  cat(c(names(par.unc),"log.prior.dens","log.likeli","log.post","P.accept"),"\n",file=fileout,sep="\t") 
  
  cat(c(par.unc,log.prior.dens,log.likeli,log.prior.dens+log.likeli,1),"\r",file=fileout,sep="\t",append=TRUE) 
  
  cat(names(res.ss),"\r",file=fileout.res,sep="\t")
  
  cat(res.ss,"\r",file=fileout.res,sep="\t",append=TRUE)
  
  # start MCMC 
  
  Covar.jump <- beta*Covar
  
  j=0
  
  par <- par.unc
  
  for ( i in 1:n.sample)
  {
    ## sample proposal parameter 
    par.new <- par + MASS::mvrnorm(1, mu=rep(0, length(par)), Sigma=Covar.jump)
    
    # calc prior dens
    
    log.prior.dens.new <- d.prior(par.new,
                                  par.prior.delta=par.prior.delta,
                                  prior.par.defs=prior.par.defs,
                                  log=TRUE)
    cat("log prior dens new:",log.prior.dens.new,"\n")
    
    if (log.prior.dens.new==-Inf) { Prob.accept <- 0 } else {
      
      # calc likelihood
      
      res.new <- f.likeli(par.new,y.names=y.names,C=C,verbose=verbose,par.prior.delta=par.prior.delta, 
                          M.taxa=M.taxa,p.obs=p.obs,p.abs=p.abs,D.drift=D.drift,
                          return.res=TRUE,...)
      
      log.likeli.new <- res.new$loglikeli
      
      cat("log likeli new:  ",log.likeli.new,"\n")
      
      cat("log post dens new: ",log.prior.dens.new+log.likeli.new,"\n")
      
      # acceptance step:
      
      Prob.accept <- min(1,exp(log.prior.dens.new+log.likeli.new-log.prior.dens-log.likeli))
      
    }
    
    if(runif(1) < Prob.accept) 
    {
      par <- par.new
      log.prior.dens <- log.prior.dens.new
      log.likeli     <- log.likeli.new
      res.ss         <- res.new$res.ss
      j=j+1
    }
    
    par.samp <- rbind(par.samp,par)
    
    cat(c(par,log.prior.dens,log.likeli,log.prior.dens+log.likeli,Prob.accept),"\r",file=fileout,sep="\t",append=TRUE)
    
    cat(res.ss,"\r",file=fileout.res,sep="\t",append=TRUE)
    
    cat("i:",i,"j:",j,"acc.rate:",j/i,"\n\n")
    
  }
  
  return(par.samp)
}
