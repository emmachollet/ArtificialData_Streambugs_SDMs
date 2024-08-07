###############################################################################
# modification 04-12-2014, AP,NS, provide acc.rate with adapt.beta    = FALSE
###############################################################################

# MCMC Metropolis-Hastings, Marjoram et al. 2003 PNAS vol 100, no.26,p. 15324-15328, Case E

mcmc.simulation.caseE <- function(par.ini,
                                  par.prior.delta,
                                  prior.par.defs,
                                  f.likeli, 
                                  d.prior, 
                                  n.sample=10,
                                  name.run="test1",
                                  Covar,
                                  beta,
                                  p.obs,
                                  p.abs,
                                  K.abs,
                                  D.drift,
                                  M.taxa     = M.taxa,
                                  y.names    = y.names,
                                  C          = TRUE,
                                  verbose    = FALSE,
                                  adapt.beta = FALSE,
                                  adapt.cov  = FALSE,
                                  name.with.pid = TRUE,
                                  ...) 
{
  
  if(name.with.pid) name.run <- paste(name.run,Sys.getpid(),sep="_")
  
  par <- par.ini
  par.samp <- matrix(par,nrow=1)
  colnames(par.samp) <- names(par)
  
  log.prior.dens <- d.prior(par,
                            par.prior.delta=par.prior.delta,
                            prior.par.defs=prior.par.defs,
                            log=TRUE)
  
  cat("log prior dens ini:",log.prior.dens,"\n")
  if(log.prior.dens==-Inf) stop ("initial prior density is -Inf, find better starting values" )
  
  res.ini <- f.likeli(par=par,y.names=y.names,C=C,verbose=verbose,par.prior.delta=par.prior.delta, 
                      M.taxa=M.taxa,p.obs=p.obs,p.abs=p.abs,K.abs=K.abs,D.drift=D.drift,
                      return.res=TRUE,...)
  
  log.likeli <- res.ini$loglikeli
  
  res.ss <- res.ini$res.ss
  
  fileout <- paste("output/post_par_sample_",name.run,".dat",sep="")
  fileout.sugpar <- paste("output/post_par_allsuggested_sample_",name.run,".dat",sep="")  
  fileout.res <- paste("output/post_res_sample_",name.run,".dat",sep="")  
  
  cat(c(names(par),"log.prior.dens","log.likeli","log.post","P.accept"),"\n",file=fileout,sep="\t") 
  cat(names(par),"\n",file=fileout.sugpar,sep="\t") 
  cat(names(res.ss),"\r",file=fileout.res,sep="\t")
  
  cat(c(par,log.prior.dens,log.likeli,log.prior.dens+log.likeli,1),"\r",file=fileout,sep="\t",append=TRUE) 
  cat(par,"\r",file=fileout.sugpar,sep="\t",append=TRUE) 
  
  cat(res.ss,"\r",file=fileout.res,sep="\t",append=TRUE)
  
  cat("log likeli ini:  ",log.likeli,"\n")
  cat("log post dens ini: ",log.prior.dens+log.likeli,"\n")
  
  
#   if(adapt.beta==TRUE)
#   {
    fileout.beta <- paste("output/beta_adapt_",name.run,".dat",sep="") 
    cat(c("beta","curr.acc.rate","tot.acc.rate","i"),"\n",file=fileout.beta,sep="\t")
    cat(c(beta,0,0,0),"\n",file=fileout.beta,sep="\t",append=TRUE)
#   }
  
  n.sample=n.sample
  
  Covar.jump <- beta*Covar
  
  j=0
  k=0
  l=0
  
  for ( i in 1:n.sample)
  {
    k=k+1
    
    ## sample proposal parameter 
    par.new <- par + MASS::mvrnorm(1, mu=rep(0, length(par)), Sigma=Covar.jump)
    
    cat(par.new,"\r",file=fileout.sugpar,sep="\t",append=TRUE)
    
    # calc prior dens
    
    log.prior.dens.new <- d.prior(par.new,par.prior.delta=par.prior.delta,prior.par.defs=prior.par.defs,log=TRUE)
    cat("log prior dens new:",log.prior.dens.new,"\n")
    
    if (log.prior.dens.new==-Inf) { Prob.accept <- 0 } else {
      
      # calc likelihood
      
      res.new <- f.likeli(par.new,y.names=y.names,C=C,verbose=verbose,par.prior.delta=par.prior.delta, 
                          M.taxa=M.taxa,p.obs=p.obs,p.abs=p.abs,K.abs=K.abs,D.drift=D.drift,
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
      l=l+1
    }
    
    par.samp <- rbind(par.samp,par)
    
    cat(c(par,log.prior.dens,log.likeli,log.prior.dens+log.likeli,Prob.accept),"\r",file=fileout,sep="\t",append=TRUE)
    
    cat(res.ss,"\r",file=fileout.res,sep="\t",append=TRUE)
    
    cat("i:",i,"j:",j,"acc.rate:",j/i,"\n")
    
    cat("k:",k,"l:",l,"current.acc.rate:",l/k,"\n\n")
    
    # update jump distance:
    
    if(k==50) 
    {      
      if(adapt.beta==TRUE)
      { 
        if(l/k > 0.6)  { beta <- 3*beta; cat("jumpdist updated: fact. 3 applied, beta:",beta,"\n") } else {
        if(l/k > 0.4)  { beta <- 1.5*beta; cat("jumpdist updated: fact. 1.5 applied, beta:",beta,"\n") } else {
        if(l/k < 0.05) { beta <- 0.5*beta; cat("jumpdist updated: fact. 0.5 applied, beta:",beta,"\n") } else {
        if(l/k < 0.12) { beta <- 0.7*beta; cat("jumpdist updated: fact. 0.7 applied, beta:",beta,"\n") } else {
        cat("jumpdist not updated, beta:",beta,"\n")
        
        } } } }

        Covar.jump <- beta*Covar
        
            }
      cat(c(beta,l/k,j/i,i),"\n",file=fileout.beta,sep="\t",append=TRUE)
      k=0
      l=0
    }
    
    # update covariance matrix and decrease correlation if necessary:
    
    if(adapt.cov ==TRUE & i%%825==0 )  # & i > 4000
    {
      Covar.new <- cov(par.samp)
      cat("determinant of Corr.new:",det(cov2cor(Covar.new)),"\n")
      
      while( det(cov2cor(Covar.new)) < 1e-10  )
      {  
        diagcn <- diag(Covar.new)
        Covar.new <- 0.99 * Covar.new
        diag(Covar.new) <- diagcn
      }
      
      Covar = Covar.new
      
      Covar.jump <- beta*Covar
      
      cat("Covar.jump updated, determinant of cor.mat:",det(cov2cor(Covar)),"\n")
      write.table(Covar.new,paste("output/Covar_last_",name.run,".dat",sep=""),sep="\t",row.names=TRUE,col.names=NA)
    }
  }
   
  return(par.samp)
}