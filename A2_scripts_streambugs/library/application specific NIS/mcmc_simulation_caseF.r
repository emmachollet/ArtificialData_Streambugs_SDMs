

# MCMC without likelihood, Marjoram et al. 2003 PNAS vol 100, no.26,p. 15324-15328, Case F

mcmc.simulation.caseF <- function(par.ini,f.dist, d.prior, r.prior, n.sample=10, eps=3,
                                  name.run="test1",Covar,beta,
                                  ...) 
{
  
  Covar.jump <- beta*Covar
  
  par <- par.ini
  par.samp <- matrix(par,nrow=1)
  colnames(par.samp) <- names(par)
  
  fileout <- paste("output/post_sample_temp",name.run,".dat",sep="")      
  
  cat(c(names(par),"log.prior.dens","dist","taxa false"),"\n",file=fileout,sep="\t")
  
  n.sample=n.sample
  
  j=0
  
  for ( i in 1:n.sample)
  {
    ## sample proposal parameter 
    par.new <- par + MASS::mvrnorm(1, mu=rep(0, length(par)), Sigma=Covar.jump)
    
    # calc dist
    
    dist.new <- f.dist(par.new,...) #...C=TRUE,verbose=FALSE,return.taxa.false
    cat("dist new:",dist.new,"\n")
    
    # calc prior dens
    
    log.prior.dens.new <- d.prior(par.new,log=TRUE)
    cat("log prior dens new:",log.prior.dens.new,"\n")
    
    # acceptance step:
    
    Prob.accept <- exp(log.prior.dens.new-log.prior.dens)
    
    if(runif(1) < Prob.accept & dist.new <= eps) 
    {
      par <- par.new
      dist <- dist.new
      log.prior.dens <- log.prior.dens.new
      j=j+1
    }
    
    par.samp <- rbind(par.samp,par)
    
    cat(c(par,log.prior.dens,dist),"\n",file=fileout,sep="\t",append=TRUE)
    
    cat("i:",i,"j:",j,"acc.rate:",j/i,"\n\n")
    
  }

  return(par.samp)
}
