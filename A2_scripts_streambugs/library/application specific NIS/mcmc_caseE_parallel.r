MCMC.caseE.parallel = function (par.unc, par.prior.delta=par.prior.delta,prior.par.defs,
                                f.likeli, d.prior, n.sample, name.run, Covar, beta, C=TRUE, verbose=TRUE,  
                                y.names   = y.names,
                                M.taxa    = M.taxa,
                                p.obs     = p.obs,
                                p.abs     = p.abs,
                                D.drift   = D.drift,
                                update    = FALSE,
                                n.chain   = 4,
                                n.cpu,
                                libraries =NULL,
                                ...)
{
  if ( !require(parallel) )
  {
    install.packages("parallel")
    library(parallel)
  }
  
  cl <- makeCluster(min(n.cpu, detectCores()))
  on.exit({
    stopCluster(cl)
    print("Cluster stopped.")
  })
  
  varlist <- unique(c(ls(), ls(envir = .GlobalEnv), ls(envir = parent.env(environment()))))
  clusterExport(cl, varlist = varlist, envir = environment())
  clusterSetRNGStream(cl)
  wd <- getwd()
  clusterExport(cl, varlist = c("libraries", "wd"), 
                envir = environment())
  
  
  if(update == FALSE)
  {
    MCMC.wrap <- function(x, ...) 
    {
      require(Matrix)
      if (!is.null(libraries)) 
        source(libraries)
      
      mcmc.simulation.caseE.parallel(
        par.unc = par.unc,
        name.run = name.run[1],
        par.prior.delta = par.prior.delta,
        prior.par.defs=prior.par.defs,
        f.likeli = f.likeli,
        d.prior = d.prior,
        n.sample = n.sample,
        Covar = Covar,
        beta = beta,
        M.taxa    = M.taxa,
        p.obs     = p.obs,
        p.abs     = p.abs,
        D.drift   = D.drift,        
        C = C,
        verbose = verbose,
        y.names = y.names,
        update = update,
        ...)
    }
    
    result <- clusterApply(cl, x=1:n.chain, fun=MCMC.wrap)
    #, par.unc = par.unc,par.prior.delta=par.prior.delta, f.likeli = f.likeli, d.prior = d.prior,
    #  n.sample = n.sample, name.run = name.run, Covar = Covar, beta = beta, p1 = p1, p2 = p2,
    # C = C, verbose = verbose, update = update
    
  } else {
    if(length(name.run) != n.chain)
    {
      warning("In updating mode, the number of indicated run names (name.run vector) outrule the number of chains (n.chain)")
    }
    
    par.list <- list()
    for(i in 1:length(name.run))
    {
      filein <- paste("output/post_par_sample_",name.run[i],".dat",sep="") 
      
      par.samp.old <- read.table(filein, h=T, sep="\t")[,names(par.unc)]
      par.old <- as.numeric(as.vector(as.matrix(par.samp.old[nrow(par.samp.old),])))
      names(par.old) <- names(par.unc)
      
      par.list[[i]] <- list(par.unc = par.old, name.run = name.run[i])
    }
    
    MCMC.wrap <- function(x, par.list,...) 
    {
      require(Matrix)
      if (!is.null(libraries)) source(libraries)
      mcmc.simulation.caseE.parallel(par.unc   = par.list[[x]][["par.unc"]],
                                     name.run  = par.list[[x]][["name.run"]],
                                     par.prior.delta = par.prior.delta,
                                     prior.par.defs  = prior.par.defs,
                                     f.likeli  = f.likeli,
                                     d.prior   = d.prior,
                                     n.sample  = n.sample,
                                     Covar     = Covar,
                                     beta      = beta,
                                     M.taxa    = M.taxa,
                                     p.obs     = p.obs,
                                     p.abs     = p.abs,
                                     D.drift   = D.drift,
                                     C         = C,
                                     verbose   = verbose,
                                     y.names   = y.names,
                                     update    = update,
                                     ...)
    }
    
    result <- clusterApply(cl=cl, x=1:length(name.run), fun=MCMC.wrap, par.list = par.list #, par.prior.delta = par.prior.delta, 
                           #                            f.likeli = f.likeli, d.prior = d.prior, 
                           #                            n.sample = n.sample, Covar = Covar, beta = beta,p1 = p1, p2 = p2, 
                           #                            C = C, verbose = verbose, update = update, 
                           #                            libraries = libraries
    )
    
    
  }
  return(result)
}
