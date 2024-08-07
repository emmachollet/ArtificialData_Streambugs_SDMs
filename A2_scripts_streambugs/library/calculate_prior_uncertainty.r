# calculate results for parameter sample:
# =======================================

ptm <- proc.time()
par.defs.ini = prior.par.defs

for(i in names(par.ini))
{
  if(as.numeric(par.defs.ini[[i]][2]) != par.ini[i])
  {
    par.defs.ini[[i]][2] = as.character(par.ini[i])
  }
}

par.samp <- generate.par.samp.matrix(n.samp=n.samp,par.unc=par.defs.ini)

res.array <- array(NA,dim=c(length(tout),length(y.names$y.names)+1,n.samp+1), #time,states,samples
                   dimnames=list(1:length(tout),c("tout",y.names$y.names),
                                 1:(n.samp+1)))
for (n in 1:ncol(par.samp) )
{
  if(n%%50 == 0) save.image(file = paste("streambugs_",name.run,".RData",sep=""))
  st <- proc.time()
  
  par.fix <- as.numeric(par.samp[,n])
  names(par.fix) <- rownames(par.samp)
  
  par.fix <- c(par.fix, par.prior.delta)
  par.fix <- back.transform.parameters(par.fix)
  
  # calc stoichiometric parameters
  
  par.stoich.out <- calc.stoich(par=as.list(par.fix),returns="parout")
  
  # assign calculated group stoichiometric parameters to the taxa:
  
  par.stoich.taxa <- assign.par.stoich(par.invtraits,par.stoich.out,y.names)
  
  # combine parameters
  
  par.fix <- c(par.fix,par.stoich.taxa)
  
  # calculate and plot results if stoich is ok, else write NAs
  
  if( round(par.stoich.out["Death_Algae_O2"],digits=10) < 0 |  
        round(par.stoich.out["Death_Invertebrates_O2"],digits=10) < 0 )
  {
    res.array[, ,n]  <- NA
    res.array[,1,n]  <- tout
    
  } else
  {
    
    # calculate results
    
    #     res.array[,,n] <- run.streambugs(y.names  = y.names, #y.names$y.names,
    #                                      times    = tout,
    #                                      par      = par.fix,
    #                                      inp      = inp,
    #                                      C        = TRUE)$res
    
    for ( i in 1:length(y.names$reaches))
    {
      reach <- y.names$reaches[i]
      ind.reach <- y.names$y.reaches==reach
      
      res.i <- run.streambugs(y.names = y.names$y.names[ind.reach],
                              times   = tout,
                              par     = par.fix,
                              inp     = inp,
                              C       = TRUE,
                              verbose=FALSE)$res
      
      
      if(i==1) res.n <- res.i else  
      { 
        if (nrow(res.n)==nrow(res.i))  res.n <- cbind(res.n,res.i) else
          warning("nrow(res.n) not equal nrow(res.i")
      } 
      
    }
    
    ind.time <- which(colnames(res.n)=="time")
    ind.time <- ind.time[-1]
    if(length(ind.time) > 0)
    {
      res.n      <- res.n[,-ind.time]
    } 
  
    if( sum( dim(res.array[,,n])==dim(res.n) )==2 ) 
    {
      res.array[,,n] <- res.n
    } else
    {
      res.array[, ,n]  <- NA
      res.array[,1,n]  <- tout
    }
    
  }
  
  cat("simulation no", n, "done in", (proc.time() - st)[1:3], "@", date(),"\n" )
}
print(proc.time()-ptm)


# remove NAs in pars and res

ind.na <- numeric(0)
for (i in 1:dim(res.array)[3])
{
  if(is.na(sum(res.array[,,i])))
  {
    ind.na <- c(ind.na,i)
  }
}    

if(length(ind.na)>0)
{  
  res.array     <- res.array[,,-ind.na ]     #time,states,samples
  par.samp.nona <- par.samp[,-ind.na]        
  #par.samp.uni.nona <- par.samp.uni[,-ind.na]
} else
{
  par.samp.nona <- par.samp
}  

#Number of simulations
paste("No of Simulations: ",dim(res.array)[3])

save.image(file = paste("streambugs_",name.run,".RData",sep=""))


# calculate quantiles of results:

res.quant.05 <- as.matrix(res.array[,,1])
res.quant.50 <- as.matrix(res.array[,,1])
res.quant.95 <- as.matrix(res.array[,,1])

for ( j in 2:ncol(res.quant.05) )
{
  for ( i in 1:nrow(res.quant.05) )
  {
    res.quant.05[i,j] <- quantile( res.array[i,j,] , probs = 0.05,na.rm=T)
    res.quant.50[i,j] <- quantile( res.array[i,j,] , probs = 0.5, na.rm=T)
    res.quant.95[i,j] <- quantile( res.array[i,j,] , probs = 0.95,na.rm=T)
  }        
}


save.image(file = paste("streambugs_",name.run,".RData",sep=""))

