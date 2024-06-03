simulation.parsamp <- function(par.samp,par.prior.delta,f.survivals,name.run="test1",tout=tout,
                               y.names=y.names,M.taxa=M.taxa,
                               p.obs=p.obs,p.abs=p.abs,D.drift=D.drift,outputfolder="output") 
{
  
  # par.samp with rows = samples and columns = parameters
  
  n.samp=ncol(par.samp)
  
  fileout.res <- paste(outputfolder,"/res_ss_sample_",name.run,".dat",sep="") 
  fileout.py  <- paste(outputfolder,"/res_py_sample_",name.run,".dat",sep="") 
  fileout.p1  <- paste(outputfolder,"/res_p1_sample_",name.run,".dat",sep="") 
  fileout.p0  <- paste(outputfolder,"/res_p0_sample_",name.run,".dat",sep="") 
  
  for (n in 1:n.samp )
  {
    
    st <- proc.time()
    if(n%%50 == 0) save.image(file = paste("streambugs_",name.run,".RData",sep=""))
    
    par <- par.samp[,n]
    
    res <- f.survivals(y.names=y.names,par=par,par.prior.delta=par.prior.delta,verbose=FALSE,return.res=TRUE,
                       M.taxa=M.taxa,p.obs=p.obs,p.abs=p.abs,D.drift=D.drift)
    
    res.D <- c(res$res.ss,res$loglikeli)
    names(res.D) <- c(names(res$res.ss),"loglikeli")
    
    if(n==1) 
    {
      cat(names(res.D),"\n",file=fileout.res,sep="\t")
      cat(names(res$p.y),"\n",file=fileout.py,sep="\t")
      cat(names(res$p.1),"\n",file=fileout.p1,sep="\t")
      cat(names(res$p.0),"\n",file=fileout.p0,sep="\t")
      res.mat <- matrix(NA,ncol=length(res.D),nrow=n.samp)
      colnames(res.mat) <- names(res.D)
    }
    
    cat(res.D,"\n",file=fileout.res,sep="\t",append=TRUE)
    cat(res$p.y,"\n",file=fileout.py,sep="\t",append=TRUE)
    cat(res$p.1,"\n",file=fileout.p1,sep="\t",append=TRUE)
    cat(res$p.0,"\n",file=fileout.p0,sep="\t",append=TRUE)
    
    res.mat[n,] <- res.D
    
    cat("simulation no", n, "done in", (proc.time() - st)[1:3], "@", date(),"\n" )
    
  }
  
  return(res.mat)
  
}