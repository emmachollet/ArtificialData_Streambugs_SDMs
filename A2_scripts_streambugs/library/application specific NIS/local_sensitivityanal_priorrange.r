


local.sensitivityanal <- function(par.ini,
                                  par.prior.delta,
                                  prior.par.defs,
                                  f.likeli, 
                                  d.prior, 
                                  n.samp   = 10,
                                  name.run = "locsens_test1",
                                  inc.perc = 1,
                                  p.obs,
                                  p.abs,
#                                   K.abs,
                                  D.drift,
                                  M.taxa   = M.taxa,
                                  y.names  = y.names,
                                  C        = TRUE,
                                  verbose  = FALSE,
                                  sd.par   = NA ,   # has priority over inc.perc
                                  ...) 
{
  
  par <- par.ini[1:length(prior.par.defs)]
  par.bt <- back.transform.parameters(par) 
  
  if(!sum(is.na(sd.par))>0)  # if not na
  {
    sd.par <- sd.par[names(par)]
    sd.par.bt <- par.bt*sqrt(exp(sd.par^2)-1) 
  }
  
  log.prior.dens <- d.prior(par,
                            par.prior.delta=par.prior.delta,
                            prior.par.defs=prior.par.defs,
                            log=TRUE)
  
  if(log.prior.dens==-Inf) stop ("initial prior density is -Inf, find better starting values" )
  
  res.ini <- f.likeli(par=par,y.names=y.names,C=C,verbose=verbose,par.prior.delta=par.prior.delta, 
                      M.taxa=M.taxa,p.obs=p.obs,p.abs=p.abs,D.drift=D.drift,
                      return.res=TRUE,...)
  
  log.likeli <- res.ini$loglikeli
  
  res.ss <- res.ini$res.ss
  res.py <- res.ini$p.y
  
  fileout <- paste("output/post_par_sample_",name.run,".dat",sep="")      
  fileout.res <- paste("output/post_res_sample_",name.run,".dat",sep="")  
  fileout.py <- paste("output/post_res_py_",name.run,".dat",sep="")
  
  cat(c(names(par),"log.prior.dens","log.likeli","log.post"),"\n",file=fileout,sep="\t") 
  cat(names(res.ss),"\r",file=fileout.res,sep="\t")
  cat(names(res.py),"\r",file=fileout.py,sep="\t")
  
  cat(c(par,log.prior.dens,log.likeli,log.prior.dens+log.likeli),"\r",file=fileout,sep="\t",append=TRUE) 
  cat(res.ss,"\r",file=fileout.res,sep="\t",append=TRUE)
  cat(res.py,"\r",file=fileout.py,sep="\t",append=TRUE)
  
  pdf(paste("output/res_sensanal_",name.run,".pdf",sep=""),width=11.69,height=16.54)#,width=7,height=7)
  par(mfrow=c(7,3),mar=c(4,4,3,1),mgp=c(2,1,0))  #c(bottom, left, top, right) ,las=1
  
  for(i in 1:length(par))
  {
    par.new <- par
    if(sum(is.na(sd.par))>0)
    {
      par.i.v <- par[i] + par[i] * sort(c((inc.perc/100)*-1:-(n.samp/2),(inc.perc/100)*1:(n.samp/2)))
    } else {
      par.i.v <- log(par.bt[i] + sort(c( (3*sd.par.bt[i]/n.samp)*-1:-(n.samp/2),(3*sd.par.bt[i]/n.samp)*1:(n.samp/2))))
    }
    log.post.i <- rep(NA,times=length(par.i.v))
    #     log.post.i[1] <- log.prior.dens+log.likeli
    
    log.prior.i <- rep(NA,times=length(par.i.v))
    #     log.prior.i[1] <- log.prior.dens
    
    log.likeli.i <- rep(NA,times=length(par.i.v))
    #     log.likeli.i[1] <- log.likeli
    
    for(j in 1:length(par.i.v))
    {
      par.new <- par
      par.new[i] <- par.i.v[j]
      
      log.prior.dens.new <- d.prior(par.new,par.prior.delta=par.prior.delta,prior.par.defs=prior.par.defs,log=TRUE)
      
      cat(names(par)[i],j,": value:  ",par.i.v[j],"\n")
      cat(names(par)[i],j,": log prior dens new:",log.prior.dens.new,"\n")
      
      res.new <- f.likeli(par.new,y.names=y.names,C=C,verbose=verbose,par.prior.delta=par.prior.delta, 
                          M.taxa=M.taxa,p.obs=p.obs,p.abs=p.abs,D.drift=D.drift,
                          return.res=TRUE,...) #
      
      log.likeli.new <- res.new$loglikeli
      
     
      cat(names(par)[i],j,": log likeli new:  ",log.likeli.new,"\n")
      
      cat(names(par)[i],j,": log post dens new: ",log.prior.dens.new+log.likeli.new,"\n")
      
      cat(c(par.new,log.prior.dens.new,log.likeli.new,log.prior.dens.new+log.likeli.new),"\r",file=fileout,sep="\t",append=TRUE)
      
      cat(res.new$res.ss,"\r",file=fileout.res,sep="\t",append=TRUE)
      cat(res.new$p.y,"\r",file=fileout.py,sep="\t",append=TRUE)
      
      log.post.i[j] <- log.prior.dens.new+log.likeli.new
      log.prior.i[j] <- log.prior.dens.new
      log.likeli.i[j] <- log.likeli.new
    }
    
    par.i.v.bt <- exp(par.i.v) 
    
    plot(par.i.v.bt,log.prior.i,main=names(par)[i],type="p",xlab="",ylab="log prior",pch=19,col=1)
    points(par.bt[i],log.prior.dens,pch=19,col=2)
    
    plot(par.i.v.bt,log.likeli.i,main=names(par)[i],type="p",xlab="",ylab="log likelihood",pch=19,col=1)
    points(par.bt[i],log.likeli,pch=19,col=2)
    
    plot(par.i.v.bt,log.post.i,main=names(par)[i],type="p",xlab="",ylab="log posterior",pch=19,col=1)
    points(par.bt[i],(log.prior.dens+log.likeli),pch=19,col=2)
  }
  
  dev.off()
  
  return()
}