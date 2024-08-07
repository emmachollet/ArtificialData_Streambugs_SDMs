convert.CSusPOM <- function(par)
{
  class(par) <- "numeric"
  
  rnames <- strsplit(names(par),split="_")
  
  CSusPOM.ind <- NULL
  Fs.ind      <- which(names(par)=="Filt_scope")
  if(length(Fs.ind)==0) stop("parameter Filt_scope not found")
  
  for (i in 1:length(rnames))
  {
    if(rnames[[i]][length(rnames[[i]])]=="CSusPOM") CSusPOM.ind <- c(CSusPOM.ind,i)  
  }
  
  if(!length(CSusPOM.ind)>0) {warning("CSusPOM not found"); return(par)} else {
    
    # convert CSusPOM to DSusPOM
    
    par[CSusPOM.ind] <- par[Fs.ind] * par[CSusPOM.ind]  
    
    names(par)[CSusPOM.ind] <- sub("CSusPOM","DSusPOM",names(par)[CSusPOM.ind])
    
    return(par) }
}

revert.CSusPOM <- function(par)
{
  rnames <- strsplit(names(par),split="_")
  
  DSusPOM.ind <- NULL
  Fs.ind      <- which(names(par)=="Filt_scope")
  if(length(Fs.ind)==0) stop("parameter Filt_scope not found")
  
  for (i in 1:length(rnames))
  {
    if(rnames[[i]][length(rnames[[i]])]=="DSusPOM") DSusPOM.ind <- c(DSusPOM.ind,i)  
  }
  
  # convert DSusPOM to CSusPOM 
  
  par[DSusPOM.ind] <- par[DSusPOM.ind] / par[Fs.ind] 
  
  names(par)[DSusPOM.ind] <- sub("DSusPOM","CSusPOM",names(par)[DSusPOM.ind])
  
  return(par)
}