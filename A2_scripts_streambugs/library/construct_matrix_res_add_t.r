construct.matrix.res.add.t <- function(res.add,t,y.names,file=NA)
  
  # t has to be a single timepoint that is part of tout
  # the function constructs a matrix with y.names as rows and 
  # additional output parameters as columns, 
  # foodweb parameters that depend on more than one y are excluded
  
{  
  if(length(t)>1) stop("error: more than one timepoint speciefied")
  
  #if( sum(dim(res.add)) < 1) res.add <- t(as.matrix(res.add))
  
  if(is.list(y.names)) y.names <- y.names$y.names
  
  rind <- which(res.add[,"time"]==t)
  
  if(sum(rind)==0) stop("time ",t," not available in res") 
  
  indfw <- grep("__",colnames(res.add))
  if(sum(indfw)>0) res.add <- res.add[,-indfw]
  
  if( sum(dim(res.add)) < 1) res.add <- t(as.matrix(res.add))
  
  #if(length(colnames(res.add)) ==0) stop ("error: colnames of res.add missing")
  
  cnames <- strsplit(colnames(res.add)[-1],split=".",fixed=TRUE)
  
  parnames <- NA
  for (i in 1:length(cnames))
  {
    parnames[i] <- cnames[[i]][1]
  }
  parnames <- sort(unique(parnames))
  
  res.add.t <- matrix(NA,
                      ncol=length(parnames),
                      nrow=length(y.names))
  
  rownames(res.add.t) <- y.names
  colnames(res.add.t) <- parnames
  
  for(j in 1:ncol(res.add.t))
  {
    for(i in 1:nrow(res.add.t))
    {
      cind <- which(colnames(res.add)==paste(parnames[j],y.names[i],sep="."))      
      
      if(sum(cind)>0)
      {
        res.add.t[i,j] <- res.add[rind,cind]
      }     
    }
  }
  
  y.names <- decode.statevarnames(y.names)
  
  if(!is.na(file)) write.table(res.add.t,file,sep="\t",row.names=TRUE,col.names=NA)
  
  return(res.add.t)
}