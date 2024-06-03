get.factors.inverts.wrong.predicted <- function(y.res,res.add.t,file.fact=NA,sort="taxa") #sort="sites"
  
{

  flim <- c("ffoodlim","fselfinh","fcurrent","ftempmax","fmicrohab","forgmicropoll","fsapro")
  
  res.add.t <- res.add.t[,flim]
  
  if(sort=="taxa")
  {
    y.over <- decode.statevarnames(y.res$y.over)
    ind.taxa.sort <- sort(y.over$y.taxa,index.return=TRUE)$ix
    y.res$y.over <- y.res$y.over[ind.taxa.sort]

    
    y.under<- decode.statevarnames(y.res$y.under)
    ind.taxa.sort <- sort(y.under$y.taxa,index.return=TRUE)$ix
    y.res$y.under <- y.res$y.under[ind.taxa.sort]
    
    y.corr <- decode.statevarnames(y.res$y.corr)
    ind.taxa.sort <- sort(y.corr$y.taxa,index.return=TRUE)$ix
    y.res$y.corr <- y.res$y.corr[ind.taxa.sort]
  } 
    
  res.add.over <- res.add.t[y.res$y.over,]
  res.add.under <- res.add.t[y.res$y.under,]
  res.add.corr  <- res.add.t[y.res$y.corr,]
    
  
  res.add.over.c <- cbind(rownames(res.add.over),res.add.over)
  colnames.over <- c("overest. taxa",colnames(res.add.over))
  
  res.add.under.c <- cbind(rownames(res.add.under),res.add.under)
  colnames.under <- c("underest. taxa",colnames(res.add.under))
  
  res.add.corr.c <- cbind(rownames(res.add.corr),res.add.corr)
  colnames.corr <- c("corr. taxa",colnames(res.add.corr))
  
  res.add.all <- rbind(colnames.under,res.add.under.c,
                       colnames.over,res.add.over.c,
                       colnames.corr,res.add.corr.c)
  
  if(!is.na(file.fact)) write.table(res.add.all,file.fact,sep="\t",row.names=FALSE,col.names=FALSE)
  
  return(res.add.all)
  
}