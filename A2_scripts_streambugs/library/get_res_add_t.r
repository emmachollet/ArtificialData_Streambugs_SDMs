get.res.add.t.invertebrates <- function(res.add.t,y.names.invertebrates,par=par)
  # get "limitation" factors for specific invertebrates
{
  y.names <- rownames(res.add.t)
  
  y.names <- decode.statevarnames(y.names)
  
  res.add.t.Invertebrates <- matrix(NA,nrow=length(y.names.invertebrates),ncol=9)
  colnames(res.add.t.Invertebrates) <- c("fselfinh", "ffoodlim","fcurrent","ftempmax","fmicrohab",
                                         "forgmicropoll","fsapro","fgrotax","fbasaltax")
  
  for(i in 1:length(y.names.invertebrates))
  {
    rind <- which(y.names$y.names==y.names.invertebrates[i])
    
    if (sum(rind)<1) stop("y.names.invertebrates not part of rownames(res.add.t)")
    
    res.add.t.Invertebrates[i,1:7] <- signif(res.add.t[rind,colnames(res.add.t.Invertebrates)[1:7]],digits=3)
    
    ind.par.fgrotax <- which(names(par)==paste(y.names$y.taxa[rind],"fgrotax",sep="_"))
    ind.par.fbasaltax <- which(names(par)==paste(y.names$y.taxa[rind],"fbasaltax",sep="_"))
    
    res.add.t.Invertebrates[i,8:9] <- signif(par[c(ind.par.fgrotax,ind.par.fbasaltax)],digits=3)
    
  }
  rownames(res.add.t.Invertebrates) <- y.names.invertebrates
  
  return(res.add.t.Invertebrates)
}

# Algae: 

get.res.add.t.algae <- function(res.add.t,y.names.algae)
  # get "limitation" factors for algae taxa
{
  y.names <- rownames(res.add.t)
  
  y.names <- decode.statevarnames(y.names)
  
  res.add.t.Algae <- matrix(NA,nrow=length(y.names.algae),ncol=7)
  colnames(res.add.t.Algae) <- c("flimI","flimnutrients","flimN","flimP","fselfshade","forgmicropoll","fsapro")
  
  for(i in 1:length(y.names.algae))
  {
    rind <- which(y.names$y.names==y.names.algae[i])
    
    if (sum(rind)<1) stop("y.names.algae not part of rownames(res.add.t)")
    
    res.add.t.Algae[i,] <- signif(res.add.t[rind,colnames(res.add.t.Algae)],digits=3)
    
  }
  rownames(res.add.t.Algae) <- y.names.algae
  
  return(res.add.t.Algae)
}