construct.observation.matrix <- function(data,taxa)
{
  # input: 
  #  data:  table with taxa in columns and sites in rows and abundance as entries
  #  taxa:  vector with names of invertebrate taxa
  # output: matrix with taxa in columns, sites/habitats in rows 
  #         and "always"/"never"/"sometimes" as entries
  
  # example: 
  # data <- read.table("input/observed_taxa.dat",sep="\t",header=TRUE)
  # taxa <- y.names$taxa  
  
  sites      <- unique(data[,"ReachID"])
  habitats   <-unique(data[,"Habitat"])
  sites.habs <- unique(paste(data[,"ReachID"],data[,"Habitat"],sep="."))
  dates      <- unique(data[,"Date"])
  

  observed.matrix <- matrix(NA,nrow=length(sites.habs),ncol=length(taxa))
  colnames(observed.matrix) <- taxa
  rownames(observed.matrix) <- sites.habs
  
  for (i in 1:length(sites.habs))
  {
    site <- unlist(strsplit(sites.habs[i],split="\\."))[1]
    hab  <- unlist(strsplit(sites.habs[i],split="\\."))[2]
  #  observed.matrix[i,"ReachID"] <- site
  #  observed.matrix[i,"Habitat"] <- hab
    
    rind <- data[,"ReachID"]==site & data[,"Habitat"]==hab
    
    
    for (j in 1:length(taxa))
    {
      cind <- which(taxa[j]==colnames(data))
      if(sum(cind)>0)
      {
        nobs <- sum (data[rind,cind] > 0 )
        
      } else {nobs <- NA}  
      
      if(is.na(nobs)) {flag="NA"} else
      {
        if( nobs == sum(rind) ) {flag="always"} else
        {
          if( nobs == 0 ) {flag="never"} else
          {
            if(nobs > 0 & nobs < sum(rind) ) {flag<-"sometimes"}  else
            {
              warning(taxa[j]," cannot be flagged")
            }
          } 
        }
      }
      observed.matrix[i,j] <- flag
    }
  }
  
  return(observed.matrix)
  
}  

check.observed <- function(observed.matrix,y.names)
{
 
  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)
  
  c.observed <- c(rep(NA,length(y.names$y.names)))
  names(c.observed) <- y.names$y.names
  
  for ( i in 1:length(c.observed))
  {
    y.substr <- unlist(strsplit(names(c.observed)[i],split="\\_"))
    
    reach     <- y.names$reaches[grep(paste(y.substr,collapse="|"),y.names$reaches)]
    hab       <- y.names$habitats[grep(paste(y.substr,collapse="|"),y.names$habitats)]
    group     <- y.names$groups[grep(paste(y.substr,collapse="|"),y.names$groups)]
    taxon     <- y.substr[c(-match(group,y.substr),-match(reach,y.substr),-match(hab,y.substr))]

    indr  <- grep(reach,rownames(observed.matrix))
    indh  <- grep(hab,rownames(observed.matrix))
    rind <- intersect(indr,indh)
    cind  <- which(taxon==colnames(observed.matrix))
    if (sum(rind)>0 & sum(cind)>0)
    {
      c.observed[i] <- observed.matrix[rind,cind]
    } else 
    {
      c.observed[i] <- NA
    }
    
  }

  return(c.observed)
  
}


