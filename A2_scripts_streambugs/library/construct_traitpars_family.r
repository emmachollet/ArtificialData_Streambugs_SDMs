construct.traitpars.family <- function(invert, cind, name.trait,method=max, #method=max,mean,min,sum,... #Achtung, auch mean normiert nicht auf summe=10                   
                                       default=1,
                                       db.traits=db.traits,
                                       db.taxonomy=db.taxonomy,... )
{
  # cind <- 3:15
  # name.trait <- "microhab"  
  # default=1 
  # invert <- "Nele"
  # invert <- "Isoperladifformis"  
  
  # matching 
  
  par <- numeric(0)
  family <- numeric(0)
  genera <- numeric(0)
  
  rindfam <- which(db.taxonomy[,"Family"] == invert)
  
  if(length(rindfam)>0)
  {
    family <- unique(db.taxonomy[rindfam,"Family"])[1]
  } else
  {
    rindgat <- which(db.taxonomy[,"Genus"] == invert)
    
    if(length(rindgat)>0)
    {
      family <- unique(db.taxonomy[rindgat,"Family"])[1]
    } else
    {
      rindspec <- which(db.taxonomy[,"Species2"] == invert)
      
      if(length(rindspec)>0)
      {
        family <- unique(db.taxonomy[rindspec,"Family"])[1]
      }
    }  
  }
  
  if (length(family)>0)  
  {
    rind   <- which(db.taxonomy[,"Family"] == family)
    genera <- unique(db.taxonomy[rind,"Genus"])
    species <- unique(db.taxonomy[rind, "Species2"])
  } 
  
  if (length(genera)==0) 
  {   
    #cat("  note:",name.trait,"for",invert,"not found \n") 
  } else  
    
  {   
    rind <- which(rownames(db.traits) %in% c(genera, species, 
                                             paste(species,"Ad.",sep=""),
                                             paste(species,"Lv.",sep="")))
    
    if ( length(rind) > 0 ) # genera matched 
    {
      
      if ( sum(is.na(db.traits[rind,cind]))/length(rind) < length(cind)   & 
            ( length(cind)==1 | sum(db.traits[rind,cind], na.rm=TRUE) > 0) )         # not all traitlevels NA or 0
      {  
        par[paste(invert,name.trait,colnames(db.traits[cind]),sep="_")] <- 0
        
        for( j in cind )  
        {          
          if( sum(is.na(db.traits[rind,j])) < length(rind) ) 
          {   
            par[paste(invert,name.trait,colnames(db.traits)[j],sep="_")] <- method(db.traits[rind,j],na.rm=TRUE)
          }  
        }  
        
        if(length(rindfam)<1) 
        {
          cat("  note:",name.trait,"for",invert,"only found on family level \n")
        }
        
      } else # taxon uniquely matched but all entries are NA
      {
        if (  sum( is.na(db.traits[rind,cind]))/length(rind)  == length(cind)    |   #xxx Problem with either spear or tachet database
                (length(cind)>1 & sum(db.traits[rind,cind], na.rm=TRUE) == 0 ))
        {
          cat("  note:",invert," only found on family level and all",name.trait,"levels are NA \n") 
        }  
      }    
    }    
  } 
  
  if( length(par) < 1 ) 
  {
    for( j in cind )
    {
      par[paste(invert,name.trait,colnames(db.traits[j]),sep="_")] <- default
    }   
    
    cat("warning:",name.trait,"for",invert,"not found,levels set to ",default," \n")
  }   
  
  return(par) 
}   


#  invert<- "Isoperladifformis"
#  invert<- "Nele"
#  invert<- "Ephydatiamuelleri"
# 
#  par <- construct.traitpars.family(invert,cind<-3:15,name.trait="microhab" )


