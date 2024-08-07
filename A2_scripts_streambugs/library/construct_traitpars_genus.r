construct.traitpars.genus <- function ( invert,cind,name.trait,method=max, #method=max,mean,min,sum,... #Achtung, auch mean normiert nicht auf summe=10
                                        db.traits=db.traits,
                                        db.taxonomy=db.taxonomy,... )
{ 
  # cind <- 3:15
  # name.trait <- "microhab"  
  
  # invert <- "Nele"
  # invert <- "Gammarus"
  
  
  # matching and set parameters
  
  par <- numeric(0)
  genus <- numeric(0)
  
  rindgat <- which(db.taxonomy[,"Genus"] == invert)
  
  if (sum(rindgat)==0)                            # invert not given at genus level
  {
    rindspec <- which(db.taxonomy[,"Species2"] == invert) 
    
    if (length(rindspec)==1) # invert given at species level and species found
    {
      genus   <- as.character(db.taxonomy[rindspec,"Genus"])
      species <- as.character(db.taxonomy[db.taxonomy[,"Genus"] %in% genus, "Species2"])
    } 
  }
  
  if (sum(rindgat)>0)  #invert given at genus level
  {
    genus   <- as.character(unique(db.taxonomy[rindgat,"Genus"]))
    species <- as.character(db.taxonomy[db.taxonomy[,"Genus"] %in% genus, "Species2"])
  }  
  
  if(length(genus)>1) 
  {
    warning("Genus ", genus ," for",invert," not unique")  # something would have gone wrong in this case!
  } 
  
  if(length(genus)<1)   # genus not matched
  { 
    if ( substr(invert,start=nchar(invert)-1,stop=nchar(invert)) != "ae" )
    {
      cat("  note: genus for",invert," not found, search on family level \n")
    }
  } else 
  {
    
    rind <- which(rownames(db.traits) %in% c(genus, species, paste(species, "Ad.", sep=""),
                                             paste(species, "Lv.", sep="")))
    
    if ( length(rind) > 0 ) # genus matched 
    {
      
      if ( sum( is.na(db.traits[rind,cind])/length(rind) ) < length(cind) & 
            (length(cind)== 1 | sum(db.traits[rind,cind], na.rm=TRUE) > 0) )          # not all traitlevels NA or 0   # xxx problem either with spear or tachet   
        
        
        {  
        
        par[paste(invert,name.trait,colnames(db.traits[cind]),sep="_")] <- 0
        
        for( j in cind )  
        {
          if( sum(is.na(db.traits[rind,j])) < length(rind)    )
          {  
            par[paste(invert,name.trait,colnames(db.traits)[j],sep="_")] <- method(db.traits[rind,j],na.rm=TRUE)
          }
        }
        
      } else # taxon uniquely matched but all entries are NA or 0
      {
        if (  sum( is.na(db.traits[rind,cind])/length(rind) ) == length(cind) |    #xxxProblem with either Spear or Tachet database
               (length(cind)>1 & sum(db.traits[rind,cind], na.rm=TRUE) == 0)  )
        {
          cat("  note:",invert,"found but all",name.trait,"levels are NA, search on family level \n") 
        }  
      }    
    }    
  }
  if( length(par) < 1 ) par=NA
  
  return(par) 
  
}  

#  invert<- "Isoperladifformis"
#  invert<- "Isoperlawieauchimmer"
#  invert<- "Nele"
#  invert<- "Ephydatia"
#  invert<- "Isoperla" 

#  parmax <- construct.traitpars.genus(invert,cind<-3:15,name.trait="microhab" )
#  parmean <- construct.traitpars.genus(invert,cind<-3:15,name.trait="microhab",method=mean )