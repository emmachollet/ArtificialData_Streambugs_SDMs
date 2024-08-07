construct.traitpars.species <- function(invert,cind,name.trait,
                                        db.traits=db.traits,
                                        db.taxonomy=db.taxonomy,
                                        method=max,...)
{
  # matching and set parameters
  
  par <- numeric(0)
  
  rind <- which(rownames(db.traits) %in% c(invert, paste(invert, "Ad.", sep=""), paste(invert, "Lv.", sep=""))) # species level found or genus matched uniquely
  
  
  #   if(sum(rind)<1) rind <- grep(invert,rownames(db.traits))
  
  if ( length(rind) == 1 ) # species level found or genus matched uniquely
  {
    
    if ( sum( is.na(db.traits[rind,cind]) ) < length(cind) &                     # not all traitlevels NA (freshwater data base)
         ( length(cind)== 1 | sum(db.traits[rind,cind], na.rm=TRUE) > 0   ) )        # not all trait levels are 0 (Tachet data base)  # xxxxx nis Problem either with spear or tachet db!
                  
    {  
      for( j in cind )  
      {          
        
        ifelse(is.na(db.traits[rind,j]),
               par[paste(invert,name.trait,colnames(db.traits[j]),sep="_")] <- 0,
               par[paste(invert,name.trait,colnames(db.traits[j]),sep="_")] <- 
                 db.traits[rind,j])
      }
    } else # taxon uniquely matched but all entries are NA
    {
      if ( sum( is.na(db.traits[rind,cind]) ) == length(cind) |
           (length(cind)> 1 &  sum(db.traits[rind,cind], na.rm=TRUE) == 0)) #all traitlevels NA or 0
      {
        cat("  note:",invert,"found but all",name.trait,"levels are NA, search on genus level \n")
      }  
    }    
  }    
  
  if ( length(rind) > 1) 
  {
    if ( sum( is.na(db.traits[rind,cind])/length(rind) ) < length(cind)  & 
           sum(db.traits[rind,cind], na.rm=TRUE) > 0)  
    {
      par[paste(invert,name.trait,colnames(db.traits[cind]),sep="_")] <- 0
      
      for( j in cind )  
      {
        if( sum(is.na(db.traits[rind,j])) < length(rind)    )
        {  
          par[paste(invert,name.trait,colnames(db.traits)[j],sep="_")] <- method(db.traits[rind,j],na.rm=TRUE)
        }
      }
    }
  }
  
  if( length(par) < 1 ) par=NA
  
  return(par) 
  
}  

#  invert<- "Isoperladifformis"
#  invert<- "Nele"
#  invert<- "Ephydatiamuelleri"
#  invert <- "Elmisaenea"
# 
#  par <- construct.traitpars.species(invert,cind<-3:15,name.trait="microhab" )