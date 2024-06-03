construct.traitpars <- function(invert,cind,name.trait,method=max, #method=max,mean,min... #Achtung, auch mean normiert nicht auf summe=10
                                db.traits=db.traits, default=NA,
                                db.taxonomy=db.taxonomy )
{
    
  # cind <- 3:15
  # name.trait <- "microhab"  
  # method="max"
  # invert <- "Nele"
  # invert <- "Gammarus"
  # invert <- "Isoperladifformis"  
  # invert <- "Isoperla"
  
  par<-numeric(0)
  
  parspecies <- construct.traitpars.species(invert,cind,name.trait,
                                            db.traits=db.traits,
                                            db.taxonomy=db.taxonomy)
  if( sum(!is.na(parspecies))>0  )
  {
    par <- parspecies
  } else
  {
    pargenus <- construct.traitpars.genus(invert,cind,name.trait,method=method,
                                          db.traits=db.traits,
                                          db.taxonomy=db.taxonomy) 
    if( sum(!is.na(pargenus))>0  )
    {
      par<-pargenus
    } else
    {
      parfamily <- construct.traitpars.family(invert,cind,name.trait,method=method,default=default,
                                              db.traits=db.traits,
                                              db.taxonomy=db.taxonomy)
      if( sum(!is.na(parfamily))>0  )
      {
        par <- parfamily
      }
      
    }  
  }   
  return(par)
}  

#  invert<- "Rhithrogenacircumtatrica"
#  invert<- "Nele"
#  invert<- "Ephydatia"
#  invert<- "Isoperla" 

#  parmax <- construct.traitpars(invert,cind=3:15, name.trait="microhab" )
#  parmean <- construct.traitpars(invert,cind=3:15, name.trait="microhab",method="mean" )
#  parmn <- construct.traitpars(invert,cind=3:15, name.trait="microhab",method="mn" )




