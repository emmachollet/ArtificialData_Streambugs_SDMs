# 27.2.2023 adapted to run based on streambugs package, see lines 17 and 55

assign.par.stoich <- function(par.invtraits,par.stoich.out,y.names,Pref.default=1)
{
  
  if ( !is.list(y.names) ) y.names <- decode.statevarnames(y.names)
  
  ind.inv       <- y.names$y.groups=="Invertebrates"
  names.inverts <- unique(y.names$y.taxa[ind.inv])
  
  ind.alg     <- y.names$y.groups=="Algae"
  names.algae <- unique(y.names$y.taxa[ind.alg])
  
    
  # construct vector with masses of inverts
  
  masses <- streambugs:::get.inpind.parval.taxaprop(par.names = c("M"),        
                                       y.names   = y.names$y.names[ind.inv],
                                       par       = par.invtraits)$parvals
  
  rownames(masses) <- y.names$y.taxa[ind.inv]
  
  mass <- rep(NA,length(names.inverts))
  names(mass) <- names.inverts
 
  for (i in 1:length(names.inverts))
  {
    invert <- names.inverts[i]
    ind <- which(invert==rownames(masses))
    m.i <- unique(masses[ind,"M"])
    if(length(m.i)==1)
    {
      mass[i] <- m.i
    }  else { stop("mass for invert ",invert," not unique")}
  }
  
  if( sum(is.na(mass)) > 1 ) warning("masses for inverts ",names(mass)[is.na(mass)]," missing")
    
  
  filt  <- "aff|pff"
  shred <- "shr|min|xyl" 
  scra  <- "gra"
  coll  <- "gat"
  pred  <- "pre"
  para  <- "par"
  other <- "oth"
  
  # inverts <- Invertebrates
  invpars <- list()
  ind <- 0
  par.stoich.cons <- numeric(0)
  ind2 <- 0
  
  # das geht wohl noch nicht für zeitabh. Inputs
  par.feedtype <- streambugs:::get.inpind.parval.taxaprop.traits(trait.names="feedingtype",
                                                    y.names=y.names,
                                                    par=par.invtraits)$feedingtype$parvals
  
  # par.feedtype <- get.inpind.parval.taxaprop.traits(trait.names="feedingtype",
  #                                                                y.names=y.names,
  #                                                                par=par.invtraits)$parvals
  
  
  for (i in 1:length(names.inverts) )
  {
    invert <- names.inverts[i]
    
    indft  <- grep(invert,rownames(par.feedtype))
    parfts <- par.feedtype[indft,]  
    if(!is.null(nrow(parfts)))
    { parft  <- colMeans(parfts)} else #evtl. noch Fehlermeldung einbauen?
    { parft <- parfts }                 
    
    #indft <- grep(paste(invert,"feedingtype",sep="_"),names(par.invtraits))
    #parft <- par.invtraits[indft]
    
    parft <- parft[parft>0]
    
    ftypes <- unlist(strsplit(names(parft),split="_"))
    indft <- grep("feedingtype",ftypes)
    ftypes <- ftypes[indft+1]
    
    if ( sum(!is.na(ftypes))>0 )
    {
      #grazer:
      if ( sum(grepl(scra,ftypes)) > 0 )  
      {
        
        for (j in 1:length(names.algae))
        {
          alga <- names.algae[j]
          
          ind <- ind+1
          invpars[ind] <- Pref.default   
          names(invpars)[ind] <- paste(invert,alga,"Pref",sep="_")  
          
          ind2 <- ind2+1
          par.stoich.cons[ind2] <-  par.stoich.out["Cons_Invertebrates_Scra_Scra"]  
          names(par.stoich.cons)[ind2] <- paste("Cons",invert,alga,invert,sep="_")
          
          ind2 <- ind2+1
          par.stoich.cons[ind2] <- par.stoich.out["Cons_Invertebrates_Scra_Algae"]   
          names(par.stoich.cons)[ind2] <- paste("Cons",invert,alga,alga,sep="_")
          
          ind2 <- ind2+1
          par.stoich.cons[ind2] <- par.stoich.out["Cons_Invertebrates_Scra_FPOM"]   
          names(par.stoich.cons)[ind2] <- paste("Cons",invert,alga,"FPOM",sep="_")
        }
        
      }
        
      #shredder:
      if ( sum(grepl(shred,ftypes)) > 0 )  
      {
        ind <- ind+1
        invpars[ind] <- Pref.default   
        names(invpars)[ind] <- paste(invert,"CPOM_Pref",sep="")  
        
        ind2 <- ind2+1
        par.stoich.cons[ind2] <-  par.stoich.out["Cons_Invertebrates_Shred_Shred"]  
        names(par.stoich.cons)[ind2] <- paste("Cons",invert,"CPOM",invert,sep="_")
        
        ind2 <- ind2+1
        par.stoich.cons[ind2] <- par.stoich.out["Cons_Invertebrates_Shred_CPOM"]   
        names(par.stoich.cons)[ind2] <- paste("Cons",invert,"CPOM","CPOM",sep="_")
        
        ind2 <- ind2+1
        par.stoich.cons[ind2] <- par.stoich.out["Cons_Invertebrates_Shred_FPOM"]   
        names(par.stoich.cons)[ind2] <- paste("Cons",invert,"CPOM","FPOM",sep="_")
      }
      
      #filterer:
      if ( sum(grepl(filt,ftypes)) > 0 )  
      {
        ind <- ind+1
        invpars[ind] <- Pref.default   
        names(invpars)[ind] <- paste(invert,"SusPOM_Pref",sep="")  
        
        ind2 <- ind2+1
        par.stoich.cons[ind2] <-  par.stoich.out["Cons_Invertebrates_Filt_Filt"]  
        names(par.stoich.cons)[ind2] <- paste("Cons",invert,"SusPOM",invert,sep="_")
        
        #ind2 <- ind2+1
        #par.stoich.cons[ind2] <- par.stoich.out["Cons_Invertebrates_Filt_SusPOM"]   
        #names(par.stoich.cons)[ind2] <- paste("Cons",invert,"SusPOM","SusPOM",sep="_")
        
        ind2 <- ind2+1
        par.stoich.cons[ind2] <- par.stoich.out["Cons_Invertebrates_Filt_FPOM"]   
        names(par.stoich.cons)[ind2] <- paste("Cons",invert,"SusPOM","FPOM",sep="_")
      }
      
      # collector-gatherer
      if ( sum(grepl(coll,ftypes)) > 0 )  
      {
        ind <- ind+1
        invpars[ind] <- Pref.default   
        names(invpars)[ind] <- paste(invert,"FPOM_Pref",sep="")  
        
        ind2 <- ind2+1
        par.stoich.cons[ind2] <-  par.stoich.out["Cons_Invertebrates_Coll_Coll"]  
        names(par.stoich.cons)[ind2] <- paste("Cons",invert,"FPOM",invert,sep="_")
        
        ind2 <- ind2+1
        par.stoich.cons[ind2] <- par.stoich.out["Cons_Invertebrates_Coll_FPOM"]   
        names(par.stoich.cons)[ind2] <- paste("Cons",invert,"FPOM","FPOM",sep="_")
      }
      
      # predators
      if ( sum(grepl(pred,ftypes)) > 0 )  
      {
        
        prey <- names.inverts[-i]
        mass_prey <- mass[prey]
        
        for ( j in 1:length(prey))
        {
          if (mass[i] > mass_prey[j])
          {
            
            ind <- ind+1
            invpars[ind] <- Pref.default   
            names(invpars)[ind] <- paste(invert,prey[j],"Pref",sep="_")  
            
            ind2 <- ind2+1
            par.stoich.cons[ind2] <-  par.stoich.out["Cons_Invertebrates_Pred_Pred"]  
            names(par.stoich.cons)[ind2] <- paste("Cons",invert,prey[j],invert,sep="_")
            
            ind2 <- ind2+1
            par.stoich.cons[ind2] <- par.stoich.out["Cons_Invertebrates_Pred_Prey"]   
            names(par.stoich.cons)[ind2] <- paste("Cons",invert,prey[j],prey[j],sep="_")
            
            ind2 <- ind2+1
            par.stoich.cons[ind2] <- par.stoich.out["Cons_Invertebrates_Pred_FPOM"]   
            names(par.stoich.cons)[ind2] <- paste("Cons",invert,prey[j],"FPOM",sep="_")
            
          }
        }
      }
      
      if ( sum(grepl(para,ftypes)) > 0 )    
      { 
        warning(invert," has feedingtype parasite, not yet implemented",sep="") 
      } 
      
      if ( sum(grepl(other,ftypes)) > 0 )
      {
        warning(invert," has other feedingtype",sep="")           
      }
    }
    
    ind2 <- ind2+1
    par.stoich.cons[ind2] <- par.stoich.out["Death_Invertebrates_Invertebrates"]   
    names(par.stoich.cons)[ind2] <- paste("Death",invert,invert,sep="_")
   
    ind2 <- ind2+1
    par.stoich.cons[ind2] <- par.stoich.out["Death_Invertebrates_FPOM"]   
    names(par.stoich.cons)[ind2] <- paste("Death",invert,"FPOM",sep="_")   
    
  }  # loop over inverts
  
  for (k in 1:length(names.algae))
  {
    alga <- names.algae[k]
    
    ind2 <- ind2+1
    par.stoich.cons[ind2] <- par.stoich.out["Death_Algae_Algae"]   
    names(par.stoich.cons)[ind2] <- paste("Death",alga,alga,sep="_")
    
    ind2 <- ind2+1
    par.stoich.cons[ind2] <- par.stoich.out["Death_Invertebrates_FPOM"]   
    names(par.stoich.cons)[ind2] <- paste("Death",alga,"FPOM",sep="_")
  }
  
  return(par.stoich.cons) 
  
}

#test <- assign.par.stoich.consumption(par.invtraits,par.stoich.out,y.names)