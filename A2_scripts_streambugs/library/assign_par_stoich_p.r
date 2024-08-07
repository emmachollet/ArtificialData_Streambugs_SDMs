#adapted to use streambugs package

assign.par.stoich <- function(par.invtraits,par.stoich.out,y.names,
                              Pref.default=1,pref.use=FALSE,CPOM.dif=FALSE,ratio_pred_prey=1)

  # pref.use: if TRUE, affinity scores for foodsources are used as preference parameters,
  #           if TRUE but no foodsources are provided a warning is issued, 
  #           preference parameters from feedingtypes are not yet implemented
  # CPOM.dif: if TRUE, two types of CPOM (CPOMa and CPOMp) are introduces with stoichiometry of animals and plants
  
  {
##################################################################################################
######## NEW FORMULATION, removed and value by default =1, if no values assigned before in par.ini
##################################################################################################

 # if(length(grep("ratio_pred_prey", names(par))) == 0) stop("The parameter vector must contain the minimum ratio between predator and prey biomass")
#  if("ratio_pred_prey_transformed" %in% names(par))
#{
#  ratio_pred_prey = exp(par["ratio_pred_prey_transformed"])
#  names(ratio_pred_prey) = "ratio_pred_prey"
#} else {
#  ratio_pred_prey = par["ratio_pred_prey"]
#}
# 

# if("ratio_pred_prey_transformed" %in% names(par.ini))
#   {
#     ratio_pred_prey = exp(par.ini[21])
#     names(ratio_pred_prey) = "ratio_pred_prey"
#   } else {
#     ratio_pred_prey = par["ratio_pred_prey"]
#   }

########################################################## 
    
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
  
  # Feeding types
  
  filt  <- c("filter.feeder","aff", "pff")
  shred <- c("shredder","shr", "min", "xyl")           
  scra  <- c("scraper", "gra")
  coll  <- c("deposit.feeder", "gat")
  pred  <- c("predator", "pre")
  pierc <- "piercer"
  para  <- c("parasite", "par")
  abs <- "absorber"
  other <- "oth"
  
  # Food sources
  
  mic    <- "microorganisms"
  det    <- "detritus"
  dpl    <- "dead.plant"
  micpl  <- "living.microphytes"
  macpl  <- "living.macrophytes"
  dinv   <- "dead.animal"
  micinv <- "living.microinvertebrates"
  macinv <- "living.macroinvertebrates"
  vert   <- "vertebrates"
  
  # inverts <- Invertebrates
  invpars <- list()
  ind <- 0
  par.stoich.cons <- numeric(0)
  ind2 <- 0
  
  # das geht wohl noch nicht für zeitabh. Inputs
  par.feedtype <-  streambugs:::get.inpind.parval.taxaprop.traits(trait.names="feedingtype",
                                                    y.names=y.names,
                                                    par=par.invtraits,inp=inp)$feedingtype$parvals
  if(length(grep("food", names(par.invtraits))) > 0)
  {
    par.foodsources <-  streambugs:::get.inpind.parval.taxaprop.traits(trait.names="food",
                                                         y.names=y.names,
                                                         par=par.invtraits)$food$parvals
  } else {
    par.foodsources <- NA
  }
  
  rownames.par.feedtype <- strsplit(rownames(par.feedtype),split="_")
  
  for (i in 1:length(names.inverts) )
  {
    invert <- names.inverts[i]
    
    # feeding types
    
    indft <- NULL
    
    for(j in 1:length(rownames(par.feedtype)))
    {
      if(invert==rownames.par.feedtype[[j]][3]) indft <- c(indft,j)      
    }
    
    if(!is.null(nrow(par.feedtype[indft,])))
    {
      uniqueparft <- apply(par.feedtype[indft,],2,unique)
      if(length(unlist(uniqueparft)) > ncol(par.feedtype)) 
      {
        warning(paste("feedingtypes are not unique for invert",
                      invert,"in different habitats and/or reaches"))
        
      }
    }
    
    parft <- par.feedtype[indft[1],] 
        
    parft <- parft[parft>0]
    
    ftypes <- unlist(strsplit(names(parft),split="_"))
    indft <- grep("feedingtype",ftypes)
    ftypes <- ftypes[indft+1]
    
    # food sources
    
    if(length(na.omit(par.foodsources)) > 0)
    {
      
      rownames.par.foodsources <- strsplit(rownames(par.foodsources),split="_")
      
      indfs <- NULL
      for(j in 1:length(rownames(par.foodsources)))
      {
        if(invert==rownames.par.foodsources[[j]][3]) indfs <- c(indfs,j)      
      }
      
      if(!is.null(nrow(par.foodsources[indfs,])))
      {
        uniqueparfs <- apply(par.foodsources[indfs,],2,unique)
        if(length(unlist(uniqueparfs)) > ncol(par.foodsources)) 
        {
          warning(paste("foodsources are not unique for invert",
                        invert,"in different habitats and/or reaches"))
          
        }
      }
      
      parfs <- par.foodsources[indfs[1],]
      
      parfs <- parfs[parfs>0]
      
      fsources <- unlist(strsplit(names(parfs), split="_"))
      indfs <- grep("food", fsources)
      fsources <- fsources[indfs+1]
      
    }  else {
      
      if(pref.use) warning("pref.use is true but no foodsources are provided")
    }  
    
    if ( sum(!is.na(ftypes))>0 )
    {
      #grazer:
      if ( sum(which(scra %in% ftypes)) > 0 )  
      {
        if(length(na.omit(par.foodsources)) <= 0)     # general case, no information on food sources (consider that they feed on all types of algae)
        {
          for (j in 1:length(names.algae))
          {
            alga <- names.algae[j]
            
            ind <- ind+1
            invpars[ind] <- Pref.default   
            names(invpars)[ind] <- paste(invert,alga,"Pref",sep="_")  
            
            temp           <- numeric(0)
            temp[1]        <- par.stoich.out["Cons_Invertebrates_Scra_Scra"]
            names(temp)[1] <- paste("Cons",invert,alga,invert,sep="_")
            temp[2]        <- par.stoich.out["Cons_Invertebrates_Scra_Algae"]
            names(temp)[2] <- paste("Cons",invert,alga,alga,sep="_")
            temp[3]        <- par.stoich.out["Cons_Invertebrates_Scra_FPOM"] 
            names(temp)[3] <- paste("Cons",invert,alga,"FPOM",sep="_")
            
            if(names(temp)[2] %in% names(par.stoich.cons)){
              if(abs(temp[2]) > abs(par.stoich.cons[names(par.stoich.cons) %in% names(temp)[2]]))
              {
                par.stoich.cons[match(names(temp), names(par.stoich.cons))] <- temp
              }
            } else {
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[1]
              names(par.stoich.cons)[ind2] <- names(temp)[1]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[2] 
              names(par.stoich.cons)[ind2] <- names(temp)[2]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[3]  
              names(par.stoich.cons)[ind2] <- names(temp)[3]
              
            }
            
          }
        } else {                                 
          # Specific case where
          # information on food sources
          # is available
          
          if ( sum(which(micpl %in% fsources)) > 0 )  # scrapers feeding on biofilm
          {
            
            alga <- "crustyAlgae"
            
            ind <- ind+1
            if(pref.use == FALSE)
            {
              invpars[ind] <- Pref.default
            } else {
              invpars[ind] <- parfs[which(names(parfs)==paste("food_",micpl,sep=""))]
            }
               
            names(invpars)[ind] <- paste(invert,alga,"Pref",sep="_")  
            
            temp           <- numeric(0)
            temp[1]        <- par.stoich.out["Cons_Invertebrates_Scra_Scra"]
            names(temp)[1] <- paste("Cons",invert,alga,invert,sep="_")
            temp[2]        <- par.stoich.out["Cons_Invertebrates_Scra_Algae"]
            names(temp)[2] <- paste("Cons",invert,alga,alga,sep="_")
            temp[3]        <- par.stoich.out["Cons_Invertebrates_Scra_FPOM"] 
            names(temp)[3] <- paste("Cons",invert,alga,"FPOM",sep="_")
            
            if(names(temp)[2] %in% names(par.stoich.cons)){
              if(abs(temp[2]) > abs(par.stoich.cons[names(par.stoich.cons) %in% names(temp)[2]]))
              {
                par.stoich.cons[match(names(temp), names(par.stoich.cons))] <- temp
              }
            } else {
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[1]
              names(par.stoich.cons)[ind2] <- names(temp)[1]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[2]  
              names(par.stoich.cons)[ind2] <- names(temp)[2]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[3]  
              names(par.stoich.cons)[ind2] <- names(temp)[3]
              
            }
          }
          
          if ( sum(which(macpl %in% fsources)) > 0 )  # scrapers feeding on living
          {                                      # macrophytes
            
            alga <- "filamentousAlgae"           # for the moment the only
            # 'macrophytes' are filamentous
            # algae
            
            ind <- ind+1
            if(pref.use == FALSE)
            {
              invpars[ind] <- Pref.default
            } else {
              invpars[ind] <- parfs[which(names(parfs)==paste("food_",macpl,sep=""))]   
            }
            names(invpars)[ind] <- paste(invert,alga,"Pref",sep="_")  
            
            temp           <- numeric(0)
            temp[1]        <- par.stoich.out["Cons_Invertebrates_Scra_Scra"]
            names(temp)[1] <- paste("Cons",invert,alga,invert,sep="_")
            temp[2]        <- par.stoich.out["Cons_Invertebrates_Scra_Algae"]
            names(temp)[2] <- paste("Cons",invert,alga,alga,sep="_")
            temp[3]        <- par.stoich.out["Cons_Invertebrates_Scra_FPOM"] 
            names(temp)[3] <- paste("Cons",invert,alga,"FPOM",sep="_")
            
            if(names(temp)[2] %in% names(par.stoich.cons)){
              if(abs(temp[2]) > abs(par.stoich.cons[names(par.stoich.cons) %in% names(temp)[2]]))
              {
                par.stoich.cons[match(names(temp), names(par.stoich.cons))] <- temp
              }
            } else {
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[1]   
              names(par.stoich.cons)[ind2] <- names(temp)[1]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[2] 
              names(par.stoich.cons)[ind2] <- names(temp)[2]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[3]  
              names(par.stoich.cons)[ind2] <- names(temp)[3]
              
            }
          }
          
        }
        
      }
      
      #shredder:
      if ( sum(which(shred %in% ftypes)) > 0 )  
      {
        if(length(na.omit(par.foodsources)) <= 0)     # general case, no information on food sources (consider that they feed on not determined CPOM)
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
        } else {                                 # Specific case where
          # information on food sources
          # is available
          
          if ( sum(which(macpl %in% fsources)) > 0 ) # shredders feeding on living
          {                                    # macrophytes
            food   <- "filamentousAlgae"       # for the moment the only
            shred2 <- "Shredm"                 # 'macrophytes' are filamentous
            # algae
            
            ind <- ind+1
            if(pref.use == FALSE)
            {
              invpars[ind] <- Pref.default
            } else {
              invpars[ind] <- parfs[which(names(parfs)==paste("food_",macpl,sep=""))]
            }
            names(invpars)[ind] <- paste(invert,food,"Pref",sep="_")  
            
            temp           <- numeric(0)
            temp[1]        <- par.stoich.out[paste("Cons_Invertebrates", shred2, shred2, sep = "_")]  
            names(temp)[1] <- paste("Cons",invert,food,invert,sep="_")
            temp[2]        <- par.stoich.out[paste("Cons_Invertebrates", shred2, food, sep="_")]   
            names(temp)[2] <- paste("Cons",invert,food,food,sep="_")
            temp[3]        <- par.stoich.out[paste("Cons_Invertebrates", shred2, "FPOM", sep="_")] 
            names(temp)[3] <- paste("Cons",invert,food,"FPOM",sep="_")
            
            if(names(temp)[2] %in% names(par.stoich.cons)){
              if(abs(temp[2]) > abs(par.stoich.cons[names(par.stoich.cons) %in% names(temp)[2]]))
              {
                par.stoich.cons[match(names(temp), names(par.stoich.cons))] <- temp
              }
            } else {
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[1]  
              names(par.stoich.cons)[ind2] <- names(temp)[1]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[2]  
              names(par.stoich.cons)[ind2] <- names(temp)[2]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[3]  
              names(par.stoich.cons)[ind2] <- names(temp)[3]
              
            }
          }
          
          if ( sum(which(dpl %in% fsources)) > 0 ) # shredders feeding on dead
          {                                   # plants
            food <- "CPOMp"                     
            shred2 <- "Shredp"                   
            
            ind <- ind+1
            if(pref.use == FALSE)
            {
              invpars[ind] <- Pref.default
            } else {
              invpars[ind] <- parfs[which(names(parfs)==paste("food_",dpl,sep=""))]
            }
            names(invpars)[ind] <- paste(invert,food,"Pref",sep="_")  
            
            temp           <- numeric(0)
            temp[1]        <- par.stoich.out[paste("Cons_Invertebrates", shred2, shred2, sep = "_")]  
            names(temp)[1] <- paste("Cons",invert,food,invert,sep="_")
            temp[2]        <- par.stoich.out[paste("Cons_Invertebrates", shred2, food, sep="_")]   
            names(temp)[2] <- paste("Cons",invert,food,food,sep="_")
            temp[3]        <- par.stoich.out[paste("Cons_Invertebrates", shred2, "FPOM", sep="_")]   
            names(temp)[3] <- paste("Cons",invert,food,"FPOM",sep="_")
            
            if(names(temp)[2] %in% names(par.stoich.cons)){
              if(abs(temp[2]) > abs(par.stoich.cons[names(par.stoich.cons) %in% names(temp)[2]]))
              {
                par.stoich.cons[match(names(temp), names(par.stoich.cons))] <- temp
              }
            } else {
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[1]    
              names(par.stoich.cons)[ind2] <- names(temp)[1]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[2]  
              names(par.stoich.cons)[ind2] <- names(temp)[2]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[3]  
              names(par.stoich.cons)[ind2] <- names(temp)[3]
              
            }
          }
          
          if ( sum(which(dinv %in% fsources)) > 0 ) # shredders feeding on dead
          {                                   # animals
            food <- "CPOMa"                     
            shred2 <- "Shreda"                   
            
            ind <- ind+1
            if(pref.use == FALSE)
            {
              invpars[ind] <- Pref.default
            } else {
              invpars[ind] <- parfs[which(names(parfs)==paste("food_",dinv,sep=""))]
            }
            names(invpars)[ind] <- paste(invert,food,"Pref",sep="_")  
            
            temp           <- numeric(0)
            temp[1]        <- par.stoich.out[paste("Cons_Invertebrates", shred2, shred2, sep = "_")]  
            names(temp)[1] <- paste("Cons",invert,food,invert,sep="_")
            temp[2]        <- par.stoich.out[paste("Cons_Invertebrates", shred2, food, sep="_")]   
            names(temp)[2] <- paste("Cons",invert,food,food,sep="_")
            temp[3]        <- par.stoich.out[paste("Cons_Invertebrates", shred2, "FPOM", sep="_")]   
            names(temp)[3] <- paste("Cons",invert,food,"FPOM",sep="_")
            
            if(names(temp)[2] %in% names(par.stoich.cons)){
              if(abs(temp[2]) > abs(par.stoich.cons[names(par.stoich.cons) %in% names(temp)[2]]))
              {
                par.stoich.cons[match(names(temp), names(par.stoich.cons))] <- temp
              }
            } else {
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[1]    
              names(par.stoich.cons)[ind2] <- names(temp)[1]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[2]  
              names(par.stoich.cons)[ind2] <- names(temp)[2]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[3]
              names(par.stoich.cons)[ind2] <- names(temp)[3]
              
            }
          }
          
        }
      }
      
      #filterer:
      if ( sum(which(filt %in% ftypes)) > 0 )  
      {
        if(length(na.omit(par.foodsources)) <= 0)
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
          
        } else {          
          if ( sum(which(det %in% fsources)) > 0 )  # filterer feeding on FPOM
          {
            ind <- ind+1
            if(pref.use == FALSE)
            {
              invpars[ind] <- Pref.default
            } else {
              invpars[ind] <- parfs[which(names(parfs)==paste("food_",det,sep=""))]
            }
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
        }
      }
      
      # collector-gatherer
      if ( sum(which(coll %in% ftypes)) > 0 )  
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
      if ( sum(which(pred %in% ftypes)) > 0 )  
      {
        
        prey <- names.inverts[-i]
        mass_prey <- mass[prey]
        
        for ( j in 1:length(prey))
        {
          if (mass[i] > ratio_pred_prey*mass_prey[j])
          {
            
            ind <- ind+1
            invpars[ind] <- Pref.default   
            names(invpars)[ind] <- paste(invert,prey[j],"Pref",sep="_")  
            
            temp           <- numeric(0)
            temp[1]        <- par.stoich.out["Cons_Invertebrates_Pred_Pred"]  
            names(temp)[1] <- paste("Cons",invert,prey[j],invert,sep="_")
            temp[2]        <- par.stoich.out["Cons_Invertebrates_Pred_Prey"]   
            names(temp)[2] <- paste("Cons",invert,prey[j],prey[j],sep="_")
            temp[3]        <- par.stoich.out["Cons_Invertebrates_Pred_FPOM"]   
            names(temp)[3] <- paste("Cons",invert,prey[j],"FPOM",sep="_")
            
            if(names(temp)[2] %in% names(par.stoich.cons)){
              if(abs(temp[2]) > abs(par.stoich.cons[names(par.stoich.cons) %in% names(temp)[2]]))
              {
                par.stoich.cons[match(names(temp), names(par.stoich.cons))] <- temp
              }
            } else {
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[1] 
              names(par.stoich.cons)[ind2] <- names(temp)[1]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[2]  
              names(par.stoich.cons)[ind2] <- names(temp)[2]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[3]  
              names(par.stoich.cons)[ind2] <- names(temp)[3]
              
            }
          }
        }
      }
      
      #piercers:
      if ( sum(which(pierc %in% ftypes)) > 0 )  
      {
        if(length(na.omit(par.foodsources)) <= 0)     # general case, no information on food sources (consider that they feed on both macrophytes and other macroinvertebrate)
        {
          prey <- names.inverts[-i]
          pierc2 <- "Pierca"
          
          for ( j in 1:length(prey))
          {
            ind <- ind+1
            invpars[ind] <- Pref.default   
            names(invpars)[ind] <- paste(invert,prey[j],"Pref",sep="_")  
            
            temp           <- numeric(0)
            temp[1]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, pierc2, sep="_")]  
            names(temp)[1] <- paste("Cons",invert,prey[j],invert,sep="_")
            temp[2]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, "Prey", sep="_")]   
            names(temp)[2] <- paste("Cons",invert,prey[j],prey[j],sep="_")
            temp[3]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, "FPOM", sep="_")]   
            names(temp)[3] <- paste("Cons",invert,prey[j],"FPOM",sep="_")
            
            if(names(temp)[2] %in% names(par.stoich.cons)){
              if(abs(temp[2]) > abs(par.stoich.cons[names(par.stoich.cons) %in% names(temp)[2]]))
              {
                par.stoich.cons[match(names(temp), names(par.stoich.cons))] <- temp
              }
            } else {
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[1]   
              names(par.stoich.cons)[ind2] <- names(temp)[1]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[2]  
              names(par.stoich.cons)[ind2] <- names(temp)[2]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[3]  
              names(par.stoich.cons)[ind2] <- names(temp)[3]
              
            }
          }
          
          food <- "filamentousAlgae"          # for the moment the only
          pierc2 <- "Piercm"                  # 'macrophytes' are filamentous
          # algae
          
          ind <- ind+1
          invpars[ind] <- Pref.default   
          names(invpars)[ind] <- paste(invert,food,"Pref",sep="_")  
          
          temp           <- numeric(0)
          temp[1]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, pierc2, sep = "_")]    
          names(temp)[1] <- paste("Cons",invert,food,invert,sep="_")
          temp[2]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, food, sep="_")]      
          names(temp)[2] <- paste("Cons",invert,food,food,sep="_")
          temp[3]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, "FPOM", sep="_")]     
          names(temp)[3] <- paste("Cons",invert,food,"FPOM",sep="_")
          
          if(names(temp)[2] %in% names(par.stoich.cons)){
            if(abs(temp[2]) > abs(par.stoich.cons[names(par.stoich.cons) %in% names(temp)[2]]))
            {
              par.stoich.cons[match(names(temp), names(par.stoich.cons))] <- temp
            }
          } else {
            ind2 <- ind2+1
            par.stoich.cons[ind2]        <-  temp[1]    
            names(par.stoich.cons)[ind2] <- names(temp)[1]
            
            ind2 <- ind2+1
            par.stoich.cons[ind2]        <-  temp[2]  
            names(par.stoich.cons)[ind2] <- names(temp)[2]
            
            ind2 <- ind2+1
            par.stoich.cons[ind2]        <-  temp[3]  
            names(par.stoich.cons)[ind2] <- names(temp)[3]
            
          }
          
        } else {
          if ( sum(which(macpl %in% fsources)) > 0 ) # piercers feeding on living
          {                                     # macrophytes
            food <- "filamentousAlgae"          # for the moment the only
            pierc2 <- "Piercm"                  # 'macrophytes' are filamentous
            # algae
            
            ind <- ind+1
            if(pref.use == FALSE)
            {
              invpars[ind] <- Pref.default
            } else {
              invpars[ind] <- parfs[which(names(parfs)==paste("food_",macpl,sep=""))]
            }
            names(invpars)[ind] <- paste(invert,food,"Pref",sep="_")  
            
            temp           <- numeric(0)
            temp[1]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, pierc2, sep = "_")]    
            names(temp)[1] <- paste("Cons",invert,food,invert,sep="_")
            temp[2]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, food, sep="_")]      
            names(temp)[2] <- paste("Cons",invert,food,food,sep="_")
            temp[3]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, "FPOM", sep="_")]     
            names(temp)[3] <- paste("Cons",invert,food,"FPOM",sep="_")
            
            if(names(temp)[2] %in% names(par.stoich.cons)){
              if(abs(temp[2]) > abs(par.stoich.cons[names(par.stoich.cons) %in% names(temp)[2]]))
              {
                par.stoich.cons[match(names(temp), names(par.stoich.cons))] <- temp
              }
            } else {
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[1]    
              names(par.stoich.cons)[ind2] <- names(temp)[1]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[2]  
              names(par.stoich.cons)[ind2] <- names(temp)[2]
              
              ind2 <- ind2+1
              par.stoich.cons[ind2]        <-  temp[3]  
              names(par.stoich.cons)[ind2] <- names(temp)[3]
              
            }
          }
          
          if ( sum(which(macinv %in% fsources)) > 0 ) # piercers feeding on living
          {                                      # macroinvertebrates
            prey <- names.inverts[-i]
            pierc2 <- "Pierca"
            
            for ( j in 1:length(prey))
            {
              ind <- ind+1
              if(pref.use == FALSE)
              {
                invpars[ind] <- Pref.default
              } else {
                invpars[ind] <- parfs[which(names(parfs)==paste("food_",macinv,sep=""))]
              }
              names(invpars)[ind] <- paste(invert,prey[j],"Pref",sep="_")  
              
              temp           <- numeric(0)
              temp[1]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, pierc2, sep="_")]  
              names(temp)[1] <- paste("Cons",invert,prey[j],invert,sep="_")
              temp[2]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, "Prey", sep="_")]   
              names(temp)[2] <- paste("Cons",invert,prey[j],prey[j],sep="_")
              temp[3]        <- par.stoich.out[paste("Cons_Invertebrates", pierc2, "FPOM", sep="_")]   
              names(temp)[3] <- paste("Cons",invert,prey[j],"FPOM",sep="_")
              
              if(names(temp)[2] %in% names(par.stoich.cons)){
                if(abs(temp[2]) > abs(par.stoich.cons[names(par.stoich.cons) %in% names(temp)[2]]))
                {
                  par.stoich.cons[match(names(temp), names(par.stoich.cons))] <- temp
                }
              } else {
                ind2 <- ind2+1
                par.stoich.cons[ind2]        <-  temp[1]   
                names(par.stoich.cons)[ind2] <- names(temp)[1]
                
                ind2 <- ind2+1
                par.stoich.cons[ind2]        <-  temp[2]  
                names(par.stoich.cons)[ind2] <- names(temp)[2]
                
                ind2 <- ind2+1
                par.stoich.cons[ind2]        <-  temp[3]  
                names(par.stoich.cons)[ind2] <- names(temp)[3]
                
              }
            }
          }
        }
      }
      
      
      if ( sum(which(para %in% ftypes)) > 0 )    
      { 
        warning(invert," has feeding type parasite, not yet implemented",sep="") 
      } 
      
      if ( sum(which(abs %in% ftypes)) > 0 )    
      { 
        warning(invert," has feeding type absorber, not yet implemented",sep="") 
      } 
      
      if ( sum(which(other %in% ftypes)) > 0 )
      {
        warning(invert," has other feedingtype",sep="")           
      }
      
      if(length(na.omit(par.foodsources)) > 0)
      {
        if ( sum(which(mic %in% fsources)) > 0 )
        {
          warning(invert," feed on microorganisms, not yet implemented",sep="")           
        }
        
        if ( sum(which(micinv %in% fsources)) > 0 )
        {
          warning(invert," feed on microinvertebrates, not yet implemented",sep="")           
        }
        
        if ( sum(which(vert %in% fsources)) > 0 )
        {
          warning(invert," feed on vertebrates, not yet implemented",sep="")           
        }
        
      }      
    }
    
    ind2 <- ind2+1
    if(CPOM.dif == FALSE)
    {
      par.stoich.cons[ind2] <- par.stoich.out["Death1_Invertebrates_Invertebrates"] 
    } else {
      par.stoich.cons[ind2] <- par.stoich.out["Death2_Invertebrates_Invertebrates"] 
    }
    
    names(par.stoich.cons)[ind2] <- paste("Death",invert,invert,sep="_")
    
    ind2 <- ind2+1
    if(CPOM.dif == FALSE)
    {
      par.stoich.cons[ind2] <- par.stoich.out["Death1_Invertebrates_FPOM"]  
      names(par.stoich.cons)[ind2] <- paste("Death",invert,"FPOM",sep="_")  
    } else {
      par.stoich.cons[ind2] <- par.stoich.out["Death2_Invertebrates_CPOMa"]  
      names(par.stoich.cons)[ind2] <- paste("Death",invert,"CPOMa",sep="_")  
    }
    
    
    
  }  # loop over inverts
  
  for (k in 1:length(names.algae))
  {
    alga <- names.algae[k]
    
    ind2 <- ind2+1
    par.stoich.cons[ind2] <- par.stoich.out["Death_Algae_Algae"]   
    names(par.stoich.cons)[ind2] <- paste("Death",alga,alga,sep="_")
    
    ind2 <- ind2+1
    par.stoich.cons[ind2] <- par.stoich.out["Death_Algae_FPOM"]   
    names(par.stoich.cons)[ind2] <- paste("Death",alga,"FPOM",sep="_")
  }
  
  
  return(par.stoich.cons) 
  
}

#test <- assign.par.stoich.consumption(par.invtraits,par.stoich.out,y.names)