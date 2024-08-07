
construct.invpars.traits <- function(Invertebrates,
                                     file.db.taxonomy="databases/db_taxonomy.dat",
                                     file.db.bodymass="databases/db_bodymass_meas_aquaplus.dat",
                                     file.db.microhab=NA,
                                     file.db.current=NA,
                                     file.db.temp=NA,
                                     file.db.sapro=NA,
                                     file.db.pH=NA,
                                     file.db.spear=NA,
                                     file.db.food=NA,
                                     file.db.feeding="databases/db_freshwaterecology_20131206.dat",
                                     method.bodymass = mean,
                                     method.microhab = mean,
                                     method.current = mean,
                                     method.temp = mean,
                                     method.sapro = mean,
                                     method.pH = mean,
                                     method.spear = min,
                                     method.food = mean,
                                     method.feeding = mean)
{
  
  ###############
  ## Functions ##
  ###############
  
  # Remove space in row names
  # -------------------------
  
  remove.space <- function (stringalt)
  {
    if(grepl(" ",stringalt))
    {
      splitstring <- unlist(strsplit(stringalt,split=" "))
      stringneu   <- paste(splitstring,collapse="")  
    } else
    {
      stringneu <- stringalt 
    }  
    return(stringneu)
  }
  
  remove.trait <- function (stringalt)
  {
    if(grepl("_",stringalt))
    {
      splitstring <- unlist(strsplit(stringalt,split="_"))
      stringneu   <- splitstring[[2]]  
    } else
    {
      stringneu <- stringalt 
    }  
    return(stringneu)
  }
  
  #######################
  ## Loading databases ##
  #######################
  
  db.taxonomy   = read.table(file.db.taxonomy,header=T,sep="\t")  
  db.taxonomy[,"Species2"] <- unlist(lapply(as.character(db.taxonomy[,"Species"]),remove.space))
  db.taxonomy <- apply(db.taxonomy, 2, as.character)
  
  db.bodymass   = read.table(file.db.bodymass,header=T,sep="\t")
  rownames(db.bodymass) <- lapply(as.character(db.bodymass[,"Taxon"]),remove.space)
    
  if(!is.na(file.db.microhab))
  {
    db.microhab           <- read.table(file.db.microhab,header=T,sep="\t")
    db.microhab           <- db.microhab[,c(grep("Taxon", colnames(db.microhab)), grep("microhab_", colnames(db.microhab)))]
    rownames(db.microhab) <- lapply(as.character(db.microhab[,"Taxon"]),remove.space)
    db.microhab           <- scale.traits(db.microhab)
    colnames(db.microhab) <- c("Taxon", paste("type", 1:(ncol(db.microhab)-1), sep=""))
#     colnames(db.microhab) <- lapply(colnames(db.microhab), remove.trait)
  }
  
  if(!is.na(file.db.current))
  {
    db.current            <- read.table(file.db.current,header=T,sep="\t")
    db.current            <- db.current[,c(grep("Taxon", colnames(db.current)), grep("current_", colnames(db.current)))]  
    rownames(db.current)  <- lapply(as.character(db.current[,"Taxon"]),remove.space)
    db.current            <- scale.traits(db.current)
    colnames(db.current)  <- c("Taxon", paste("class", 1:(ncol(db.current)-1), sep=""))
#     colnames(db.current)  <- lapply(colnames(db.current), remove.trait)
  }
  
  if(!is.na(file.db.temp))
  {
    db.temp               <- read.table(file.db.temp,header=T,sep="\t")
    db.temp               <- db.temp[,c(grep("Taxon", colnames(db.temp)), grep("temp_", colnames(db.temp)))]  
    rownames(db.temp)     <- lapply(as.character(db.temp[,"Taxon"]),remove.space)
    db.temp               <- scale.traits(db.temp)
    colnames(db.temp)     <- c("Taxon", paste("class", 1:(ncol(db.temp)-1), sep=""))    
#     colnames(db.temp)     <- lapply(colnames(db.temp), remove.trait)  
  }
  
  if(!is.na(file.db.sapro))
  {
    db.sapro              <- read.table(file.db.sapro,header=T,sep="\t")
    db.sapro              <- db.sapro[,c(grep("Taxon", colnames(db.sapro)), grep("sapro_", colnames(db.sapro)))]
    rownames(db.sapro)    <- lapply(as.character(db.sapro[,"Taxon"]),remove.space)
    db.sapro              <- scale.traits(db.sapro)
    colnames(db.sapro)    <- c("Taxon", paste("class", 0:(ncol(db.sapro)-2), sep=""))   # 0 should refer to xeno, 1 to oligo etc.
#     colnames(db.sapro)    <- lapply(colnames(db.sapro), remove.trait)
  }
  
  if(!is.na(file.db.pH))
  {
    db.pH            <- read.table(file.db.pH,header=T,sep="\t")
    db.pH            <- db.pH[,c(grep("Taxon", colnames(db.pH)), grep("pH_", colnames(db.pH)))] 
    rownames(db.pH)  <- lapply(as.character(db.pH[,"Taxon"]),remove.space)
    db.pH            <- scale.traits(db.pH)
    colnames(db.pH)  <- c("Taxon", paste("class", 1:(ncol(db.pH)-1), sep=""))
#     colnames(db.pH)       <- lapply(colnames(db.pH), remove.trait)
  }
  
  if(!is.na(file.db.spear))
  {
    db.spear      = read.table(file.db.spear,header=T,sep="\t")
    rownames(db.spear)    <- lapply(as.character(db.spear[,1]),remove.space)
  }
  
  if(!is.na(file.db.food))
  {
    db.food            <- read.table(file.db.food,header=T,sep="\t")
    db.food            <- db.food[,c(grep("Taxon", colnames(db.food)), grep("food_", colnames(db.food)))]
    rownames(db.food)  <- lapply(as.character(db.food[,"Taxon"]),remove.space)
    colnames(db.food)  <- lapply(colnames(db.food), remove.trait)
    db.food            <- scale.traits(db.food)
  }
  
  if(!is.na(file.db.feeding))
  {
    db.feeding            <- read.table(file.db.feeding,header=T,sep="\t")
    db.feeding            <- db.feeding[,c(grep("Taxon", colnames(db.feeding)), grep("feedingtype_", colnames(db.feeding)))]
    rownames(db.feeding)  <- lapply(as.character(db.feeding[,"Taxon"]),remove.space)
    colnames(db.feeding)  <- lapply(colnames(db.feeding), remove.trait)
    db.feeding            <- scale.traits(db.feeding)
  }
  
  
  par.traits <- NULL
  par.M      <- NULL
  par.orgmicropolltolval <- NULL
  
  ##############
  ## Bodymass ##
  ##############
  
  if(!is.na(file.db.bodymass))
  {
    par.M <- numeric(0)
    cind <- which(colnames(db.bodymass)=="mass_mg")
    par.M.default <- signif(median(db.bodymass[,cind],na.rm=TRUE),digits=3) #
    
    for (i in 1:length(Invertebrates))
    {
      parM <- construct.traitpars(invert=Invertebrates[i],
                                  cind=cind,
                                  method=method.bodymass,
                                  name.trait="M",
                                  default=par.M.default,
                                  db.traits=db.bodymass,
                                  db.taxonomy=db.taxonomy)/1000  #gDM!
      
      #delete "_mean" from the name
      nstop <- gregexpr("_M_",names(parM))[[1]][1]+1  
      names(parM) <- substring(names(parM),1,nstop)
      
      if(length(parM)>0) par.M <- c(par.M,parM)
    }
    
  }
  
  ################################
  ## Fitting to the environment ##
  ################################
  
  # Microhabitat:
  # -------------
  
  if(!is.na(file.db.microhab))
  {
    par.microhab <- numeric(0)
    
    for (i in 1:length(Invertebrates))
    {
      par <- construct.traitpars(invert=Invertebrates[i],
                                 cind=2:ncol(db.microhab),
                                 method=method.microhab,
                                 name.trait="microhabtolval",
                                 default=1,
                                 db.traits=db.microhab,
                                 db.taxonomy=db.taxonomy)
      
      if(length(par)>0) par.microhab <- c(par.microhab,par)
    }
    
    par.traits <- c(par.traits,par.microhab)
    
  }
  
  
  # Current tolerance:
  # ------------------
  
  if(!is.na(file.db.current))
  {
    par.currenttol <- numeric(0)
        
    for (i in 1:length(Invertebrates))
    {
      par <- construct.traitpars(Invertebrates[i],
                                 cind=2:ncol(db.current),
                                 method=method.current,
                                 name.trait="currenttolval",
                                 default=1, 
                                 db.traits=db.current,
                                 db.taxonomy=db.taxonomy)
      
      if(length(par)>0) par.currenttol <- c(par.currenttol,par)
    }
    
    par.traits <- c(par.traits,par.currenttol)
    
  }
  
  #################################
  ## Physico-chemical tolerances ##
  #################################
  
  # Temperature tolerance:
  # ----------------------
  
  if(!is.na(file.db.temp))
  {
    
    par.tempmaxtol <- numeric(0)
        
    for (i in 1:length(Invertebrates))
    {
      par <- construct.traitpars(Invertebrates[i],
                                 cind=2:ncol(db.temp),
                                 method=method.temp,
                                 name.trait="tempmaxtolval",
                                 default=1,
                                 db.traits=db.temp,
                                 db.taxonomy=db.taxonomy)
      
      if(length(par)>0) par.tempmaxtol <- c(par.tempmaxtol,par)
    }
    
    par.traits <- c(par.traits,par.tempmaxtol)
    
  }
  
  # Saprobity:
  # ----------
  if(!is.na(file.db.sapro))
  {
    par.saprotolval <- numeric(0)
    
     for (i in 1:length(Invertebrates))
    {
      par <- construct.traitpars(Invertebrates[i],
                                 cind=2:ncol(db.sapro),
                                 method=method.sapro,
                                 name.trait="saprotolval",
                                 default=1,
                                 db.traits=db.sapro,
                                 db.taxonomy=db.taxonomy)
      
      if(length(par)>0) par.saprotolval <- c(par.saprotolval,par)
    }
    
    par.traits <- c(par.traits,par.saprotolval)
    
  }
  # pH
  # --
  
  if(!is.na(file.db.pH))
  {
    par.acidtol <- numeric(0)
    
    for (i in 1:length(Invertebrates))
    {
      par <- construct.traitpars(Invertebrates[i],
                                 cind=2:ncol(db.pH),
                                 method=method.pH,
                                 name.trait="acidtolval",
                                 default=1,
                                 db.traits=db.pH,
                                 db.taxonomy=db.taxonomy)
      
      if(length(par)>0) par.acidtol <- c(par.acidtol,par)
    }  
    
    names(par.acidtol) <- sub(pattern = "ind.1", replacement = "ind", x=names(par.acidtol))
    
    par.traits <- c(par.traits,par.acidtol)
    
  }
  
  
  # Spear:
  # ------ 
  
  if(!is.na(file.db.spear))
  {
    par.orgmicropolltolval <- numeric(0)
    
    for (i in 1:length(Invertebrates))
    {
      parSpear <- construct.traitpars(invert=Invertebrates[i],
                                      cind=13,
                                      method=method.spear,
                                      name.trait="orgmicropolltolval",
                                      default=0,
                                      db.traits=db.spear,
                                      db.taxonomy=db.taxonomy)
      
      names(parSpear) <- paste(Invertebrates[i],"orgmicropolltolval",sep="_")
      
      parSpear2        <- c(NA, NA)
      names(parSpear2) <- c(paste(Invertebrates[i],"orgmicropolltolval", "class1",sep="_"),
                            paste(Invertebrates[i],"orgmicropolltolval", "class2",sep="_"))
      
      if(parSpear == 0) # tolerant species, could live in conditions with or without pesticides
      {
        parSpear2[grep("class1", names(parSpear2))] = 1
        parSpear2[grep("class2", names(parSpear2))] = 1
      }
      if(parSpear == 1) # sensitive species, could live in conditions without pesticides
      {
        parSpear2[grep("class1", names(parSpear2))] = 1
        parSpear2[grep("class2", names(parSpear2))] = 0
      }
      
      if(length(parSpear2)>0) par.orgmicropolltolval <- c(par.orgmicropolltolval,parSpear2)
    }
    
  }
  
  ##############
  ## Food web ##
  ##############
  
  # Food sources  
  # ------------
  
  if(!is.na(file.db.food))
  {
    par.fs <- numeric(0)
    
    for (i in 1:length(Invertebrates))
    {
      par <- construct.traitpars(Invertebrates[i],
                                 cind=2:ncol(db.food),
                                 method=method.food,
                                 name.trait="food",
                                 default=1,
                                 db.traits=db.food,
                                 db.taxonomy=db.taxonomy)
      
      if(length(par)>0) par.fs <- c(par.fs,par)
    }  
    
    par.traits <- c(par.traits,par.fs)
    
  }
  # Feeding type  
  # ------------
  
  if(!is.na(file.db.feeding))
  {
    
    par.ft <- numeric(0)
    
    for (i in 1:length(Invertebrates))
    {
      par <- construct.traitpars(invert=Invertebrates[i],
                                 cind=2:ncol(db.feeding),
                                 method=method.feeding,
                                 name.trait="feedingtype",
                                 default=1,
                                 db.traits=db.feeding,
                                 db.taxonomy=db.taxonomy)
      
      if(length(par)>0) par.ft <- c(par.ft,par)
    }
    
    par.traits <- c(par.traits,par.ft)
    
  }
  
  ## Collect all parameters
  ## ----------------------
  
  if (length(par.traits)>0)  par.traits <- par.traits

  if (length(par.M)>0)       par.traits <- c(par.traits, par.M)
  
  if (length(par.orgmicropolltolval)>0)   par.traits <- c(par.traits, par.orgmicropolltolval)
  
  
  return(par.traits)
}

#invpartraits <- construct.invpars.traits(Invertebrates)