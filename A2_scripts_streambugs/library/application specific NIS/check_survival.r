# 27.2.2023 adapted to run based on streambugs package, see lines 15 and 22


check.survival <- function(res,times,par,Dens.crit=1e-10,D.crit.na=1e-10)                                                                                  
  
  # for testing: res=res.R;times=tout;par=par.fix;Dens.crit=D.crit;D.crit.na=1e-10
    
{
  
  if(length(times)>nrow(res)) stop("length(times) larger than nrow(res), simulations not successful")
  
  sys.def <- streambugs.get.sys.def(y.names=colnames(res)[-1],par=par,inp=inp)
  y.names <- sys.def$y.names
  
  par.envcond.w <- streambugs:::get.inpind.parval.envcond.reach(
    par.names = c("w"),
    y.names   = y.names,
    par       = par,
    inp       = inp,
    required  = c("w"))
  
  par.envcond.fA <- streambugs:::get.inpind.parval.envcond.habitat(
    par.names = c("fA"),
    y.names   = y.names,
    par       = par,
    inp       = inp,
    defaults  = c(fA=1))

  t <- c(ceiling(length(times)*4/5):length(times))
  c.survival <- c(rep(NA,length(y.names$y.names)))
  names(c.survival) <- y.names$y.names
  lind.POM <- y.names$y.groups == "POM"
  c.survival[lind.POM]<- "survived"
  
  # Construct matrix D.mod with results converted from mass per unit length into mass per unit 
  # surface area:
  
  D.mod <- matrix(NA,nrow=length(t),ncol=length(y.names$y.names))
  colnames(D.mod)<-colnames(res)[-1]
  rownames(D.mod)<-times[t]
  
  for ( j in 1:length(t))
  {
    # update w and fA:
    
    inpvals <- streambugs:::interpolate.inputs(inp,res[j,1])
    
    w  <- streambugs:::streambugs.update.envcond.reach(par.envcond.w,inpvals)[,"w"]
    fA <- streambugs:::streambugs.update.envcond.habitat(par.envcond.fA,inpvals,y.names$ind.fA)[,"fA"]
    
    # convert results in mass per unit length into mass per unit 
    # surface area:
    
    D.mod[j,] <- res[t[j],-1]/(w*fA)
  }
  
  # compare modelled density with Density threshold for "survival"
  
  for (i in 1:length(y.names$y.names) )
  {
    if ( !lind.POM[i] ) # alle y die nicht POM sind
    {
      D.i <- mean(D.mod[,y.names$y.names[i]])
      
      if(length(Dens.crit)==1) {D.crit<-Dens.crit}
      if(length(Dens.crit)>1)
      {
        ind <- match(y.names$y.names[i],names(Dens.crit))
        ifelse (is.na(ind), D.crit <- D.crit.na, D.crit <- Dens.crit[ind]) 
      }
      
      ifelse ( D.i < D.crit, c.survival[i] <- "extinct", c.survival[i] <- "survived" )
    }
  }
  
  return(c.survival)
}

#####################

construct.survival.matrix <- function(c.survival,y.names)
  # return a matrix with rows = reach/habitat combinations and cols=taxa
{
  if(!is.list(y.names)) y.names <- decode.statevarnames(y.names)
  
  reaches <- y.names$y.reaches
  habs <- y.names$y.habitats
  reach.hab <- unique(paste(y.names$y.reaches,y.names$y.habitats,sep="."))
  
  survival.matrix <- matrix (NA,nrow=length(reach.hab),ncol=length(y.names$taxa))
  colnames(survival.matrix) <- y.names$taxa
  rownames(survival.matrix) <- reach.hab
  for ( i in 1:length(reach.hab))
  {
    for (j in 1:length(y.names$taxa))
    {
      reach <- unlist(strsplit(reach.hab[i],split="\\."))[1]
      hab   <- unlist(strsplit(reach.hab[i],split="\\."))[2]
      taxon <- y.names$taxa[j]
      
      indr  <- grep(reach,names(c.survival))
      indh  <- grep(hab,names(c.survival))
      indtaxa  <- grep(taxon,names(c.survival))
      indt <- NULL
      for (k in indtaxa)
      {
        splitname <- unlist(strsplit(names(c.survival)[k],split="_"))
        if(sum(taxon==splitname)==1) indt <- c(indt,k)
      }
      
      ind   <- intersect(intersect(indr,indh),indt)
      if(length(ind)>1) warning(taxon," not unique")
      survival.matrix[i,j] <- c.survival[ind]
    }
  }
  return(survival.matrix)
}
  
  



############################################

# survival_prob <- function(res_array,pars,Dens_crit=1e-10,Dcritna=1e-10 )                                                                                    #         cos(2*pi*(times-pars["t.max"])))       # degC
# {
#   states <- unlist(dimnames(res_array)[2])[-1]
#   times <- res_array[,1,1]
#   n_parsamp <- length(res_array[1,1,])
#   t <- c(ceiling(length(times)*4/5):length(times))
#   p_survival <- c(rep(NA,length(states)))
#   names(p_survival) <- states
#   
#   if(length(Dens_crit)==1)
#   {
#    for (i in 1:length(states) )
#    {
#      surv_ind <- rep(NA,n_parsamp)
#      for (j in 1:n_parsamp )
#      {
#        ifelse(sum(res_array[t,states[i],j]/(pars["L"]*pars["w"]),na.rm = T)/length(t) > Dens_crit,
#               surv_ind[j]<-1,surv_ind[j]<-0)   #hier na.rm=T oder so
#      }
#      p_survival[i]<-  round(sum(surv_ind)/n_parsamp,digits=3)
#    }
#   }
#   
#   if(length(Dens_crit)>1)
#   {   
#    for (i in 1:length(states) )
#    {
#      ind <- match(names(state)[i],names(Dens_crit))
#      if(is.na(ind)) Dcrit=Dcritna else Dcrit <- Dens_crit[ind]
#      
#      surv_ind <- rep(NA,n_parsamp)
#      for (j in 1:n_parsamp )
#      {
#        ifelse(sum(res_array[t,states[i],j]/(pars["L"]*pars["w"]),na.rm = T)/length(t) > Dcrit,
#               surv_ind[j]<-1,surv_ind[j]<-0)   #hier na.rm=T oder so
#      }
#      p_survival[i]<-  round(sum(surv_ind)/n_parsamp,digits=3)
#    }
#   }
#  
#   return(p_survival)
# }

############################################

# check_survival_prob <- function(p_survival,Threshold=0.5)
# {
#   c_survival <- rep(NA,length(p_survival))
#   names(c_survival) <- names(p_survival)
#   for (i in 1:length(p_survival) )
#   {
#     ifelse ( p_survival[i] < Threshold,
#              c_survival[i]<- "extinct",  c_survival[i]<- "survived"  )
#   }
#   
#   return(c_survival)
# }

############################################

# count_foodwebs <- function(res_array,pars,Dens_crit=1e-10,Dcritna=1e-10) 
# {
#   states <- c("SusPOM",unlist(dimnames(res_array)[2])[-1]  )
#   times <- res_array[,1,1]
#   n_parsamp <- length(res_array[1,1,])
#   n_not_na   <- n_parsamp
#   t <- c(ceiling(length(times)*4/5):length(times))
#   fwstrings <- matrix(NA,nrow=1,ncol=2)
#   colnames(fwstrings) <- c("fwstrings","counts")
#  
#   surv_ind_na <- c("S",rep(NA,length=length(states)-1))
#   nastring <- paste(surv_ind_na,collapse="")
#   
#   for (n in 1:n_parsamp)
#   {
#     surv_ind <- rep(NA,length(states))
#     surv_ind[1] <- "S"
#     
#     if(length(Dens_crit)==1)
#     {
#       for (i in 2:length(states) )
#       {
#         surv_ind[i] <- ifelse(sum(res_array[t,states[i],n]/(pars["L"]*pars["w"]))/length(t) > Dens_crit,
#                               "S","E") 
#       }               
#     }
#     
#     if(length(Dens_crit)>1)
#     {
#       for (i in 2:length(states) )
#       {
#         ind <- match(names(state)[i],names(Dens_crit))
#         if(is.na(ind)) Dcrit=Dcritna else Dcrit <- Dens_crit[ind]
#         
#         surv_ind[i] <- ifelse(sum(res_array[t,states[i],n]/(pars["L"]*pars["w"]))/length(t) > Dcrit,
#                               "S","E") 
#       }               
#     }
#     
#     fwstring  <- paste(surv_ind,collapse="")             
#     
#     namatch <- (nastring==fwstring)
#     
#     if(namatch) {n_not_na <- n_not_na-1}
#     
#     if(!namatch)
#     {
#      rindmatch <- grep(fwstring,fwstrings[,1])
#                 
#      if (sum(rindmatch) > 0)
#      {
#       fwstrings[rindmatch,2] <- as.numeric(fwstrings[rindmatch,2])+1
#      }
#      if (sum(rindmatch) == 0)
#      {
#       fwstrings <- rbind(fwstrings, c(fwstring,1))
#      }
#     }
#   }
#   
#    
#   
#   fwstrings[,2] <- as.numeric(fwstrings[,2])/n_not_na
#   fwstrings <- fwstrings[-1,]
#   if( is.matrix(fwstrings) ) 
#   { 
#     o <- order(as.numeric(fwstrings[,2]),decreasing=TRUE)
#     fwstrings <- fwstrings[o,]
#   }
#   return(fwstrings)
# }

############################################

paste.survival <- function(c.survival)              #c_survival=survivals
{
  survivors <- rep(NA,sum(c.survival=="survived"))
  
  ind <- (c.survival == "survived")
  survivors <- names(c.survival[ind])
  cat("survived: ",survivors,fill=T)
  
  ind2 <- (c.survival == "extinct")
  extincts <- names(c.survival[ind2])
  cat("\ngot extinct: ",extincts,fill=T)
}

############################################
# 
# paste_factors_extincts <- function(c_survival,pars,res)
# {
#    #c_survival=survivals; pars=pars; res=res
#    
#    
#    ind        <- (c_survival == "extinct")
#    extincts   <- names(c_survival[ind])
#    if(length(extincts)==0) cat("no extincts \n") else 
#    {
#      f_extincts <- matrix(data=NA,nrow=length(extincts),ncol=7)
#      colnames(f_extincts) <- c("ext taxa","fcurrent","fmicrohab","ftemp","fspear","fsapro","foodlim")
#      
#      p <- strsplit(names(pars),split="_")
#      par1_state <- rep(NA,length(p))
#      par2_name  <- rep(NA,length(p))
#      par3  <- rep(NA,length(p))
#      par4  <- rep(NA,length(p))
#      for ( i in 1:length(p) )
#      {
#         par1_state[i] <- p[[i]][1]
#         par2_name[i]  <- p[[i]][2]
#         par3[i]       <- p[[i]][3]
#         par4[i]       <- p[[i]][4]
#      }
#      A   <- pars["L"]*pars["w"] 
#      
#      for (i in 1:length(extincts) )
#      {
#         invert <- extincts[i]
#         f_extincts[i,1] <- invert
#         f_extincts[i,2] <- pars[paste(invert,"_fcurrent",sep="")]
#         f_extincts[i,3] <- pars[paste(invert,"_fmicrohab",sep="")]
#         f_extincts[i,4] <- pars[paste(invert,"_ftemp",sep="")]
#         f_extincts[i,5] <- pars[paste(invert,"_fspear",sep="")]
#         f_extincts[i,6] <- pars[paste(invert,"_fsapro",sep="")]
#      
#         foodind <- which((par1_state == invert) & (par2_name=="K"))
#         foods   <- unique(par3[foodind])
#         sumfood <- 0
#         if(length(foods)>0)
#         {
#           for (j in 1:length(foods))
#           {
#             food <- foods[j]
#             
#             if(food=="I")      cfood <- pars["I0"]*(1-pars["shade"]) else  
#             if(food=="SusPOM") cfood <- pars["D_SusPOM"]             else             
#                                cfood <- res[nrow(res),food]/A 
#             
#             sumfood <- sumfood + cfood
#           }
#           K     <- pars[paste(invert,"_K_",foods[1],sep="")]        #nur erstes K genommen
#           f_extincts[i,7] <- foodlim <- round(sumfood/(K+sumfood),digits=2)  
#         }
#      }
#      return(f_extincts)
#    }
# }
# 
# ############################################
# 
# paste_factors_survivers <- function(c_survival,pars,res)
# {
#    ind        <- (c_survival == "survived")
#    survivers   <- names(c_survival[ind])
#    f_survivers <- matrix(data=NA,nrow=length(survivers),ncol=7)
#    colnames(f_survivers) <- c("surv taxa","fcurrent","fmicrohab","ftemp","fspear","fsapro","foodlim")
#    
#    p <- strsplit(names(pars),split="_")
#    par1_state <- rep(NA,length(p))
#    par2_name  <- rep(NA,length(p))
#    par3  <- rep(NA,length(p))
#    par4  <- rep(NA,length(p))
#    for ( i in 1:length(p) )
#    {
#       par1_state[i] <- p[[i]][1]
#       par2_name[i]  <- p[[i]][2]
#       par3[i]       <- p[[i]][3]
#       par4[i]       <- p[[i]][4]
#    }
#    A   <- pars["L"]*pars["w"] 
#    
#    for (i in 1:length(survivers) )
#    {
#       invert <- survivers[i]
#       f_survivers[i,1] <- invert
#       f_survivers[i,2] <- pars[paste(invert,"_fcurrent",sep="")]
#       f_survivers[i,3] <- pars[paste(invert,"_fmicrohab",sep="")]
#       f_survivers[i,4] <- pars[paste(invert,"_ftemp",sep="")]
#       f_survivers[i,5] <- pars[paste(invert,"_fspear",sep="")]
#       f_survivers[i,6] <- pars[paste(invert,"_fsapro",sep="")]
#       
#       foodind <- which((par1_state == invert) & (par2_name=="K"))
#       foods   <- unique(par3[foodind])
#       sumfood <- 0
#       if(length(foods)>0)
#       {
#         for (j in 1:length(foods))
#         {
#           food <- foods[j]
#           
#           if(food=="I")      cfood <- pars["I0"]*(1-pars["shade"]) else 
#           if(food=="SusPOM") cfood <- pars["D_SusPOM"]             else    
#                              cfood <- res[nrow(res),food]/A 
#           
#           sumfood <- sumfood + cfood
#         }
#         K     <- pars[paste(invert,"_K_",foods[1],sep="")]        #nur erstes K genommen
#         f_survivers[i,7] <- foodlim <- round(sumfood/(K+sumfood),digits=2)  
#       }
#    }
#    return(f_survivers)
# }


  