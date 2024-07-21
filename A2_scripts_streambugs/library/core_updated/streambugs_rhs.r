################################################################
#
# streambugs 1.1dev
# =================
#
# -----------------------------------------
# right-hand side of differential equations
# -----------------------------------------
#
# creation:      06.07.2012 
# modifications: 30.04.2013
#
################################################################



# ==============================================================
# R implementation of the function to define the right-hand side 
# of the odes
# ==============================================================


rhs.streambugs <- function(t,y,parms,add.out=FALSE)
{
  if ( !exists("streambugs.debug") ) streambugs.debug <- 0

if (add.out==TRUE) res.add <- NA
  
  if(sum(is.na(y))>0) stop ("One or more state variables are NA: ", names(y[is.na(y)])," at time t=",t,"\n",sep="  ")
  
  # update time.dependent parameters:
  # ---------------------------------
  
  sys.def <- streambugs.updatepar(parms,t)   # parms is sys.def at initial time
  
  # facilitate usage of some of the parameter arrays:
  # -------------------------------------------------
  
  y.names          <- sys.def$y.names
  
  par.global       <- sys.def$par.global$parvals
  
  par.global.envtraits <- sys.def$par.global.envtraits
  
  par.envcond      <- cbind(sys.def$par.envcond.reach$parvals,
                            sys.def$par.envcond.habitat$parvals)
  
  par.envcond.habitat.group <- sys.def$par.envcond.habitat.group
  
  par.input        <- sys.def$par.input$parvals
  
  par.taxaprop     <- sys.def$par.taxaprop.direct$parvals
  
  par.taxaprop.traits <- sys.def$par.taxaprop.traits
    
  y.input <- par.input[,"Input"] * par.envcond[,"w"] * par.envcond[,"fA"]

  par.trait.cond <- get.par.trait.cond(par.taxaprop.traits       = par.taxaprop.traits,
                                       par.envcond               = par.envcond,
                                       par.envcond.habitat.group = par.envcond.habitat.group,
                                       par.global.envtraits      = par.global.envtraits,
                                       par                       = par.global)
  
  
  
  # make local copy of states to handle negative states
  # ---------------------------------------------------
  
  y.original <- y
  
  if ( streambugs.debug > 0 )
  { 
    if (sum(is.na(y))>0)
    {
      ind.na <- which(is.na(y))
      
      y.na <- y.names$y.names[ind.na]
      
      cat("state variable ",y.na," is na at t=",t) 
    }
  }
  
  y <- ifelse(y<0,0,y)
  
  
  # calculate basal metabolic rates (gDM/m2/a; will be NA for dead material):    
  # -------------------------------------------------------------------------
  
  r.basal <- par.taxaprop[,"i0"]*
             par.taxaprop[,"fbasaltax"]*
             (par.taxaprop[,"M"]^(par.taxaprop[,"b"]-1))/(par.global["M0"]^par.taxaprop[,"b"])*
             exp(-par.taxaprop[,"Ea"]/(par.global["kBoltzmann"]*par.envcond[,"T"]))/
             par.taxaprop[,"EC"]*
             y/(par.envcond[,"fA"]*par.envcond[,"w"])
  
if (add.out==TRUE) res.add <- c(res.add,r_basal=r.basal)
  
  # calcuate right-hand side of odes: 
  # ---------------------------------
  
  # input
  
  dydt <- y.input
  names(dydt) <- names(y)
  
  # processes
  
  for ( i in 1:length(y) )
  {
    w.fA <- par.envcond[i,"w"] * par.envcond[i,"fA"]
    
if ( streambugs.debug > 0 )
{
if ( t==0 ) 
{
  cat("i = ",i,", w = ",par.envcond[i,"w"]," fA = ",par.envcond[i,"fA"],"\n",sep="")
  cat("i = ",i,", r.basal = ",r.basal[i],"\n",sep="")
}
}
    
    # drift induced by floods:
    # ------------------------
    if ( length(sys.def$par.proc.taxon[[i]][["Drift"]]) > 0 )
    {
      
      tau <- par.envcond[i,"tau"]
      if ( is.na(tau) )
      {
        warning("Stoichiometry defines drift for component ",i,
                " but parameter \"tau\" is not defined",
                " for component ",i," of state vector")
      }
      
      taucrit <- par.envcond[i,"taucrit"]
      if ( is.na(tau) )
      {
        warning("Stoichiometry defines drift for component ",i,
                " but parameter \"taucrit\" is not defined",
                " for component ",i," of state vector")
      }
      
      if ( !is.na(tau) & !is.na(taucrit) )
      {
        
        # get parameter values:
        
        par.proc1 <- sys.def$par.proc.taxon[[i]][["Drift"]]$parvals
        
        # calculate rate:
        
        r.drift <- 0
        
        if(tau >= taucrit)
        {
          r.drift <- par.proc1["cdet"] * y[i]/w.fA * (tau-taucrit)^2  # gDM/m2/a
        }

if (add.out==TRUE) 
{
 names(r.drift) <- names(y[i])
 res.add <- c(res.add,r.drift=r.drift)
}
        
        
if ( streambugs.debug > 0 )
{
 if ( t==0 ) 
 {
  cat("i = ",i,", r.drift = ",r.drift,"\n",sep="")
 }
}
        
        # calculate contributions to rate of change:
        
        stoich <- sys.def$par.proc.taxon[[i]][["Drift"]]$stoich
        for ( j in 1:ncol(stoich) )
        {
          if ( !is.na(stoich[1,j]) ) 
          {
            dydt[stoich[1,j]] <- dydt[stoich[1,j]] + stoich[2,j]*r.drift*w.fA
          }
        } 
      } 
    }
    
    
    # mineralization:
    # ---------------
    
    if ( length(sys.def$par.proc.taxon[[i]][["Miner"]]) > 0 )
    {
      # get parameter values:
      
      par.proc1 <- sys.def$par.proc.taxon[[i]][["Miner"]]$parvals
      
      # calculate rate:

      r.miner <- par.proc1["kminer"]*y[i]/w.fA   # gDM/m2/a

if (add.out==TRUE) 
{
  names(r.miner) <- names(y[i])
  res.add <- c(res.add,r_miner=r.miner)
}
  

if ( streambugs.debug > 0 )
{
if ( t==0 ) 
{
cat("i = ",i,", r.miner = ",r.miner,"\n",sep="")
}
}
      
      # calculate contributions to rate of change of all affected taxa:
      
      stoich <- sys.def$par.proc.taxon[[i]][["Miner"]]$stoich
      for ( j in 1:ncol(stoich) )
      {
        if ( !is.na(stoich[1,j]) ) 
        {
          dydt[stoich[1,j]] <- dydt[stoich[1,j]] + stoich[2,j]*r.miner*w.fA
        }
      }
    }

    # respiration:
    # ------------
    
    if ( length(sys.def$par.proc.taxon[[i]][["Resp"]]) > 0 )
    {
      # get parameter values:
      
      par.proc1 <- sys.def$par.proc.taxon[[i]][["Resp"]]$parvals
      
      # calculate rate:

      r.resp <- par.proc1["fresp"]*r.basal[i]  # gDM/m2/a

if (add.out==TRUE) 
{
  names(r.resp) <- names(y[i])
  res.add <- c(res.add,r_resp=r.resp)
} 
  
  
      
if ( streambugs.debug > 0 )
{
if ( t==0 ) 
{
cat("i = ",i,", r.resp = ",r.resp,"\n",sep="")
}
}
      
      # calculate contributions to rate of change of all affected taxa:
      
      stoich <- sys.def$par.proc.taxon[[i]][["Resp"]]$stoich
      for ( j in 1:ncol(stoich) )
      {
        if ( !is.na(stoich[1,j]) ) 
        {
          dydt[stoich[1,j]] <- dydt[stoich[1,j]] + stoich[2,j]*r.resp*w.fA
        }
      }
    }
    
    # death:
    # ------
    
    if ( length(sys.def$par.proc.taxon[[i]][["Death"]]) > 0 )
    {
      # get parameter values:
      
      par.proc1 <- sys.def$par.proc.taxon[[i]][["Death"]]$parvals
           
      fsapro     <- par.trait.cond[rownames(par.trait.cond)==names(y[i]),"fsapro" ] 
      forgmicropoll <- par.trait.cond[rownames(par.trait.cond)==names(y[i]),"forgmicropoll" ]

if (add.out==TRUE)       
{
  names(fsapro)     <- names(y[i])
  names(forgmicropoll) <- names(y[i])
  res.add <- c(res.add,fsapro=fsapro,forgmicropoll=forgmicropoll)
} 
      
      
      
      
if ( streambugs.debug > 0 )
{
if ( t==0 ) 
{
cat("i = ",i,", fsapro = ",fsapro,", forgmicropoll = ",forgmicropoll,"\n",sep="")
}
}
         
      # calculate rate:

      r.death <- par.proc1["fdeath"]*r.basal[i]*fsapro*forgmicropoll  # gDM/m2/a

if (add.out==TRUE) 
{
  names(r.death)     <- names(y[i])
  res.add <- c(res.add,r_death=r.death)
}  

if ( streambugs.debug > 0 )
{
if ( t==0 ) 
{
cat("i = ",i,", r.death = ",r.death,"\n",sep="")
}
}

      # calculate contributions to rate of change of all affected taxa:
      
      stoich <- sys.def$par.proc.taxon[[i]][["Death"]]$stoich
      for ( j in 1:ncol(stoich) )
      {
        if ( !is.na(stoich[1,j]) ) 
        {
          dydt[stoich[1,j]] <- dydt[stoich[1,j]] + stoich[2,j]*r.death*w.fA
        }
      }
    }

    # primary production:
    # -------------------
    
    if ( length(sys.def$par.proc.taxon[[i]][["Prod"]]) > 0 )
    {
      # get parameter values:
      
      par.proc1 <- sys.def$par.proc.taxon[[i]][["Prod"]]$parvals
      
      fgrotax <- par.proc1["fgrotax"]
      
      I0 <- par.envcond[i,"I0"]
      if ( is.na(I0) )
      {
        warning("Stoichiometry defines Production for component ",i,
                " but parameter \"I0\" is not defined",
                " for component ",i," of state vector")
      }
      CP <- par.envcond[i,"CP"]
      if ( is.na(CP) )
      {
        warning("Stoichiometry defines Production for component ",i,
                " but parameter \"CP\" is not defined",
                " for component ",i," of state vector")
      }
      CN <- par.envcond[i,"CN"]
      if ( is.na(CN) )
      {
        warning("Stoichiometry defines Production for component ",i,
                " but parameter \"CN\" is not defined",
                " for component ",i," of state vector")
      }
      fshade <- par.envcond[i,"fshade"]
      if ( is.na(fshade) )
      {
        warning("Stoichiometry defines Production for component ",i,
                " but parameter \"fshade\" is not defined",
                " for component ",i," of state vector")
      }
      
      # calculate rate:
      
      flimI         <- I0/(par.proc1["KI"]+I0)
      flimP         <- CP/(par.proc1["KP"]+CP)
      flimN         <- CN/(par.proc1["KN"]+CN)
      flimnutrients <- min(flimP,flimN)
      fselfshade    <- par.proc1["hdens"]/(par.proc1["hdens"]+y[i]/w.fA)
      
      r.prod <- par.proc1["fprod"]*
                fgrotax*
                flimI*
                flimnutrients*
                fselfshade*
                (1-fshade)*
                r.basal[i]   # gDM/m2/a
      

if (add.out==TRUE)
{
  names(flimI) <- names(flimP) <- names(flimN) <- names(flimnutrients) <- 
    names(fselfshade) <- names(r.prod) <- names(y[i])
  
  res.add <- c(res.add,flimI=flimI,flimP=flimP,flimN=flimN,
               flimnutrients=flimnutrients,fselfshade=fselfshade,r_prod=r.prod)
}   

if ( streambugs.debug > 0 )
{
if ( t==0 ) 
{
cat("i = ",i,", r.prod = ",r.prod,"\n",sep="")
}
}
      
      # calculate contributions to rate of change of all affected taxa:
      
      stoich <- sys.def$par.proc.taxon[[i]][["Prod"]]$stoich
      for ( j in 1:ncol(stoich) )
      {
        if ( !is.na(stoich[1,j]) ) 
        {
          dydt[stoich[1,j]] <- dydt[stoich[1,j]] + stoich[2,j]*r.prod*w.fA
        }
      }
    }

    # consumption:
    # ------------
    
    if ( length(sys.def$par.proc.web[[i]][["Cons"]]) > 0 )
    {
      # get parameter values:
      
      par.proc1 <- sys.def$par.proc.web[[i]][["Cons"]]$parvals
      taxa2     <- names(sys.def$par.proc.web[[i]][["Cons"]]$taxa2)
      
      fcurrent  <- par.trait.cond[rownames(par.trait.cond)==names(y[i]),"fcurrent"]
      ftempmax  <- par.trait.cond[rownames(par.trait.cond)==names(y[i]),"ftempmax"]
      fmicrohab <- par.trait.cond[rownames(par.trait.cond)==names(y[i]),"fmicrohab"]

        
if ( streambugs.debug > 0 )
{
if ( t==0 ) 
{
cat("i = ",i,", fcurrent = ",fcurrent,", ftempmax = ",ftempmax,", fmicrohab = ",fmicrohab,"\n",sep="")
}
}
      Kdens   <- par.proc1["hdens"] * fcurrent * ftempmax * fmicrohab
      fgrotax <- par.proc1["fgrotax"]
      Kfood   <- par.proc1["Kfood"]
      q       <- par.proc1["q"] 
      
      # get preferences and food densities:
      
      y.food <- rep(0,length(taxa2))   #XXX mod 26.4.
      pref   <- rep(NA,length(taxa2))
      for ( k in 1:length(taxa2) )
      {
        par.proc2   <- sys.def$par.proc.web[[i]][["Cons"]]$taxa2[[k]]$parvals
        stoich      <- sys.def$par.proc.web[[i]][["Cons"]]$taxa2[[k]]$stoich
        if ( taxa2[k] == "SusPOM" )
        {
          y.food[k] <- par.envcond[i,"DSusPOM"]*w.fA
          if ( is.na(y.food[k]) )
          {
            warning("Stoichiometry defines Consumption of \"DSusPOM\" for component ",i,
                    " but parameter \"DSusPOM\" is not defined",
                    " for component ",i," of state vector")
          }
        }
        else
        {           
          ind.food  <- stoich[1,match(taxa2[k],colnames(stoich))]
          y.food[k] <- y[ind.food]
        }
        pref[k] <- par.proc2["Pref"]
      }
      if( sum(is.na(y.food)) > 0 ) 
      {
        ind.na.food <- which(is.na(y.food))  #xxx added 26.4.2013
        stop("One or more consumption food sources: ",taxa2[ind.na.food]," of state variable ",i," is NA at time",t, sep=" ")
      }
      sum.food <- sum(y.food)
      sum.food.pref <- sum(y.food*pref)
      
      # consumption rate without preference factor 
      
      fselfinh <- Kdens/(Kdens+y[i]/w.fA)
      ffoodlim <- ifelse(sum.food>=0,(sum.food/w.fA)^q/(Kfood^q+(sum.food/w.fA)^q),(sum.food/w.fA)^q/(Kfood^q))
      
      r.cons.tot <- par.proc1["fcons"]*
                    fgrotax*                                
                    fselfinh*
                    ffoodlim*
                    r.basal[i]  # gDM/m2/a

if (add.out==TRUE) 
{
  names(fcurrent) <- names(ftempmax) <- names(fmicrohab) <- names(y[i])
  names(sum.food) <- names(sum.food.pref) <-  names(y[i])
  names(fselfinh) <- names(ffoodlim) <-  names(r.cons.tot) <-names(y[i])
  
  res.add <- c(res.add,fcurrent=fcurrent,ftempmax=ftempmax,fmicrohab=fmicrohab,
               sum_food=sum.food,sum_food_pref=sum.food.pref,
               fselfinh=fselfinh,ffoodlim=ffoodlim,r_cons_tot=r.cons.tot)
}
      
if ( streambugs.debug > 0 )
{
if ( t==0 ) 
{
  cat("i = ",i,", r.cons.tot = ",r.cons.tot,
      ", sum.food = ",sum.food,", sum.food.pref = ",sum.food.pref,"\n",sep="")
}
}
      
      # loop over all food types:
      
      for ( k in 1:length(taxa2) )
      {
        # calculate rate:
        
        if(sum.food.pref>0) 
        {
          fpref <- y.food[k]*pref[k]/sum.food.pref
        } else (fpref=0)
        
        
        r.cons <- r.cons.tot*fpref  # gDM/m2/a

if (add.out==TRUE) 
{
  names(fpref) <- names(r.cons) <- paste(names(y[i]),taxa2[k],sep="__")
  res.add <- c(res.add,fpref=fpref,r_cons=r.cons)
}
        
        # calculate contributions to rate of change of all affected taxa:
        
        stoich <- sys.def$par.proc.web[[i]][["Cons"]]$taxa2[[k]]$stoich
        for ( j in 1:ncol(stoich) )
        {
          if ( !is.na(stoich[1,j]) ) 
          {
            dydt[stoich[1,j]] <- dydt[stoich[1,j]] + stoich[2,j]*r.cons*w.fA  
          }
        }  
      } 
    }    
    
    # fish predation:
    # ---------------
    
    if ( length(sys.def$par.proc.web[[i]][["FishPred"]]) > 0 )
    {
      # get parameter values:
      
      par.proc1 <- sys.def$par.proc.web[[i]][["FishPred"]]$parvals
      taxa2     <- names(sys.def$par.proc.web[[i]][["FishPred"]]$taxa2)
      
      DFish <- par.envcond[i,"DFish"]
      cfish <- par.proc1["cfish"]
      Kfood <- par.proc1["Kfood"]
      q     <- par.proc1["q"] 
      
      # get preferences and food densities:
      
      y.food <- rep(NA,length(taxa2))
      pref   <- rep(NA,length(taxa2))
      for ( k in 1:length(taxa2) )
      {
        par.proc2 <- sys.def$par.proc.web[[i]][["FishPred"]]$taxa2[[k]]$parvals
        stoich    <- sys.def$par.proc.web[[i]][["FishPred"]]$taxa2[[k]]$stoich
        if ( taxa2[k] == "SusPOM" )
        {
          y.food[k] <- par.envcond[i,"DSusPOM"]*w.fA
          if ( is.na(y.food[k]) )
          {
            warning("Stoichiometry defines fish consumption of \"DSusPOM\" for component ",i,
                    " but parameter \"DSusPOM\" is not defined",
                    " for component ",i," of state vector")
          }
        }
        else
        {           
          ind.food  <- stoich[1,match(taxa2[k],colnames(stoich))]
          y.food[k] <- y[ind.food]
        }
        pref[k] <- par.proc2["Pref"]
      }
      if( sum(is.na(y.food)) > 0 ) 
      {
        stop("One or more fish predation food sources for state variable ",i," are NA")
      }
      sum.food <- sum(y.food)
      sum.food.pref <- sum(y.food*pref)
      
      # consumption rate without preference factor 
      
      ffoodlim <- (sum.food/w.fA)^q/(Kfood^q+(sum.food/w.fA)^q)
      
      r.fishpred.tot <- DFish/10000 * cfish*365.25 * ffoodlim  # gDM/m2/a
      
if (add.out==TRUE) 
{
  names(sum.food) <- names(sum.food.pref) <-  names(y[i])
  names(ffoodlim) <- names(r.fishpred.tot) <- names(y[i])
  
  res.add <- c(res.add,sum_food_fish=sum.food,sum_food_pref_fish=sum.food.pref,
               ffoodlim_fish=ffoodlim,r_fishpred_tot=r.fishpred.tot)
}
      
if ( streambugs.debug > 0 )
{
  if ( t==0 ) 
  {
    cat("i = ",i,", r.fishpred.tot = ",r.fishpred.tot,
        ", sum.food = ",sum.food,", sum.food.pref = ",sum.food.pref,"\n",sep="")
  }
}
      
      # loop over all food types:
      
      for ( k in 1:length(taxa2) )
      {
        # calculate rate:
        
        if(sum.food.pref>0) 
        {
          fpref <- y.food[k]*pref[k]/sum.food.pref
        } else (fpref=0)
        
        r.fishpred <- r.fishpred.tot*fpref  # gDM/m2/a
        
if (add.out==TRUE) 
{
  names(fpref) <- names(r.fishpred) <- paste(names(y[i]),names(y[k]),sep="__")
  res.add <- c(res.add,fpref_fish=fpref,r_fishpred=r.fishpred)
}
        
        # calculate contributions to rate of change of all affected taxa:
        
        stoich <- sys.def$par.proc.web[[i]][["FishPred"]]$taxa2[[k]]$stoich
        for ( j in 1:ncol(stoich) )
        {
          if ( !is.na(stoich[1,j]) ) 
          {
            dydt[stoich[1,j]] <- dydt[stoich[1,j]] + stoich[2,j]*r.fishpred*w.fA  
          }
        }  
      } 
    }    
  }
if ( streambugs.debug > 0 )
{
if ( t==0 ) 
{
for ( i in 1:length(dydt) )
{
cat("i = ",i," dydt = ",dydt[i],"\n",sep="")
}
}
}
  
  # handle negative states: 
  
  dydt <- ifelse( y.original<0 , as.numeric(-5*y.original), dydt )
  
  dydt <- ifelse( is.na(y.original) | is.nan(y.original), 1e-20, dydt )
  
  error <- ifelse( is.na(y.original) | is.nan(y.original), 1, 0 )
  
  if(sum(error)>0)
  {
    ind <- grep(1,error)
    cat(date()," : ",names(error)[ind]," : ",y.original[ind],"\n",file="states_errors.log",append=T)
  }       
  
  if(sum(is.na(dydt))>0) stop ("One or more dydt are NA: ", paste(names(dydt)[is.na(dydt)],y[is.na(dydt)],sep=" ")," at time t=",t,"\n")
  
  if(sum(is.na(y))>0) stop ("One or more state variables are NA: ", names(y[is.na(y)])," at time t=",t,"\n")
  
  y <- y.original
  
  if (add.out==FALSE)   {return(list(dydt))} else  {return(list(res.add))}  
}

