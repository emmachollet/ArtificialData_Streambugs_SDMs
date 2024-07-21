get.par.trait.cond <- function(par.taxaprop.traits,
                               par.envcond,
                               par.envcond.habitat.group,
                               par.global.envtraits,
                               par)
{
  
  par.trait.cond <- matrix(NA,nrow=nrow(par.envcond),ncol=5)
  colnames(par.trait.cond) <- c("fsapro","forgmicropoll","fcurrent","ftempmax","fmicrohab")
  rownames(par.trait.cond) <- rownames(par.envcond)
  
  # saprobic conditions: 
  # --------------------
  
  s.x <- par.global.envtraits$saprowqclassval$parvals
  
  if ( length(s.x) > 0 &  sum( is.na(s.x)) == 0 )
  {
    for (i in 1:nrow(par.envcond))  # loop over statevariables
    {
      par.trait.cond[i,"fsapro"] <- 1
      s.y <- par.taxaprop.traits$saprotolval$parval[i,]      
      if( sum(is.na(s.y)) == 0 )
      {
        if (!is.na(par.envcond[i,"saprowqclass"]))  
        {
          
          # # linear interpolation:
          # 
          # f.sap.i.untrans <- approx(s.x,s.y,xout=par.envcond[i,"saprowqclass"],rule=2)$y
          
          # XXX ECR ####
          # polynomial interpolation:
          
          f.sap.i.untrans <- pracma::pchip(xi = s.x, yi = s.y, x = par.envcond[i,"saprowqclass"])
          
          # exponential transformation:
          
          par.trait.cond[i,"fsapro"] <- exp.transform(f.sap.i.untrans,
                                                      intercept=par["fsapro_intercept"],
                                                      curv=par["fsapro_curv"])
        } else warning("saprowqclass of", i ,"is NA"   )        
      }  else {
        if(grepl("_Invertebrates",rownames(par.envcond)[i]))
        { 
          warning("not all sapro trait classes in par.taxaprop found for ",rownames(par.envcond)[i])
        }
      }
    }
  } else { 
    par.trait.cond[,"fsapro"] <- 1
    warning( "limitation by saprobic conditions not considered, no midvalues of trait classes provided")
  }
  
  
  # organic toxicants:   
  # ------------------  
  
  o.x <- par.global.envtraits$orgmicropollTUval$parvals
  
  if (length(o.x) > 0 & sum( is.na(o.x)) == 0 ) 
  {
    for (i in 1:nrow(par.envcond))  # loop over statevariables
    {
      par.trait.cond[i,"forgmicropoll"] <- 1
      o.y <- par.taxaprop.traits$orgmicropolltolval$parval[i,]      
      if( sum(is.na(o.y)) == 0 )
      {
        if (!is.na(par.envcond[i,"orgmicropollTU"]))  
        {
          # linear interpolation:

          f.o.i.untrans <- approx(o.x,o.y,xout=par.envcond[i,"orgmicropollTU"],rule=2)$y

          # # XXX ECR ####
          # keep linear interpolation for organic micropollutants because only 2 classes (and minimum 3 are needed for polynomial interpolation)
          # # polynomial interpolation:
          # 
          # f.o.i.untrans <- pracma::pchip(xi = o.x, yi = o.y, x = par.envcond[i,"orgmicropollTU"])
          # 
          # exponential transformation:
          
          par.trait.cond[i,"forgmicropoll"] <- exp.transform(f.o.i.untrans,
                                                             intercept=par["forgmicropoll_intercept"],
                                                             curv=par["forgmicropoll_curv"])
          
        } else warning("orgmicropollTU of", i ,"is NA"   )    
      }  else {
        if(grepl("_Invertebrates",rownames(par.envcond)[i]))
        {
          warning("not all orgmircopoll trait classes in par.taxaprop found for",rownames(par.envcond)[i] )
        }
      }
    }
  } else { 
    par.trait.cond[,"forgmicropoll"] <- 1
    warning( "limitation by organic micropollutants not considered, no midvalues of trait classes provided")
  }
  
  # current preference: 
  # -------------------
  
  c.x <- par.global.envtraits$currentmsval$parvals
  
  if ( length(c.x) > 0 &  sum( is.na(c.x)) == 0 )
  {
    for (i in 1:nrow(par.envcond))  # loop over statevariables
    {
      par.trait.cond[i,"fcurrent"] <- 1
      c.y <- par.taxaprop.traits$currenttolval$parval[i,]      
      if( sum(is.na(c.y)) == 0 )
      {
        if (!is.na(par.envcond[i,"currentms"]))  
        {
          
          # # linear interpolation:
          # 
          # f.cur.i.untrans <- approx(c.x,c.y,xout=par.envcond[i,"currentms"],rule=2)$y
          
          # XXX ECR ####
          # polynomial interpolation:
          
          f.cur.i.untrans <- pracma::pchip(xi = c.x, yi = c.y, x = par.envcond[i,"currentms"])
          
          # exponential transformation:
          
          par.trait.cond[i,"fcurrent"] <- exp.transform(f.cur.i.untrans,
                                                        intercept=par["fcurrent_intercept"],
                                                        curv=par["fcurrent_curv"])
          
        } else warning("currentms of", i ,"is NA"   )
      }  else {
        if(grepl("_Invertebrates",rownames(par.envcond)[i]))
        {
          warning("not all current trait classes in par.taxaprop found for",rownames(par.envcond)[i] )
        }
      }
    }
  } else {
    par.trait.cond[,"fcurrent"] <- 1
    warning( "limitation by current conditions not considered, no midvalues of trait classes provided")
  }
  
  
  # temperature tolerance:   
  # ----------------------
  
  t.x <- par.global.envtraits$tempmaxKval$parvals
  
  if ( length(t.x) > 0 &  sum( is.na(t.x)) == 0 )
  {    
    for (i in 1:nrow(par.envcond))  # loop over statevariables
    {
      par.trait.cond[i,"ftempmax"] <- 1
      t.y <- par.taxaprop.traits$tempmaxtolval$parval[i,]      
      if( sum(is.na(t.y)) == 0 )  # all affinity scores are available
      {
        if (!is.na(par.envcond[i,"tempmaxK"]))  
        {
          
          # # linear interpolation:
          # 
          # f.temp.i.untrans <- approx(t.x,t.y,xout=par.envcond[i,"tempmaxK"],rule=2)$y
          
          # XXX ECR ####
          # polynomial interpolation:
          
          f.temp.i.untrans <- pracma::pchip(xi = t.x, yi = t.y, x = par.envcond[i,"tempmaxK"])
          
          # exponential transformation:
          
          par.trait.cond[i,"ftempmax"] <- exp.transform(f.temp.i.untrans,
                                                        intercept=par["ftempmax_intercept"],
                                                        curv=par["ftempmax_curv"])
          
        } else warning("tempmaxK of ", i ," is NA"   )        
      }  else { 
        if(grepl("_Invertebrates",rownames(par.envcond)[i]))
        {
          warning("not all temperature trait classes in par.taxaprop found")
        }
      }
    }
  } else {
    par.trait.cond[,"ftempmax"] <- 1
    warning( "limitation by temperature conditions not considered, no midvalues of trait classes provided")
  }
  
  # microhabitat: 
  # -------------
  
  par.trait.cond[,"fmicrohab"] <- 1
  
  names.w <- colnames(par.envcond.habitat.group$microhabaf$parvals)
  
  if(length(names.w)>0)
  {
    for (i in 1:nrow(par.envcond))  # loop over state variables
    {
      par.trait.cond[i,"fmicrohab"] <- 1
      
      # get areal fraction of microhab types
      
      af.i <- par.envcond.habitat.group$microhabaf$parvals[i,]
      if(sum(!is.na(af.i))>0)
      {
        if( round(sum(af.i,na.rm=TRUE),digits=5) > 1) warning(paste("sum of areal fractions of microhabitat types > 1:",
                                                                    sum(af.i,na.rm=TRUE),"for",i))
        if( round(sum(af.i,na.rm=TRUE),digits=5) < 1) warning(paste("sum of areal fractions of microhabitat types < 1:",
                                                                    sum(af.i,na.rm=TRUE),"for",i))
        
        # get trait parameters (affinity scores) of statevariable i
        
        m.trait.untrans <- par.taxaprop.traits$microhabtolval$parval[i,]    
        
        if( sum(!is.na(m.trait.untrans)) > 0  ) # some affinity scores are available
        { 
          # exponential transformation of affinity scores:
          
          m.trait.trans <- exp.transform(m.trait.untrans,
                                         intercept=par["fmicrohab_intercept"],
                                         curv=par["fmicrohab_curv"])
          
          # calculate weighted average of transformed affinity scores weighted by areal fraction of microhab type
          
          par.trait.cond[i,"fmicrohab"] <- sum( af.i * m.trait.trans,na.rm=TRUE)
        }
        
      } else { 
        
        if(grepl("_Invertebrates",rownames(par.envcond)[i]))
        {
          warning(paste("areal fractions of microhabitat types not available for ",rownames(par.envcond)[i] ))
        }
      }
    } # loop over state variables
  } else warning( "limitation by microhabitat conditions not considered, no areal fractions provided")  
  
  return(par.trait.cond)
}

