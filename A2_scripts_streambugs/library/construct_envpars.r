construct.envpars <- function (Reaches,Habitats,env.data)
{
  envpar <- list()
  
  if(length(Reaches)!=length(Habitats)) stop("No of Reaches",length(Reaches),"not equal to",
                                             "No of Habitats", length(Habitats))
  
  for (i in 1:length(Reaches))
  {
      # i <- 19
    reach <- Reaches[i]
    habitat <- Habitats[i]
    
    # XXX ecr 19.02.24: changed condition to get exact match of ReachID, 
    # because previous conditions led to multiple sites matching ReachID, and then problem in parameter sampling
    rind <- intersect(which(reach == env.data[,"ReachID"]), grep(habitat,env.data[,"Habitat"])) 
    # rind <- intersect(grep(reach,env.data[,"ReachID"]), grep(habitat,env.data[,"Habitat"]))
    
    envpar[paste(reach,"w",sep="_")]               <- get.par.dist.def(parname="w",env.data,rind)  # m
    envpar[paste(reach,"L",sep="_")]               <- get.par.dist.def(parname="L",env.data,rind)  # m
    
    # mean Temperature: 
    
    #envpar[paste(reach,habitat,"TdegC",sep="_")]   <- get.par.dist.def(parname="Temp_average",env.data,rind) #degC
    
    temp.C <- get.par.dist.def(parname="Temp_average",env.data,rind)     # degC
    
    if(length(temp.C[[1]])>1) 
    {
      temp.K <- temp.C
      
      dist <- temp.K[[1]][1]
      
      if(dist!="Lognormal"|dist!="lognormal"|dist!="Normal"|dist!="normal"|
         dist=="LognormalTrunc"|dist=="lognormaltrunc"|dist=="NormalTrunc"|dist=="normaltrunc")
      {
        temp.K[[1]][2] <- as.numeric(temp.C[[1]][2]) + 273.15 
        if ( length(temp.K[[1]]) > 3 )  # convert min and max for truncated distributions  XXX PR
        {
          temp.K[[1]][4] <- as.numeric(temp.C[[1]][4]) + 273.15 
          temp.K[[1]][5] <- as.numeric(temp.C[[1]][5]) + 273.15 
        }
        # K
      } else {
        
        if(dist=="Uniform"|dist=="uniform") 
        {
          temp.K[[1]][2] <- as.numeric(temp.C[[1]][2]) + 273.15                         # K
          temp.K[[1]][3] <- as.numeric(temp.C[[1]][3]) + 273.15                         # K
        } else {
          
          warning("Problem converting temperature from Celsius to Kelvin: converstion of distribution",dist,"not yet implemented in construct_envpars.r")
        }  
      }
      
      envpar[paste(reach,habitat,"T",sep="_")] <- temp.K 
      
    } else {
      envpar[paste(reach,habitat,"T",sep="_")] <- as.numeric(temp.C[[1]]) + 273.15 # K
    }
    
    
    envpar[paste(reach,habitat,"I0",sep="_")]      <- get.par.dist.def(parname="I0",env.data,rind)    # W/m2
    envpar[paste(reach,habitat,"fshade",sep="_")]  <- get.par.dist.def(parname="shade",env.data,rind) # -
    envpar[paste(reach,habitat,"CP",sep="_")]      <- get.par.dist.def(parname="C_P",env.data,rind)   # mgP/L
    envpar[paste(reach,habitat,"CN",sep="_")]      <- get.par.dist.def(parname="C_N",env.data,rind)   # mgN/L
    
    envpar[paste(reach,habitat,"CSusPOM",sep="_")] <- get.par.dist.def(parname="C_SusPOM",env.data,rind) #   gDM/m3  
    
    temp.max.C <- get.par.dist.def(parname="tempmaxC",env.data,rind)     # degC
    
    if(length(temp.max.C[[1]])>1)  # rather
    {
      temp.max.K <- temp.max.C
      
      dist <- temp.max.K[[1]][1]
      
      if(dist!="Lognormal"|dist!="lognormal"|dist!="Normal"|dist!="normal"|
         dist=="LognormalTrunc"|dist=="lognormaltrunc"|dist=="NormalTrunc"|dist=="normaltrunc")
      {
        temp.max.K[[1]][2] <- as.numeric(temp.max.C[[1]][2]) + 273.15 
        if ( length(temp.max.K[[1]]) > 3 )  # convert min and max for truncated distributions  XXX PR
        {
          temp.max.K[[1]][4] <- as.numeric(temp.max.C[[1]][4]) + 273.15 
          temp.max.K[[1]][5] <- as.numeric(temp.max.C[[1]][5]) + 273.15 
        }
        # K
      } else {
        
        if(dist=="Uniform"|dist=="uniform") 
        {
          temp.max.K[[1]][2] <- as.numeric(temp.max.C[[1]][2]) + 273.15                         # K
          temp.max.K[[1]][3] <- as.numeric(temp.max.C[[1]][3]) + 273.15                         # K
        } else {
          
          warning("converstion of distribution",dist,"not yet implemented in construct_envpars.r")
        }  
      }
      
      envpar[paste(reach,habitat,"tempmaxK",sep="_")] <- temp.max.K 
      
    } else {
      envpar[paste(reach,habitat,"tempmaxK",sep="_")] <- as.numeric(temp.max.C[[1]]) + 273.15 # K
    }
    
    envpar[paste(reach,habitat,"currentms",sep="_")]      <- get.par.dist.def(parname="currentms",env.data,rind)      # m/s
    envpar[paste(reach,habitat,"orgmicropollTU",sep="_")] <- get.par.dist.def(parname="orgmicropollTU",env.data,rind) # TU
    envpar[paste(reach,habitat,"saprowqclass",sep="_")]   <- get.par.dist.def(parname="saprowqclass",env.data,rind)   # wqclass 
    
    cind <- grep("microhabaf_type",colnames(env.data))
    for(j in 1:length(cind))
    {
      envpar[paste(reach,habitat,colnames(env.data)[cind[j]],sep="_")]   <- get.par.dist.def(parname=colnames(env.data)[cind[j]],env.data,rind)  # microhabitat types
    }
    
    envpar[paste("CPOM",reach,habitat,"Input",sep="_")] <- get.par.dist.def(parname="Lit_Inp",env.data,rind) # gDM/m2/a
    
  }
  
  return(envpar)
  
}