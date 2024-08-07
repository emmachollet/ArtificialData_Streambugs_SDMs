get.par.dist.def <- function (parname, env.data,rind)
  
  #  parname="Temp_average"
  #  parname="current" 
  #  parname="orgpollut"
  #  get.par.dist.def(parname, env.data,rind)
  
{
  cind<-NULL
  cnames <- strsplit(colnames(env.data),split="__")
  
  for ( i in 1:length(cnames))
  {
    if(cnames[[i]][1]==parname) cind<-c(cind,i)
  }
    
  distdef <- list()
  dist <- NA
  
  if(length(cind)==1)
  {
    if(parname==colnames(env.data)[cind])
    {
      mean <- NULL
      mean <- env.data[rind,cind]
      if(is.null(mean)) {mean <- NA; warning("mean for ",parname," with ",dist," distr. not found")}
      distdef <- list(c(mean))       
    }
    
  } else
  {
    
    if(length(cind)>1)
    {
      
      dist <- env.data[rind, paste(parname,"dist",sep="__")]
      
      if(!is.na(dist))
      {
        
        if(dist=="Uniform"|dist=="uniform") 
        {
          min <- env.data[rind, paste(parname,"min",sep="__")]
          max <- env.data[rind, paste(parname,"max",sep="__")]
          
          if(is.null(min)) {min <- NA; warning("min for ",parname," with uniform distr not found")}
          if(is.null(max)) {max <- NA; warning("max for ",parname," with uniform distr not found")}
          
          distdef <- list(c(dist,min,max))
          
        } else
        {
          if(dist=="Lognormal"|dist=="lognormal"|dist=="Normal"|dist=="normal")
          {
            mean <- env.data[rind, paste(parname,"mean",sep="__")]
            
            abssd <- env.data[rind, paste(parname,"abssd",sep="__")]
            
            if (is.null(abssd))
            {
              relsd <- env.data[rind, paste(parname,"relsd",sep="__")]
              if (!is.null(relsd) & !is.null(mean))
              {
                abssd <- as.numeric(relsd)*as.numeric(mean)
              }
            }
            
            if(is.null(mean))  {mean  <- NA; warning("mean for ",parname, " with ",dist," distr. not found")}
            if(is.null(abssd)) {abssd <- NA; warning("sd for ",parname," with ",dist," distr. not found")}
            
            distdef <- list(c(dist,mean,abssd))
            
          } else
          {
            if(dist=="LognormalTrunc"|dist=="lognormaltrunc"|dist=="NormalTrunc"|dist=="normaltrunc")
            {
              mean <- env.data[rind, paste(parname,"mean",sep="__")]
              
              abssd <- env.data[rind, paste(parname,"abssd",sep="__")]
              
              if (is.null(abssd))
              {
                relsd <- env.data[rind, paste(parname,"relsd",sep="__")]
                if (!is.null(relsd) & !is.null(mean))
                {
                  abssd <- as.numeric(relsd)*as.numeric(mean)
                }
              }
              
              min <- env.data[rind, paste(parname,"min",sep="__")]
              max <- env.data[rind, paste(parname,"max",sep="__")]
              
              if(is.null(mean))  {mean  <- NA; warning("mean for ",parname," with ",dist," distr. not found")} 
              if(is.null(abssd)) {abssd <- NA; warning("sd for ",parname," with ",dist," distr. not found")}
              if(is.null(min))   {min   <- 0} #; warning("min for ",parname," with ",dist," distr. not found, set to 0")}
              if(is.null(max))   {max   <- 1e30} #; warning("max for ",parname," with ",dist," distr. not found, set to 1e30")}
              
              distdef <- list(c(dist,mean,abssd,min,max))
            } else
            {
              if(dist=="Discrete"|dist=="discrete")
              {
                nc <- (length(cind)-1)/2
                class <- NULL
                probs <- NULL
                for (i in 1:nc)
                {
                  class <- c(class,env.data[rind,paste(parname,"__class",i,sep="")])
                  probs <- c(probs,env.data[rind,paste(parname,"__prob",i,sep="")])
                } 
                if(is.null(class)) {class <- NA; warning("classes for ",parname," with ",dist," distr. not found")}
                if(is.null(probs)) {probs <- NA; warning("probabilities for ",parname," with ",dist," distr. not found")}
                
                distdef <- list(c(dist,class,probs))
              }
            }
          } 
        }        
      }    else warning("dist not given but more than one parameter for ",parname," exists")          
    }
  }
  
  return(distdef)
}