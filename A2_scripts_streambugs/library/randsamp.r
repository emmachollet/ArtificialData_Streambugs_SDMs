# randsamp Nele Schuwirth Version 17.1.2011

# Load libraries:
# ---------------

if ( !require(msm) )
{
   install.packages("msm")
   library(msm)
}

if ( !require(triangle) )
{
   install.packages("triangle")
   library(triangle)
}

#sampling from probability distributions
#---------------------------------------
# arguments:
   # n:               sample size
   # dist:            distribution type: "Normal" or "Lognormal" or...
   # dist.par:        distribution parameters:
   # Normal:          mean, sd
   # NormalTrunc:     mean, sd, min, max
   # Lognormal:       mean, sd             #! Achtung, anders als in R standardmässig!
   # LognormalTrunc:  mean, sd, min, max   #! Achtung, anders als in R standardmässig!
   # Triangle:        min, mode, max
   # Beta:            shape1, shape2, min, max
   # Unif:            min, max
   # Exponential:     mean                 (the function rexp() uses the rate=1/mean)
# output:   samp:     random sample from the defined distribution

randsamp<-function(n=1000,dist="Normal",dist.par=c(0,1),file=NA)
{
 if(dist == "Unif" )
  {
   samp <- runif(n,min=dist.par[1],max=dist.par[2])
  }
 else
  {if (dist == "Normal" )
  {
   samp <- rnorm(n,mean=dist.par[1],sd=dist.par[2])
  }
 else
  {if(dist == "NormalTrunc" )
  {
   samp <- rtnorm(n, mean=dist.par[1],sd=dist.par[2], lower=dist.par[3],
                  upper=dist.par[4])
  }
 else
  {if(dist == "Lognormal" )
  {
    mean  <-dist.par[1]
    sd    <-dist.par[2]
    sdlog <- sqrt(log(1+sd^2/mean^2))
    meanlog <- log(mean) - 0.5*sdlog^2
    samp <- rlnorm(n,meanlog=meanlog,sdlog=sdlog)
  }
 else
 {if(dist == "LognormalTrunc" )
  {
    mean  <-dist.par[1]
    sd    <-dist.par[2]
    sdlog <- sqrt(log(1+sd^2/mean^2))
    meanlog <- log(mean) - 0.5*sdlog^2
    
    samp <- rlnorm(n,meanlog=meanlog,sdlog=sdlog)
    
    if(min(samp) < dist.par[3] | max(samp) > dist.par[4]   )
    { 
      repeat
      {ind1 <- samp >= dist.par[3] 
       ind2 <- samp <= dist.par[4]
       ind <- ind1&ind2
       samp <- samp[ind]
       samp <- c(samp,rlnorm(sum(ind==F),meanlog=meanlog,sdlog=sdlog))     
       if( min(samp) >= dist.par[3] & max(samp) <= dist.par[4] ) {break}
      } 
    }   
  }
 else
 {if(dist == "Beta" )
  {
   samp <- rbeta(n,shape1=dist.par[1],shape2=dist.par[2])
  }
 else
 {if(dist == "Triangle")
  {
   samp <- rtriangle(n, a=dist.par[1], b=dist.par[3], c=dist.par[2])
  }
  else
  {if(dist == "Exponential")
   {  
    samp <- rexp(n,rate=1/dist.par[1])
   }
   else    {
         stop("randsamp: unknown distribution type")
      }
  }}}}}}}

 return(samp)
}



#samp <- randsamp(n=1000,dist="LognormalTrunc",dist.par=c(4,1.4,3,5),file=NA)
