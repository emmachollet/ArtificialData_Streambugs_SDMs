################################################################################
#                                                                              #
# Systems analysis R library                                                   #
# ==========================                                                   #
#                                                                              #
# Peter Reichert, Eawag, reichert@eawag.ch, last modification Sept 23, 2012    #  
#  # with small mod. in line 940 on Mar 04, 2013 NIS                           #        
#  # modification in plot.chains, postsamp can be a list Jun 6,2013            #  
#                                                                              #
################################################################################

# Overview of functions provided by this library
# (see individual descriptions of functions for more details):
#
# Data handling
# -------------
# sysanal.decode                  Decode a result code used by plotting funct.
# sysanal.read.distdef            Read distribution definition
#
# Plotting:
# ---------
# sysanal.plot.res                Plot results of a simulation
# sysanal.plot.vars               Plot several sets of results
# sysanal.contourpdf              Contour a 2d probability density function
# sysanal.plot.chains             Plot components of a Markov chain
# sysanal.plot.margs              Plot marginals of a sample
#
# Sampling:
# ---------
# sysanal.randsamp                Sample from a mv. normal or lognormal dist.
#                                 or from independent 1d distributions
# sysanal.markovchain.metropolis  Calculate a Metropolis Markov chain of a pdf
# sysanal.jump.dist               Calculate a jump dist. for a Metropolis MC
#
# Model formulation:
# ------------------
# sysanal.calcpdf                 Calculate density of univariate distributions
# sysanal.calcpdf_mv              Calculate density of multivariate distrib.
# sysanal.loglikeli               Calculate log likelihood of simple model
# sysanal.logposterior            Calculate log posterior of simple model
#
# Sensitivity analysis:
# ---------------------
# sysanal.sens.loc                Calculate local sensitivities of det. model
# sysanal.sens.var.samp           Calculate variance-based sens. of det. model
#
# Identifiability analysis:
# -------------------------
# sysanal.comb                    Calculate combinations of subsets of indices
# sysanal.collind                 Calculate collinearity index
# sysanal.ident                   Calculate identifiability measures
#
# Output transformation:
# ----------------------
# sysanal.boxcox                  Box-Cox transformation
# sysanal.boxcox.deriv            derivative of Box-Cox transformation
# sysanal.boxcox.inv              Inverse of Box-Cox transformation
#
# Multivariate regression:
# ------------------------
# sysanal.gnls                    Generalized nonlinear least squares regression
# sysanal.gnls.diag               Diagnostics for gnls
# sysanal.gnls.test               Statistical tests for gnls
# sysanal.gnls.predict            Prediction and confidende intervals for gnls
# sysanal.confint                 Estimation of confidence intervals
#
# Residual diagnostics:
# ---------------------
# sysanal.resid.diag              Residual diagnostics plots
# sysanal.resid.diag.boxcox       Residual diagnostics plots with Box-Cox trans.
#
# Miscellaneous:
# --------------
# sysanal.hessian                 Calculate the Hessian of a real-valued funct.
#
################################################################################

# sysanal.decode
# ==============

# purpose:
# decode a simple model layout code into variable names and value of input 
# variable
# (the function assumes that the code consists of the name of the output 
# variable and the value of the input variable separated by an underline 
# character; the variable name may still contain underline characters)

# arguments:
# rescode:     vector of result codes or vector of results with codes as 
#              names of the elements (the function assumes that the code 
#              consists of the name of the output variable and the value of the 
#              input variable separated by an underline character; the variable 
#              name may still contain underline characters)

# output:
# data frame with variables:
# var:         vector of character strings of variables
# val:         vector of corresponding values of input variable
# the row names of the data frame contain the original codes

sysanal.decode <- function(L)
{
  if ( is.numeric(L) ) codes <- names(L)
  else                 codes <- L
  var <- character(length(codes))
  val <- numeric(length(codes))
  for ( i in 1:length(codes))
  {
    s <- strsplit(as.character(codes[i]),split="_")
    n <- length(s[[1]])
    if ( n < 2 ) stop("error in sysanal.decode: illecal result code")
    var[i] <- paste(s[[1]][1:(n-1)],collapse="_")
    val[i] <- as.numeric(s[[1]][n])
  }
  L.decoded <- data.frame(var=var,val=val)
  L.decoded$var <- as.character(L.decoded$var)
  rownames(L.decoded) <- codes
  return(L.decoded)
}

################################################################################

sysanal.read.distdef <- function(file,sep="\t")
{
  distdef.tab <- read.delim(file,header=F,sep=sep)
  distdef.tab <- cbind(distdef.tab,rep(NA,nrow(distdef.tab)))
  distdef <- list()
  for ( i in 1:nrow(distdef.tab) )
  {
    max.ind <- min(match(NA,t(distdef.tab[i,])),match("",t(distdef.tab[i,])),na.rm=T)-1
    distdef[[i]] <- as.vector(t(distdef.tab[i,][2:max.ind]))
  }
  names(distdef) <- as.character(distdef.tab[,1])
  return(distdef)
}

################################################################################

sysanal.package <- function(package)
{
  if ( is.na(match(package,installed.packages()[,1])) )
  {
    print(paste("Do you agree to install the package \"",
                package,
                "\"?",
                sep=""))
    if ( menu(choices=c("yes","no")) == 1 )
    {
      install.packages(package)
    }
    else
    {
      stop(paste("Package \"",package,"\" not found",sep=""))
    }
  }
  library(package,character.only=TRUE)
}

################################################################################

# sysanal.plot.res
# ================

# purpose:
# plot results provided as a numeric vector with result codes

# arguments:
# res:         vector of results named by result codes 
#              (see sysanal.decode for an explanation of result codes)
# xlim:        optional limits of the x axis
# ylim:        optional limits of the y axis
# markers:     if TRUE plot markers instead of lines
# header:      optional header of the plot
# xlab:        optional label of the x axis
# ylab:        optional label of the y axis
# pos:         position of legend (only if more than one variable encoded)

# output:
# plot of all variables as a function of the independent variable
# (note that variable names and values of the independent variable are
# encoded in the names of the components of the result vector)

sysanal.plot.res <- function(res,xlim=NA,ylim=NA,
                             markers=F,header="",
                             xlab="",ylab="",pos="topright")
{
  codes <- sysanal.decode(res)
  varnames <- unique(codes$var)
  if ( is.na(xlim[1]) )      xlim <- range(codes$val)
  if ( is.na(ylim[1]) )      ylim <- range(res)
  if ( nchar(ylab[1]) == 0 ) ylab <- paste(varnames,collapse=", ")
  plot(numeric(0),numeric(0),type="n",xlim=xlim,ylim=ylim,
       xlab=xlab,ylab=ylab,main=header)
  if ( markers )
  {
    for ( i in 1:length(varnames) )
    {
      ind <- codes$var == varnames[i]
      points(codes$val[ind],res[ind],pch=i)
    }
    if ( length(varnames) > 1 )
    {
      legend(x=pos,legend=varnames,pch=1:length(varnames))
    }
  }
  else
  {
    for ( i in 1:length(varnames) )
    {
      ind <- codes$var == varnames[i]
      lines(codes$val[ind],res[ind],lty=i)
    }
    if ( length(varnames) > 1 )
    {
      legend(x=pos,legend=varnames,lty=1:length(varnames))
    }
  }
}

################################################################################

# sysanal.plot.vars
# =================

# purpose:
# plot variables provided as a data frame or matrix with result codes given 
# by the row names

# arguments:
# vars:        matrix or data frame with variables and result codes as row names 
#              (see sysanal.decode for an explanation of result codes)
# ncol:        optional number of columns of sub-panels of the plot
# mar:         optional specification of margins in the form 
#              c(bottom,left,top,right)
# ylim:        optional named (by variable name) list of limits of the y axes
# markers:     if TRUE plot markers instead of lines
# header:      optional named (by variable name) list of headers of the plots
# xlab:        optional label of the x axis
# ylab:        optional label of the y axis
# pos:         position of legend (only if more than one variable encoded)

# output:
# plot of all variables as a function of the independent variable
# (note that variable names and values of the independent variable are
# encoded in the names of the components of the result vector)

sysanal.plot.vars <- function(vars,ncol=NA,mar=NA,
                              ylim=list(),markers=F,
                              headers=list(),xlab="",ylab="",pos="topright")
{
  nvar <- ncol(vars)
  if ( is.na(ncol) ) nc <- floor(sqrt(nvar))
  nr <- ceiling(nvar/nc)
  marg <- mar
  if ( is.na(marg[1]) ) marg <- c(4.5,4.0,2.5,1.0) # c(bottom, left, top, right)
  par.def <- par(no.readonly=T)
  par(mfrow=c(nr,nc),xaxs="i",yaxs="i",mar=marg) 
  for ( i in 1:nvar )
  {
    name <- colnames(vars)[i]
    ylim.i <- c(min(vars[,i]),max(vars[,i]))
    if ( ylim.i[1] == ylim.i[2] ) 
    {
      ylim.i[1] <- 0.5*ylim.i[1]
      ylim.i[2] <- 1.5*ylim.i[2]
    }
    if ( length(ylim) > 0 )
    {
      ind <- match(name,names(ylim))
      if ( !is.na(ind) )
      {
        ylim.i <- ylim[[ind]]
      }
    }
    header.i <- name
    if ( length(headers) > 0 )
    {
      ind <- match(name,names(headers))
{
        if ( !is.na(ind) )
        {
          header.i <- headers[[ind]]
        }
      }
    }
    res <- as.numeric(vars[,i])
    names(res) <- rownames(vars)
    sysanal.plot.res(res,
                     header=header.i,markers=markers,
                     xlab=xlab,ylab=ylab,pos=pos,
                     xlim=NA,ylim=ylim.i)
  }
  par(par.def)
}

################################################################################

# Contour the probability density function of a bivariate normal or 
# lognormal distribution

sysanal.contourpdf <- function(calcpdf_mv,norm=T,xlim=c(-3,3),ylim=c(-3,3),
                               levels=c(0.05,0.5,0.95),res=20,lty="solid",
                               ...)
{
  # -----------------------------------------------------------------------
  # This function plots contour lines of a normalized or unnormalized
  # bivariate probability density function. If the function is not 
  # normalized, the integral over the domain specified by xlim and ylim is  
  # used for normalization (this leads to incorrect results if this domain 
  # does not contain most of the distribution).
  #
  # Arguments:
  # calcpdf_mv: function to calculate the log pdf values at a set of
  #             locations specified by its first argument. Further 
  #             arguments specified under ... will be passed to this 
  #             function
  # norm:       TRUE if the probability density is normalized,
  #             FALSE if integration over the domain given by xlim and ylim
  #             should be used for normalization
  # xlim:       bounds of the integration range for the first variable
  # ylim:       bounds of the integration range for the second variable 
  # levels:     vector of probabilities to be contained in the contour line
  # res:        resolution of grid used to countour
  #             (number of points in each dimension)
  # lty:        line type of contour lines
  #
  # Return Value:
  # integral of the probability density over the given range at the given
  # resolution
  #
  #                                        Peter Reichert    Feb.  02, 2005
  #                                        last modification March 31, 2008
  # -----------------------------------------------------------------------
  
  dx.grid <- (xlim[2]-xlim[1])/res
  dy.grid <- (ylim[2]-ylim[1])/res
  x.grid <- seq(xlim[1]+dx.grid/2,xlim[2]-dx.grid/2,len=res)
  y.grid <- seq(ylim[1]+dy.grid/2,ylim[2]-dy.grid/2,len=res)
  xarray <- vector()
  for ( i in 1:res ) xarray = c(xarray,rep(x.grid[i],res))
  yarray <- rep(y.grid,res)
  z.sample <- cbind(xarray,yarray)
  logpdf <- calcpdf_mv(z.sample,...)
  if ( norm == F) logpdf <- logpdf-max(logpdf,na.rm=T)
  pdf <- ifelse(is.na(logpdf),0,exp(logpdf))
  pdf.sort <- sort(pdf,decreasing=T)
  integral.cum <- pdf.sort*dx.grid*dy.grid
  for ( i in 2:(res*res) ) 
  {
    integral.cum[i] <- integral.cum[i-1]+integral.cum[i]
  }
  integral <- integral.cum[res*res]
  if ( norm == F) integral.cum <- integral.cum/integral
  logpdflevels <- vector()
  index <- 1
  levelssort <- sort(levels)
  for ( i in 1:(res*res) )
  {
    if ( integral.cum[i] > levelssort[index] )
    {
      logpdflevels[index] <- log(pdf.sort[i])
      index <- index + 1
      if ( index > length(levelssort) ) break
    }
  }
  logpdf <- matrix(logpdf,nrow=res,ncol=res,byrow=TRUE)
  contour(x=x.grid,y=y.grid,z=logpdf,levels=logpdflevels,
          add=TRUE,drawlabels=FALSE,lty=lty)
  return(integral)
}

################################################################################

# function to plot markov chains:
# -------------------------------

sysanal.plot.chains <- function(postsamp=list(),ncol=NA,nrow=NA,mar=NA,
                                ylim=list(),
                                titles=list(),xlab="chain index",ylab=list())
{
  if(!is.list(postsamp)) postsamp <- list(postsamp)
  
  nvar <- ncol(postsamp[[1]])
  
  if ( !is.na(ncol) ) nc=ncol else nc <- floor(sqrt(nvar))
  if ( !is.na(nrow) ) nr=nrow else nr <- ceiling(nvar/nc)
  
  marg <- mar
  if ( is.na(marg[1]) ) marg <- c(2.0,2.0,2.5,0.5) # c(bottom, left, top, right)
  par.def <- par(no.readonly=T)
  par(mfrow=c(nr,nc),xaxs="i",yaxs="i",mar=marg) 
  
  for ( i in 1:nvar )
  {
    name <- colnames(postsamp[[1]])[i]
    data <- postsamp[[1]][,i]
    
    for (j in 1:length(postsamp) )
    {
      data <- c(data,postsamp[[j]][,i])
    }
    
    ylim.i <- c(min(data),max(data))
    if ( ylim.i[1] == ylim.i[2] ) 
    {
      ylim.i[1] <- 0.5*ylim.i[1]
      ylim.i[2] <- 1.5*ylim.i[2]
    }
    if ( length(ylim) > 0 )
    {
      ind <- match(name,names(ylim))
      if ( !is.na(ind) )
      {
        ylim.i <- ylim[[ind]]
      }
    }
    title.i <- name
    if ( length(titles) > 0 )
    {
      ind <- match(name,names(titles))
      {
        if ( !is.na(ind) )
        {
          title.i <- titles[[ind]]
        }
      }
    }
    ylab.i <- name
    if ( length(ylab) > 0 )
    {
      ind <- match(name,names(ylab))
      if ( !is.na(ind) )
      {
        ylab.i <- ylab[[ind]]
      }
    }
    plot(postsamp[[1]][,i],ylim=ylim.i,type="l",main=title.i,xlab=xlab,ylab=ylab.i)
    if(length(postsamp)>1)
    {
      for (j in 2:length(postsamp))
      {
        lines(postsamp[[j]][,i],col=j)
      }
    }
  }
  par(par.def)
}

################################################################################

# function to plot prior and posterior marginals:
# -----------------------------------------------

sysanal.plot.margs <- function(postsamp,pridist=list(),val=NA,ncol=NA,nrow=NA,mar=NA,
                               xlim=list(),ymax=list(),
                               titles=list(),xlab=list(),ylab=list(),
                               adjust=1,lty=NA,col.post=NA,col.pri=NA)
{
  # transform samples to a list of all samples to plot the marginals of:
  
  postsamp.list <- list()
  if ( is.data.frame(postsamp) )
  {
    postsamp.list[[1]] <- postsamp
  } else
  {
    if ( is.matrix(postsamp) )
    {
      postsamp.list[[1]] <- postsamp
    } else
    {
      if ( is.list(postsamp) )  # be careful: a data frame is a list!
      {
        postsamp.list <- postsamp
      } else
      {
        stop("sysanal.plot.margs: postsamp is of illegal type")
      }
    }
  }
  nsamp <- length(postsamp.list)
  
  # transform prior definitions to a list of priors to plot:
  
  pridist.list <- list()
  if ( length(pridist) > 0 )
  {
    if ( length(names(pridist)) != length(pridist) )
    {
      pridist.list <- pridist
    }     else
    {
      pridist.list[[1]] <- pridist
    }
  }
  npri <- length(pridist.list)
  
  # get all variable names of the sample: 
  
  var <- colnames(postsamp.list[[1]])
  if ( nsamp > 1 )
  {
    for ( j in 2:nsamp ) 
    {
      var <- c(var,colnames(postsamp.list[[j]]))
    }
    var <- unique(var)
  }
  nvar <- length(var)
  
  # define layout of plot panels:
  
  if ( !is.na(ncol) ) nc <- ncol   else                nc <- floor(sqrt(nvar))
  if ( !is.na(nrow) ) nr <- nrow   else                nr <- ceiling(nvar/nc)
  marg <- mar
  if ( is.na(marg[1]) ) marg <- c(2.0,2.0,2.5,0.5) # c(bottom, left, top, right)
  
  # define line types:
  
  lty.loc <- lty
  if ( is.na(lty.loc[1]) ) lty.loc <- 1
  if ( length(lty.loc) < nsamp+npri ) lty.loc <- 1:(nsamp+npri)
  
  # plot marginals of samples and priors:
  
  par.def <- par(no.readonly=T)
  par(mfrow=c(nr,nc),xaxs="i",yaxs="i",mar=marg) # c(bottom, left, top, right)
  for ( i in 1:nvar )
  {
    name <- var[i]
    data.min    <- NA
    data.max    <- NA
    density.max <- NA
    marg.samp <- as.list(rep(NA,nsamp))
    for ( j in 1:nsamp )
    {
      ind <- match(var[i],colnames(postsamp.list[[j]]))
      if ( !is.na(ind) )
      {
        data.min       <- min(c(data.min,postsamp.list[[j]][,ind]),na.rm=TRUE)
        data.max       <- max(c(data.max,postsamp.list[[j]][,ind]),na.rm=TRUE)
        marg.samp[[j]] <- density(postsamp.list[[j]][,ind],adjust=adjust)
        density.max    <- max(c(density.max,marg.samp[[j]]$y),na.rm=TRUE)
      }
    }
    
    xlim.i <- c(data.min,data.max) # xxx ev multiply with 0.9 and 1.1
    if ( xlim.i[1] == xlim.i[2] ) 
    {
      xlim.i[1] <- 0.5*xlim.i[1]  
      xlim.i[2] <- 1.5*xlim.i[2]  
    }
    if ( length(xlim) > 0 )
    {
      ind <- match(name,names(xlim))
      if ( !is.na(ind) )
      {
        xlim.i <- xlim[[ind]]
      }
    }
    
    if ( is.na(density.max) ) ylim.i <- c(0,1) else ylim.i <- c(0,1.1*density.max)
    if ( length(ymax) > 0 )
    {
      ind <- match(name,names(ymax))
      if ( !is.na(ind) )
      {
        ylim.i <- c(0,ymax[[ind]])
      }
    }
    
    title.i <- name
    if ( length(titles) > 0 )
    {
      ind <- match(name,names(titles))
      {
        if ( !is.na(ind) )
        {
          title.i <- titles[[ind]]
        }
      }
    }
    
    xlab.i <- name
    if ( length(xlab) > 0 )
    {
      ind <- match(name,names(xlab))
      if ( !is.na(ind) )
      {
        xlab.i <- xlab[[ind]]
      }
    }
    
    ylab.i <- "f"
    if ( length(ylab) > 0 )
    {
      ind <- match(name,names(ylab))
      if ( !is.na(ind) )
      {
        ylab.i <- ylab[[ind]]
      }
    }
    
    marg.pri <- list()
    if ( npri > 0 )
    {
      x <- seq(xlim.i[1],xlim.i[2],length=101)
      y <- numeric(0)
      for ( j in 1:npri )
      {
        marg.pri[[j]] <- list()
        if ( length(pridist.list[[j]]) > 0 )
        {
          ind <- sysanal.multmatch(name,names(pridist.list[[j]]))
          if ( !is.na(ind[1]) )
          {
            for ( k in 1:length(ind) )
            {
              y <- sysanal.calcpdf(x,pridist.list[[j]][[ind[k]]])
            }
          }
        }
        marg.pri[[j]]$x <- x
        marg.pri[[j]]$y <- y
      }
    }
    
    # plot marginals of single variable:
    
    # plot frame:
    
    plot(numeric(0),numeric(0),main=title.i,xlab=xlab.i,ylab=ylab.i,
         xlim=xlim.i,ylim=ylim.i,type="n")
    
    # plot areas:
    
    if ( npri > 0 & !is.na(col.pri) )
    {
      for ( j in 1:npri )
      {
        if ( !is.na(marg.pri[[j]])[[1]] )
        {
          n <- length(marg.pri[[j]]$x)
          polygon(c(marg.pri[[j]]$x[1],marg.pri[[j]]$x,marg.pri[[j]]$x[n],marg.pri[[j]]$x[1]),
                  c(0,marg.pri[[j]]$y,0,0),
                  border=NA,col=col.pri)
        }
      }
    }
    if ( !is.na(col.post) )
    {
      for ( j in 1:nsamp )
      {
        if ( !is.na(marg.samp[[j]])[[1]] )
        {
          n <- length(marg.samp[[j]]$x)
          polygon(c(marg.samp[[j]]$x[1],marg.samp[[j]]$x,marg.samp[[j]]$x[n],marg.samp[[j]]$x[1]),
                  c(0,marg.samp[[j]]$y,0,0),
                  border=NA,col=col.post)
        }
      }
    }
    
    # plot lines:
    
    abline(h=ylim.i[1])
    abline(h=ylim.i[2])
    abline(v=xlim.i[1])
    abline(v=xlim.i[2])
    if ( !is.na(val[var[i]]) ) abline(v=val[var[i]])
    if ( npri > 0 )
    {
      for ( j in 1:npri )
      {
        if ( !is.na(marg.pri[[j]])[[1]] )
        {
          lines(marg.pri[[j]],lty=lty.loc[nsamp+j])
        }
      }
    }
    for ( j in 1:nsamp )
    {
      if ( !is.na(marg.samp[[j]])[[1]] )
      {
        lines(marg.samp[[j]],lty=lty.loc[j])
      }
    }
  }
  par(par.def)
}


sysanal.multmatch <- function(x,table,incomparables=NULL)
{
  table.loc <- table
  inds <- numeric(0)
  while (TRUE)
  {
    ind <- match(x,table.loc,incomparables)
    if ( is.na(ind) )
    {
      if ( length(inds) == 0 ) return(ind)
      else                     return(inds)
    }
    inds <- c(inds,ind)
    table.loc[ind] <- paste(x,"_",sep="")
  }  
}

################################################################################

# generate a random sample from a multivariate normal or lognormal distribution 
# or from a product of independent 1d marginals

sysanal.randsamp <- function(sampsize=1,dist="Normal",mean=0,sd=1,cor=NA,
                             distdef=NA,file=NA)
{
  # -----------------------------------------------------------------------
  # This function generates a random sample from a multivariate
  # normal or lognormal distribution.
  # This function is a simplified version of the program "randsamp"
  # available as part of the package UNCSIM at http://www.uncsim.eawag.ch
  #
  # Arguments:
  # sampsize:   sample size
  # mean:       vector of means
  # sd:         vector of standard deviations
  # cor:        correlation matrix of the distribution
  # dist:       distribution type: "Normal", "Lognormal" or "Indep".
  # distdef:    definition of 1d marginals for "Indep" instead of (mean,sd,cor)
  #
  # Return Value:
  # List of:
  # mean:       vector of means
  # sd:         vector of standard deviations
  # corr:       correlation matrix of the distribution
  # sampsize:   sample size
  # sample:     data frame of parameter samples (each row corresponds 
  #             to a draw)
  # logpdf:     vector of values of log probability density at the sample
  #             points
  #
  #                                        Peter Reichert    Dec.  29, 2003
  #                                        last modification Oct.  29, 2011
  # -----------------------------------------------------------------------
  
  R <- 0
  if ( dist == "normal"    | dist == "Normal" | 
         dist == "lognormal" | dist == "Lognormal")
  {
    # consistency checks and initializations:
    
    numpar <- length(mean)
    if ( length(sd) != numpar )
    {
      stop("sysanal.randsamp: mean and sd do not have the same length")
    }
    R <- diag(rep(1,numpar))
    if ( is.matrix(cor) ) R <- cor
    if ( nrow(R) != numpar || ncol(R) != numpar )
    {
      stop("sysanal.randsamp: illegal dimension of correlation matrix")
    }
    
    # calculate sample from multivariate uniform distribution:
    samp <- runif(sampsize*numpar)
    dim(samp) <- c(sampsize,numpar)
    
    # transform sample to multivariate normal or lognormal and calculate
    # logarithm of probability density function:
    logpdf <- numeric(sampsize)
    if ( dist == "normal" | dist == "Normal" )
    {
      # calculate transformation matrix and transform the sample:
      sigma <- diag(sd) %*% R %*% diag(sd)
      sigma.inv = solve(sigma)
      det.sigma = det(sigma)
      A <- t(chol(sigma))
      for ( i in 1:sampsize )
      {
        samp[i,]  <- A %*% qnorm(samp[i,]) + mean
        v <- samp[i,]-mean
        v <- matrix(v,nrow=numpar,ncol=1) # ensure that v is a column vector
        logpdf[i] <- -numpar/2*log(2*pi) - 0.5*log(det.sigma) -
          0.5 * t(v) %*% sigma.inv %*% (v)
      }
    }
    else # dist == "lognormal" | dist == "Lognormal"
    {
      # parameters of the log of the variable, calculate transformation 
      # matrix and transform the sample:
      sdlog    <- sqrt(log(1+sd*sd/(mean*mean)))
      meanlog  <- log(mean) - sdlog*sdlog/2
      if ( numpar > 1 )
      {
        ln.sigma <- log( 1 + diag(sqrt(exp(sdlog*sdlog)-1)) %*%
                           R %*% diag(sqrt(exp(sdlog*sdlog)-1)) )
      }
      else
      {
        ln.sigma <- as.matrix( log( 1 + sqrt(exp(sdlog*sdlog)-1)^2 ) )
      }
      ln.sigma.inv = solve(ln.sigma)
      det.ln.sigma = det(ln.sigma)
      ln.A <- t(chol(ln.sigma))
      for ( i in 1:sampsize )
      {
        log.samp.i <- ln.A %*% qnorm(samp[i,]) + meanlog
        samp[i,]  <- exp(log.samp.i)
        v <- log.samp.i-meanlog
        v <- matrix(v,nrow=numpar,ncol=1) # ensure that v is a column vector
        logpdf[i] <- -numpar/2*log(2*pi) - 0.5*log(det.ln.sigma) - 
          log(prod(samp[i,])) - 0.5 * t(v) %*% ln.sigma.inv %*% v
      }
    }
    
    # collect results:
    colnames(samp) <- names(mean)
    samp <- data.frame(samp)
    res <- list(
      mean       = mean,
      sd         = sd,
      cor        = R,
      sampsize   = sampsize,
      sample     = samp,
      logsamppdf = logpdf
    )
  }
  else   # dist != normal,Normal,lognormal,Lognormal      
  {  
    if ( dist == "indep" | dist == "Indep" )
    {
      logpdf <- 0
      samp <- matrix(NA,nrow=sampsize,ncol=length(distdef))
      colnames(samp) <- names(distdef)
      samp <- data.frame(samp)
      for ( j in 1:ncol(samp) )
      {
        distpar <- distdef[[j]]
        dist.found <- F
        if ( !dist.found )
        {
          if ( distpar[1] == "Uniform" | 
                 distpar[1] == "uniform" )
          {
            # uniform distribution; parameters are min and max
            min <- as.numeric(distpar[2])
            max <- as.numeric(distpar[3])
            samp[,j] <- runif(sampsize,min=min,max=max)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Normal" | 
                 distpar[1] == "normal" )
          {
            # normal distribution; parameters are mean and sd:
            mean <- as.numeric(distpar[2])
            sd   <- as.numeric(distpar[3])
            samp[,j] <- rnorm(sampsize,mean=mean,sd=sd)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "NormalTrunc" | 
                 distpar[1] == "normaltrunc" )
          {
            # truncated normal distribution; parameters are mean, sd, min and max
            # of untruncated normal distribution
            mean <- as.numeric(distpar[2])
            sd   <- as.numeric(distpar[3])
            min  <- as.numeric(distpar[4])
            max  <- as.numeric(distpar[5])
            cdf.min <- pnorm(min,mean=mean,sd=sd)
            cdf.max <- pnorm(max,mean=mean,sd=sd)
            samp[,j] <- runif(sampsize,min=cdf.min,max=cdf.max)
            samp[,j] <- qnorm(samp[,j],mean=mean,sd=sd)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Lognormal" | 
                 distpar[1] == "lognormal" )
          {
            # lognormal distribution; parameters are mean and sd:
            mean    <- as.numeric(distpar[2])
            sd      <- as.numeric(distpar[3])
            sdlog   <- sqrt(log(1+sd^2/mean^2))
            meanlog <- log(mean) - 0.5*sdlog^2
            samp[,j] <- rlnorm(sampsize,meanlog=meanlog,sdlog=sdlog)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "LognormalTrunc" | 
                 distpar[1] == "lognormaltrunc" )
          {
            # truncated lognormal distribution; parameters are mean, sd, min and max
            # of untruncated lognormal distribution
            mean    <- as.numeric(distpar[2])
            sd      <- as.numeric(distpar[3])
            sdlog   <- sqrt(log(1+sd^2/mean^2))
            meanlog <- log(mean) - 0.5*sdlog^2
            min     <- as.numeric(distpar[4])
            max     <- as.numeric(distpar[5])
            cdf.min <- plnorm(min,meanlog=meanlog,sdlog=sdlog)
            cdf.max <- plnorm(max,meanlog=meanlog,sdlog=sdlog)
            samp[,j] <- runif(sampsize,min=cdf.min,max=cdf.max)
            samp[,j] <- qlnorm(samp[,j],meanlog=meanlog,sdlog=sdlog)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Inv" | 
                 distpar[1] == "inv" )
          {
            # inverse distribution (f(x) prop. 1/x); 
            # parameters are min and max:
            min     <- as.numeric(distpar[2])
            max     <- as.numeric(distpar[3])
            log.min <- log(min)
            log.max <- log(max)
            samp[,j] <- runif(sampsize,min=0,max=1)
            samp[,j] <- exp(samp[,j]*(log.max-log.min)+log.min)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Exponential" | 
                 distpar[1] == "exponential" )
          {
            # exponential distribution; parameter is mean:
            mean <- as.numeric(distpar[2])
            samp[,j] <- runif(sampsize,min=0,max=1)
            samp[,j] <- -mean * log(1-samp[,j])
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Delta" | 
                 distpar[1] == "delta" )
          {
            # delta distribution; parameter is mean:
            mean     <- as.numeric(distpar[2])
            samp[,j] <- rep(mean,sampsize)   ### xxx von Nele geändert rep(sampsize,mean)
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          if ( distpar[1] == "Discrete" | 
                 distpar[1] == "discrete" )
          {
            # discrete distribution; parameters are probabilities:
            n <- length(distpar[-1])/2
            probs <- as.numeric(distpar[1+n+(1:n)])
            probs <- probs/sum(probs)
            probs.lowerbounds <- rep(NA,n)
            probs.lowerbounds[1] <- 0
            if ( n > 1 )
            {
              for ( i in 2:n ) probs.lowerbounds[i] <- 
                probs.lowerbounds[i-1] + probs[i-1]
            }
            samp.unif <- runif(sampsize,min=0,max=1)
            ind <- rep(NA,sampsize)
            for ( i in 1:sampsize )
            {
              ind[i] <- sum(samp.unif[i] > probs.lowerbounds)
            }
            samp[,j] <- distpar[1+ind]
            dist.found <- T
          }
        }
        if ( !dist.found )
        {
          stop(paste("Distribution",distpar[1],"not yet implemented"))
        }
      }
    }
    else
    {
      stop(paste("sysanal.randsamp: unknown distribution type:",dist))
    }
  }
  
  # collect results:
  res <- list(
    mean       = mean,
    sd         = sd,
    cor        = R,
    sampsize   = sampsize,
    sample     = samp,
    logsamppdf = logpdf
  )
  
  # write results:
  if ( !is.na(file) )
  {
    write.table(data.frame(samp,logsamppdf=logpdf),file=file,
                col.names=TRUE,row.names=FALSE,sep="\t")
  }
  
  # return results:
  return(res)
}


################################################################################

# generate a random sample from a uniform distribution in a ball of given radius

# see introduction of Harman, R. and Lacko, V., On decompositional algorithms
# for uniform sampling from n-spheres and n-balls, Journal of Multivariate 
# Analysis 101, 2297-2304, 2010 for details and references 

sysanal.randsamp.ball <- function(sampsize=1,dim=1,radius=1)
{
  res <- matrix(data=rnorm(sampsize*dim),nrow=sampsize)
  norm <- function(x) { return(sqrt(sum(x^2))) }
  rnorms <- apply(res,1,norm)
  res <- diag(radius/rnorms) %*% res
  res <- diag(runif(sampsize)^(1/dim)) %*% res
  return(res) 
}

################################################################################

# Calculate a Metropolis Markov Chain sample of a distribution

sysanal.markovchain.metropolis <- function(log.pdf,z.ini,
                                           jump.sd,jump.cor=0,
                                           sampsize,thin=1,...)
{
  # -----------------------------------------------------------------------
  # This function calculates a Markov Chain of a probability distribution
  # of a vector of continuous random variables using the Metropolis
  # algorithm with a normal jump distribution with given standard
  # deviations and correlation structure.
  # The log of the probability density of the distribution must be 
  # specified as a function log.pdf(z,...) where z is the vector of values
  # for which the density has to be evaluated.
  #
  # Arguments:
  # log.pdf      function "log.pdf(z,...)" that calculates the log of the
  #              probability density function of the vector of random 
  #              variables that are to be sampled.
  #              log.pdf must be given a single vector z at which the
  #              density has to be evaluated. Additional arguments from
  #              the call to calc.markovchain.metropolis will be passed on.
  #              If the probability density is zero, NA must be returned.
  # z.ini        vector of values at which the chain is to be started.
  # jump.sd      vector of standard deviations of the jump distribution.
  # jump.cor     correlation matrix of jump distribution or NA if all
  #              correlations are zero.
  # sampsize     sample size (length of the chain)
  # thin         factor with which to thin storage of results (thin=n: 
  #              only each nth result is returned; this saves memory)
  # ...          additional parameters passed to "log.pdf"
  #
  # Return Value:
  # List with the following elements:
  # z            sample as a matrix with sample points in its rows.
  # log.pdf      vector with log pdf values of the sample.
  # reject.freq  rejection frequency of the jumps.
  # error        error message (empty string if no error occurred).
  #
  #                                first version:         Dec.  08, 2007 PR
  #                                add parameter "thin":  March 26, 2008 PR
  #                                minor modification:    March 28, 2009 PR
  # -----------------------------------------------------------------------
  
  # set up and initialize arrays:
  
  returnsize         <- floor(sampsize/thin)+1
  z.sample           <- matrix(data=NA,nrow=returnsize,ncol=length(z.ini))
  colnames(z.sample) <- names(z.ini)
  log.pdf.sample     <- rep(NA,returnsize)
  reject.freq        <- 1
  error              <- ""
  
  # calculate Cholesky decomposition of variance-covariance matrix:
  
  R <- diag(rep(1,length(jump.sd)))
  if ( is.matrix(jump.cor) ) R <- jump.cor
  if ( (nrow(R)!=length(jump.sd)) | (ncol(R)!=length(jump.sd)) )
  {
    error <- paste("sysanal.markovchain.metropolis:",
                   "illegal dimension of correlation matrix")
  }
  if ( nchar(error) == 0 )
  {
    sigma <- diag(jump.sd) %*% R %*% diag(jump.sd)
    A.chol <- try(t(chol(sigma)),silent=FALSE)
    if( inherits(A.chol,"try-error") )
    {
      error <- paste("sysanal.markovchain.metropolis:",
                     "unable to calculate Cholesky decomposition of variance-covariance matrix")
    }
  }
  
  # initialize Markov chain:
  
  if ( nchar(error) == 0 )
  {
    z.current         <- z.ini
    log.pdf.current   <- log.pdf(z.current,...)
    z.sample[1,]      <- z.current
    log.pdf.sample[1] <- log.pdf.current
    if ( is.na(log.pdf.sample[1]) )
    {
      error <- paste("sysanal.markovchain.metropolis:",
                     "probability density is zero at initial point of chain")
    }
  }
  if ( nchar(error) == 0 )
  {
    num.accept <- 0
    for ( i in 2:returnsize )
    {
      for ( j in 1:thin )
      {
        # calculate suggested new sample point:
        
        jump.unif <- runif(length(z.ini),min=0,max=1)
        jump <- A.chol %*% qnorm(jump.unif)
        z.suggested <- z.current + as.vector(jump)
        
        # calculate log pdf at suggested new sample point:
        
        log.pdf.suggested <- log.pdf(z.suggested,...)
        
        # accept new point with probability r=pdf.suggested/pdf.prev
        
        accept <- FALSE
        if ( is.finite(log.pdf.suggested) )
        {
          if ( log.pdf.suggested > log.pdf.current )
          {
            accept <- TRUE
          }
          else
          {
            r <- exp(log.pdf.suggested-log.pdf.sample[i-1])
            if ( runif(n=1,min=0,max=1) <= r )
            {
              accept <- TRUE
            }
          }
        }
        if ( accept == TRUE )
        {
          z.current       <- z.suggested
          log.pdf.current <- log.pdf.suggested
          num.accept      <- num.accept+1
        }
        reject.freq <- ((i-2)*thin+j-num.accept)/((i-2)*thin+j)
      }
      z.sample[i,]      <- z.current
      log.pdf.sample[i] <- log.pdf.current
    }
  }
  
  # collect and return results:
  
  res <- list(z           = z.sample,
              log.pdf     = log.pdf.sample,
              reject.freq = reject.freq,
              error       = error)
  return(res)
}

################################################################################

# Calculate an importance sample of a distribution

sysanal.impsamp <- function(log.pdf,z,z.log.pdf,...)
{
  w <- rep(NA,nrow(z))
  log.pdf.values <- w
  for ( i in 1:nrow(z) )
  {
    par <- z[i,]
    names(par) <- colnames(z)
    log.pdf.values[i] <- log.pdf(par,...)
  }
  log.pdf.max <- max(log.pdf.values,na.rm=TRUE)
  w <- exp(log.pdf.values-log.pdf.max-z.log.pdf)
  w <- ifelse ( is.na(w),0,w )
  w <- w/sum(w)
  ess <- sum(w)^2/sum(w^2)
  return(list(z=z,w=w,ess=ess,
              log.pdf.dist=log.pdf.values,log.pdf.samp=z.log.pdf))
}


################################################################################

# function to improve jump distribution;

sysanal.jump.dist <- function(postsamp,fact.sd,fract.burnin=0,fact.cor=1,
                              plot=F)
{
  ind.end   <- nrow(postsamp)
  ind.start <- as.integer(fract.burnin*(ind.end-1))+1
  sd.postsamp <- apply(postsamp,2,sd)
  postsamp.local <- postsamp[ind.start:ind.end,sd.postsamp!=0]
  if ( is.vector(postsamp.local) )
  {
    postsamp.local <- as.matrix(postsamp.local,nrow=length(postsamp.local))
    colnames(postsamp.local) <- colnames(postsamp)[sd.postsamp!=0]
    postsamp.local <- data.frame(postsamp.local)
  }
  sd  <- fact.sd*apply(postsamp.local,2,sd)
  cor <- NA
  if ( ncol(postsamp.local) > 1 )
  {
    corr <- cor(postsamp.local)
    if ( plot )
    {
      image(x=1:ncol(corr),y=1:nrow(corr),z=abs(corr),zlim=c(0,1),
            col=grey((100:0)/100),xlab="variable index",ylab="variable index",
            main="structure of correlation matrix")
      abline(v=0.5)
      abline(h=0.5)
      abline(v=ncol(corr)+0.5)
      abline(h=nrow(corr)+0.5)
    }
    corr <- fact.cor*corr
    diag(corr) <- rep(1,nrow(corr))
    try(chol(corr)) # test if Cholesky factorization works 
    # (important for subsequent sampling)
  }
  return(list(sd=sd,cor=corr))
}

################################################################################

# function to calculate probability densities of univariate distributions:
# ------------------------------------------------------------------------

sysanal.calcpdf <- function(x,distpar,log=FALSE)
{
  if ( distpar[1] == "Uniform" | distpar[1] == "uniform" )
  {
    # uniform distribution; parameters are min and max
    min <- as.numeric(distpar[2])
    max <- as.numeric(distpar[3])
    return(dunif(x,min=min,max=max,log=log))
  }
  if ( distpar[1] == "Normal" | distpar[1] == "normal" )
  {
    # normal distribution; parameters are mean and sd:
    mean <- as.numeric(distpar[2])
    sd   <- as.numeric(distpar[3])
    return(dnorm(x,mean=mean,sd=sd,log=log))
  }
  if ( distpar[1] == "NormalTrunc" | distpar[1] == "normaltrunc" )
  {
    # truncated normal distribution; parameters are mean, sd, min and max
    mean <- as.numeric(distpar[2])
    sd   <- as.numeric(distpar[3])
    min  <- as.numeric(distpar[4])
    max  <- as.numeric(distpar[5])
    fact <- 1/(pnorm(q=max,mean=mean,sd=sd)-pnorm(q=min,mean=mean,sd=sd))
    if ( !log )
    {
      return(ifelse(x<min|x>max,0,fact*dnorm(x,mean=mean,sd=sd)))
    }
    else
    {
      return(ifelse(x<min|x>max,-Inf,
                    log(fact)+dnorm(x,mean=mean,sd=sd,log=TRUE)))
    }
  }
  if ( distpar[1] == "Lognormal" | distpar[1] == "lognormal")
  {
    # lognormal distribution; parameters are mean and sd;
    # R parameters are mean and sd of the log of the random variable
    mean    <- as.numeric(distpar[2])
    sd      <- as.numeric(distpar[3])
    sdlog   <- sqrt(log(1+sd^2/mean^2))
    meanlog <- log(mean) - 0.5*sdlog^2
    return(dlnorm(x,meanlog=meanlog,sdlog=sdlog,log=log))
  }
  if ( distpar[1] == "LognormalTrunc" | distpar[1] == "lognormaltrunc" )
  {
    # truncated lognormal distribution; parameters are mean, sd, min and max;
    # R parameters are mean and sd of the log of the random variable
    mean    <- as.numeric(distpar[2])
    sd      <- as.numeric(distpar[3])
    sdlog   <- sqrt(log(1+sd^2/mean^2))
    meanlog <- log(mean) - 0.5*sdlog^2
    min     <- as.numeric(distpar[4])
    max     <- as.numeric(distpar[5])
    fact <- 1/(plnorm(q=max,meanlog=meanlog,sdlog=sdlog)-plnorm(q=min,meanlog=meanlog,sdlog=sdlog))
    if ( !log )
    {
      return(ifelse(x<min|x>max,0,fact*dlnorm(x,meanlog=meanlog,sdlog=sdlog)))
    }
    else
    {
      return(ifelse(x<min|x>max,-Inf,                                                ########### xxxnis
                    log(fact)+dlnorm(x,meanlog=meanlog,sdlog=sdlog,log=TRUE)))
    }
  }
  if ( distpar[1] == "Inv" | distpar[1] == "inv" )
  {
    # inverse distribution; parameters are min and max:
    min <- as.numeric(distpar[2])
    max <- as.numeric(distpar[3])
    if ( !log )
    {
      return(ifelse(x<min|x>max,0,1/(log(max/min)*x)))   
    }
    else
    {
      return(ifelse(x<min|x>max,-Inf,-log(log(max/min)) - log(x)))   
    }
  }
  if ( distpar[1] == "Exponential" | distpar[1] == "exponential" )
  {
    # exponential distribution; parameter is mean:
    mean <- as.numeric(distpar[2])
    if ( !log )
    {
      return(ifelse(x<0,0,1/mean*exp(-x/mean)))   
    }
    else
    {
      return(ifelse(x<0,NA,-log(mean)-x/mean))   
    }
  }
  stop(paste("Distribution",dist,"not yet implemented"))
}

################################################################################

# Calculate the logarithm of the probability density function of a multivariate
# normal or lognormal distribution or of a product of independent marginals

sysanal.calcpdf_mv <- function(z,dist="normal",mean=0,sd=1,cor=0,
                               cor.inv=NA,log=TRUE,distdef=NA,file=NA)
{
  # -----------------------------------------------------------------------
  # This function calculates the logarithm of the probability density 
  # function of a multivariate normal or lognormal distribution or of 
  # independent 1d distributions.
  #
  # Arguments:
  # z:          vector, matrix or data frame at which the logarithm of the 
  #             probability density function has to be evaluated
  # dist:       distribution type: "Normal", "Lognormal" or "Indep".
  # mean:       vector of means
  # sd:         vector of standard deviations
  # cor:        correlation matrix of the distribution
  # cor.inv:    inverse of correlation matrix of the distribution
  #             (alternative input to cor, saves computation time for
  #             repeated calls)
  # log:        if TRUE returns log of pdf, otherwise pdf
  # distdef:    distribution definition for independent 1d distributions
  #
  # Return Value:
  # probability density (or log of) for all sample points in x
  #
  #                                        Peter Reichert    Jan.  01, 2004
  #                                        last modification Oct.  29, 2011
  # -----------------------------------------------------------------------
  
  # consistency checks and initializations:
  mean <- as.vector(mean)
  sd <- as.vector(sd)
  if ( length(sd) != length(mean) )
  {
    stop("sysanal.calcpdf_mv: illegal dimension of standard deviations")
  }
  if ( is.vector(z) )
  {
    len <- length(z)
    names <- names(z)
    z <- as.matrix(z)
    dim(z) <- c(1,len)
    if ( length(names) == len ) colnames(z) <- names
  }
  numpar <- ncol(z)
  R <- diag(rep(1,numpar))
  if ( is.matrix(cor) ) R <- cor
  if ( nrow(R) != numpar || ncol(R) != numpar )
  {
    stop("sysanal.calcpdf_mv: illegal dimension of correlation matrix")
  }
  
  # calculate logarithm of probability density function:
  sampsize <- nrow(z)
  logpdf <- numeric(sampsize)
  if ( dist == "normal" | dist == "Normal" )
  {
    # multivariate normal distribution:
    n <- length(sd)
    sigma <- diag(sd,nrow=n,ncol=n) %*% R %*% diag(sd,nrow=n,ncol=n)
    if ( is.matrix(cor.inv) )
    {
      R.inv <- cor.inv
    }
    else
    {
      R.inv <- solve(R)
    }
    det.R.inv <- det(R.inv)
    if ( det.R.inv > 0 )
    {
      for ( i in 1:sampsize )
      {
        v <- as.matrix(z[i,]-mean,nrow=numpar,ncol=1)/
          as.matrix(sd,nrow=numpar,ncol=1)
        logpdf[i] <- -numpar/2*log(2*pi) + 0.5*log(det.R.inv) - log(prod(sd)) -
          0.5 * t(v) %*% R.inv %*% v
      }
    }
    else
    {
      logpdf <- rep(NA,sampsize)
    }
  }
  else
  {
    if ( dist == "lognormal" | dist == "Lognormal" )
    {
      # multivariate lognormal distribution:
      sdlog    <- sqrt(log(1+sd*sd/(mean*mean)))
      meanlog  <- log(mean) - sdlog*sdlog/2
      if ( numpar > 1 )
      {
        n <- length(sdlog)
        ln.sigma <- log( 1 + diag(sqrt(exp(sdlog*sdlog)-1),nrow=n,ncol=n) %*% 
                           R %*% diag(sqrt(exp(sdlog*sdlog)-1),nrow=n,ncol=n) )
      }
      else
      {
        ln.sigma <- as.matrix( log( 1 + sqrt(exp(sdlog*sdlog)-1)^2 ) )
      }
      ln.sigma.inv = solve(ln.sigma)
      det.ln.sigma = det(ln.sigma)
      if ( det.ln.sigma > 0 )
      {
        for ( i in 1:sampsize )
        {
          if ( min(z[i,]) <= 0 )
          {
            logpdf[i] <- NA 
          }
          else
          {
            v <- log(z[i,])-meanlog
            v <- matrix(v,nrow=numpar,ncol=1) # ensure that v is a column vector
            logpdf[i] <- -numpar/2*log(2*pi) - 0.5*log(det.ln.sigma) - 
              log(prod(z[i,])) - 0.5 * t(v) %*% ln.sigma.inv %*% v
          }
        }
      }
      else
      {
        logpdf <- rep(NA,sampsize)
      }
    }
    else
    {
      if ( dist == "indep" | dist == "Indep" ) 
      {
        for ( i in 1:sampsize )
        {
          logpdf[i] <- 0
          for ( j in 1:ncol(z) )
          {
            ind <- j
            if ( length(colnames(z)) == ncol(z) )
            {
              ind <- match(colnames(z)[j],names(distdef))
              if ( is.na(ind) )
              {
                stop(paste("error in calcpdf_mv:",
                           "variable",colnames(z)[j],"not found"))
              }
            }
            logpdf[i] <- logpdf[i] + sysanal.calcpdf(z[i,j],
                                                     distdef[[ind]],
                                                     log=TRUE)
          }
        }
      }
      else
      {
        stop(paste("sysanal.calcpdf_mv: unknown distribution type:",
                   dist))
      }
    }
  }
  
  # write results:
  if ( !is.na(file) )
  {
    write.table(data.frame(z,logpdf=logpdf),file=file,
                col.names=TRUE,row.names=FALSE,sep="\t")
  }
  
  # return result:
  if ( log )
  {
    return(logpdf)
  }
  else
  {
    return(exp(logpdf))
  }
}

################################################################################

# Function implementing a simple standard likelihood function for a model
# provided by a deterministic function. This function provides the density of
# a multivariate normal distribution centered at the results of the 
# deterministic model.
# If there exists a component of the parameter vector labelled "sd_rel" 
# instead of the provided standard deviations, "error.sd" the standard 
# deviations "error.sd * par["sd_rel"]" are used. This makes it possible to
# estimate a common factor of a given variance-covariance structure.
# This likelihood implementation serves as a template for implementing more
# specific likelihood functions. It is called by "sysanal.logposterior".

sysanal.loglikeli <- function(par,model,L,y,
                              error.sd=NA,error.cor.inv=NA,...)
{
  par.local <- par; if (!is.matrix(par)) par.local = t(as.matrix(par))
  sd.rel <- rep(1,nrow(par.local))
  ind.sd_rel <- match("sd_rel",colnames(par.local))
  if ( !is.na(ind.sd_rel) ) sd.rel <- par.local[,ind.sd_rel]
  res <- rep(NA,nrow(par.local))
  for ( i in 1:nrow(par.local) )
  {
    mean <- model(par.local[i,],L=L,...)
    sd <- error.sd; if ( is.na(sd[1]) ) sd <- rep(1,length(mean))
    if ( is.na(error.cor.inv[1]) )
    {
      res[i] <- 0;
      for ( i in 1:length(mean) )
      {
        res[i] <- res[i] + dnorm(y[i],mean[i],sd[i]*sd.rel[i],log=TRUE)
      }
    }
    else 
    {
      res[i] <- sysanal.calcpdf_mv(z=y,dist="normal",
                                   mean=mean,sd=sd*sd.rel[i],
                                   cor.inv=error.cor.inv)
    }
  }
  return(res)
}

################################################################################

# Function implementing a simple standard posterior probability density (up to
# a multiplicative normalizing constant). It is based on a likelihood 
# definition provided by the user and a multivariate normal or
# lognormal distribution at the prior.
# This log posterior implementation serves as a template for implementing more
# specific posteriors.

sysanal.logposterior <- function(par,model,L,y,
                                 prior.dist="lognormal",prior.mean=1,prior.sd=1,
                                 prior.cor=NA,prior.def=NA,
                                 loglikeli=sysanal.loglikeli,
                                 ...)
{
  logprior  <- sysanal.calcpdf_mv(z=par,dist=prior.dist,mean=prior.mean,
                                  sd=prior.sd,cor=prior.cor,distdef=prior.def)
  loglik <- rep(NA,length(logprior))
  if ( length(logprior) == 1 )
  {
    if ( !is.na(logprior) )
    {
      loglik <- loglikeli(par=par,model=model,L=L,y=y,...)
    }
  }
  else
  {
    ind <- !is.na(logprior)
    if ( sum(ind) > 0 )
    {
      loglik[ind] <- loglikeli(par=par[ind,],model=model,L=L,y=y,...)
    }
  }
  return(logprior+loglik)
}


################################################################################

# sysanal.sens.loc
# ================

# purpose:
# calculate the local sensitivity matrix of a model
# (matrix of partial derivatives of model results with respect to parameters)

# arguments:
# par:         named parameter vector; passed to the model as separate arguments
# model        function representing the model; typically returns a vector,
#              but could also return a matrix
#              (cf. sysanal.model)
# par.inc:     increments of parameters used to approximate the derivatives
# ...          further arguments are passed to model

# output:
# matrix of partial derivatives of model results (rows) 
# with respect to parameters (columns)
# or 3 dimensional array of partial derivatives if model output is a matrix

sysanal.sens.loc <- function(par,model,par.inc=0.01*par,...)
{
  if ( length(par) != length(par.inc) ) 
    stop("*** error in sysanal.sens.loc: unequal length of par and par.inc")
  if ( min(abs(par.inc)) == 0 )
    stop("*** error in sysanal.sens.loc: elements of par.inc must not be zero")
  res.par <- model(par,...)
  V = NA
  if ( is.vector(res.par) )
  {
    V <- matrix(NA,nrow=length(res.par),ncol=length(par))
    colnames(V) <- names(par)
    rownames(V) <- names(res.par)
    for( j in 1:length(par) )
    {
      par.j <- par
      par.j[j] <- par.j[j] + par.inc[j]
      V[,j] <- ( model(par.j,...) - res.par ) / par.inc[j]
    }
  }
  else
  {
    if ( is.matrix(res.par) )
    {
      V <- array(NA,dim=c(dim(res.par),length(par)),
                 dimnames=list(rownames(res.par),
                               colnames(res.par),
                               names(res.par)))
      for ( j in 1:length(par) )
      {
        par.j <- par
        par.j[j] <- par.j[j] + par.inc[j]
        V[,,j] <- ( model(par.j,...) - res.par ) / par.inc[j]
      }
    }
  }
  return(V)
}

# version of the same function that calls a model with explicitly listed 
# arguments (this is primarily for internal use when internally calling nls)
# model.name is the the name of the function with explicit parameter arguments

sysanal.sens.loc.explicitpars <- function(par,model.name,par.inc=0.01*par,x)
{
  if ( length(par) != length(par.inc) ) 
    stop("*** error in sysanal.sensfun: unequal length of par and par.inc")
  if ( min(abs(par.inc)) == 0 )
    stop("*** error in sysanal.sensfun: elements of par.inc must not be zero")
  args <- as.list(par)
  if ( length(x) > 0 ) args <- c(args,x)
  res.par <- do.call(model.name,args=args)
  V <- matrix(NA,nrow=length(res.par),ncol=length(par))
  colnames(V) <- names(par)
  rownames(V) <- names(res.par)
  for( j in 1:length(par) )
  {
    par.j <- par
    par.j[j] <- par.j[j] + par.inc[j]
    args <- as.list(par.j)
    if ( length(x) > 0 ) args <- c(args,x)
    res.j <- do.call(model.name,args=args)
    V[,j] <- ( res.j - res.par ) / par.inc[j]
  }
  return(V)
}

################################################################################

# sysanal.sens.var.samp
# ---------------------

# function to calculate the first order variance based sensitivity
# coefficients from a parameter sample and the corresponding result sample

# Call:      sysanal.sens.var.samp(parsamp,ressamp,nbin=NA,
#                                  method="smooth",order=2,
#                                  bandwidth=1,span=0.75,sd.rel=0.05,
#                                  plot=F)					
# ------------------------------------------------------------------										
#													
# parsamp	   matrix containing the parameter sample (each row corresponds to
#            a sample point)						
#													
# ressamp	   vector (of only one result per parameter sample point) or matrix 
#            of results corresponding to the parameter sample (each row 
#            provides the results corresponding to the parameter values in 
#            the same row of parsamp)									
#													
# nbin 	     number of quantile intervals for which the conditional means	
#		         are calculated, default is the square root of the sample size								
#
# method     "smooth", "loess", "glkerns", "lokerns", "lpepa", "lpridge"	
#            routine to be used for smoothing
#
# order		   order of local regression polynomial or of kernel
#
# bandwidth  method="lpepa" or method="lpdidge" only: bandwidth of the
#            smoothing algorithm
#
# span       method="loess" only: fration of points used for local regression 
#
# sd.rel     method="smooth" only: standard deviation of normal distribution of 
#            smoothing algorithm relative to the 99% quantile interval
#
#	plot       logical variable indicating if a scatter plot of the relationship
#            between parameter and model output should be plotted								
#													
# Output:												
# -------												
#													
# List with two elements:	
#								
# var	       vector with the total variance of each column (each model	 
#            output) of the ressamp matrix						
#													
# var.cond.E matrix with the variance of the conditional expected value	
#			       of each parameter (columns) and each model output (rows)	

sysanal.sens.var.samp <- function(parsamp,ressamp,nbin=NA,
                                  method="smooth",order=2,
                                  bandwidth=1,span=0.75,sd.rel=0.1,
                                  plot=F)
{
  sysanal.package("lpridge")
  sysanal.package("lokern")
  
  if ( is.vector(ressamp) ) ressamp <- as.matrix(ressamp,ncol=1)
  npar  <- ncol(parsamp)
  nsamp <- nrow(parsamp)
  nres  <- ncol(ressamp)
  if ( is.na(nbin) ) nbin <- ceiling(sqrt(nsamp))
  if ( nrow(parsamp ) != nrow(ressamp)) 
  {
    stop ("ressamp and parsamp do not have the same row length")
  }
  
  var_k <- rep(NA,nres)
  names(var_k) <- colnames(ressamp)
  for ( k in 1:nres ) var_k[k] <- var(ressamp[,k])
  
  var_q <- data.frame(matrix(NA,nrow=nres,ncol=npar))
  colnames(var_q) <- colnames(parsamp)
  rownames(var_q) <- colnames(ressamp)
  
  for ( i in 1:npar )
  {
    q  <- quantile(parsamp[,i],probs=(((1:nbin)-0.5)/nbin),
                   na.rm=FALSE,names=TRUE,type=7)
    for ( k in 1:nres )
    {
      if ( method == "lpepa" )
      {
        mean_res <- lpepa(parsamp[,i],ressamp[,k],bandwidth=bandwidth,
                          x.out=q,order=order)$est
      }
      else
      {
        if ( method == "lpridge" )
        {
          mean_res <- lpridge(parsamp[,i],ressamp[,k],bandwidth=bandwidth,
                              x.out=q,order=order)$est
        }
        else
        {
          if ( method == "glkerns" )
          {
            mean_res <- glkerns(parsamp[,i],ressamp[,k],
                                x.out=q,korder=order,hetero=T)$est
          }
          else
          {
            if ( method == "lokerns" )
            {
              mean_res <- lokerns(parsamp[,i],ressamp[,k],
                                  x.out=q,korder=order,hetero=T)$est
            }
            else
            {
              if ( method == "loess" )
              {
                #mean_res <- lowess(parsamp[,i],ressamp[,k],f,
                #                   x.out=q,korder=order,hetero=T)$est
                res <- 
                  loess(y~x,
                        data=data.frame(x=parsamp[,i],y=ressamp[,k]),
                        span=span,degree=order)
                mean_res <- predict(res,newdata=data.frame(x=q))
              }
              else
              {
                if ( method == "smooth" )
                {
                  if ( order == 2 ) m <- "quadratic"
                  else              m <- "linear"
                  mean_res <- sysanal.smooth(parsamp[,i],ressamp[,k],
                                             sigma=((q[nbin]-q[1])*sd.rel),
                                             newx=q,
                                             method=m)$y
                }
                else
                {
                  stop(paste("calc.var.sens: method",
                             method,
                             "not implemented"))
                }
              }
            }
          }
        }
      }
      
      var_q[k,i] <- var(mean_res)
      if ( plot )
      {
        plot(parsamp[,i],ressamp[,k],pch=19,cex=0.2,
             xlab=colnames(parsamp)[i],ylab="Y")
        lines(q,mean_res,lwd=3,col="red")
      }
    }
  }
  
  return (list(var=var_k,var.cond.E=var_q))
}

################################################################################

# calculation of index combinations (auxiliary function used in ident)
# ====================================================================

sysanal.comb <- function(n,p)
{
  # -----------------------------------------------------------------------
  # This function calculates all combination of subsets of length p out
  # of n indices.
  #
  # Arguments:
  # n:   number of indices.
  # p:   length of subset of indices.
  #
  # Return Value:
  # matrix with subsets of length p as rows.
  #
  #                                         Peter Reichert    Dec. 27, 2002
  # -----------------------------------------------------------------------
  
  # check input:
  if ( p > n ) stop("comb: illeal arguments (p>n)")
  
  # initialize array and auxiliary variables:
  num.comb <- choose(n,p)
  comb <- matrix(nrow=num.comb,ncol=p)
  ind <- 1:p
  pointer <- p
  
  # calculate index combinations:
  for ( i in 1:num.comb )
  {
    comb[i,] <- ind
    ind[pointer] <- ind[pointer]+1
    if ( ind[pointer] > n )
    {
      while ( pointer > 1 )
      {
        pointer <- pointer-1
        if ( ind[pointer] < n-(p-pointer) )
        {
          ind[pointer] <- ind[pointer]+1
          for ( j in (pointer+1):p )
          {
            ind[j] <- ind[pointer]+j-pointer
          }
          pointer <- p
          break
        }
      }
    }
  }
  
  # return results:
  return(comb)
}

################################################################################

# calculation of collinearity index (auxiliary function used in ident)
# ====================================================================

sysanal.collind <- function(sen.scaled)
{
  # -----------------------------------------------------------------------
  # This function calculates the collinearity index from a scaled 
  # sensitivity matrix.
  #
  # Arguments:
  # sen:       matrix of model sensitivities (scaled partial derivatives
  #            of model outcomes with respect to model parameters:
  #            delta.par/scale dy/dpar); the columns of sen refer to
  #            different model parameters, the rows to different model
  #            outcomes.
  #
  # Return Value:
  # collinearity index (real number).
  #
  #                                         Peter Reichert    Dec. 27, 2002
  # -----------------------------------------------------------------------
  
  # normalize sensitivity functions:
  num.par <- ncol(sen.scaled)
  norms <- numeric(num.par)
  for ( i in 1:num.par )
  {
    norms[i] <- sqrt( sen.scaled[,i] %*% sen.scaled[,i] )
  }
  sen.norm <- sen.scaled %*% diag( 1/norms )
  
  # calculate collinearity index:
  collind <- 1/sqrt(min(eigen( t(sen.norm) %*% sen.norm )$values))
  
  # return result:
  return(collind)
}

################################################################################

# calculation of identifiability measures
# =======================================

sysanal.ident <- function(sen,delta.par=0,scale=0,max.subset.size=0)
{
  # -----------------------------------------------------------------------
  # This function calculates a parameter sensitivity ranking and 
  # collinearity indices for a series of parameter combinations
  # based on liniear sensitivity functions of a model, parameter
  # uncertainty ranges and scale factors of model results.
  # This function is a simplified version of the program "ident"
  # available at http://www.ident.eawag.ch
  #
  # Arguments:
  # sen:       matrix of model sensitivities (partial derivatives
  #            of model outcomes with respect to model parameters:
  #            dy/dpar); the columns of sen refer to different 
  #            model parameters, the rows to different model outcomes.
  # delta.par: model parameter uncertainty ranges (length equal to 
  #            the number of columns of sen); if zero, then all ranges
  #            are assumed to be unity.
  # scale:     scaling factors of model results (if zero, then all 
  #            scaling factors are assumed to be unity).
  #
  # Return Value:
  # List of delta.msqr, collind.
  #
  #                                         Peter Reichert    Dec. 27, 2002
  # -----------------------------------------------------------------------
  
  # determine number of model parameters:
  num.out <- nrow(sen)
  num.par <- ncol(sen)
  names.par <- colnames(sen)
  if ( length(names.par) != num.par ) names.par <- paste("par",1:num.par,sep="")
  if ( max.subset.size == 0 ) max.subset.size <- min(num.par,4)
  
  # apply parameter uncertainty ranges and scale factors if available:
  sen.scaled <- sen
  if ( length(delta.par) == num.par ) sen.scaled <- sen.scaled %*% diag(delta.par)
  if ( length(scale)     == num.out ) sen.scaled <- diag(1/scale) %*% sen.scaled
  
  # calculate sensitivity ranking:
  delta.msqr <- numeric(num.par)
  names(delta.msqr) <- names.par
  for ( i in 1:num.par )
  {
    delta.msqr[i] <- sqrt( t(sen.scaled[,i]) %*% sen.scaled[,i] ) / sqrt(num.out)
  }
  res <- list(delta.msqr=delta.msqr)
  
  if ( max.subset.size > 1 )
  {
    for ( i in 2:min(max.subset.size,num.par) )
    {
      ind <- sysanal.comb(num.par,i)
      collind <- numeric(nrow(ind))
      par.set <- matrix(nrow=nrow(ind),ncol=i)
      colnames(par.set) <- paste("par.",1:i,sep="")
      for ( j in 1:nrow(ind) )
      {
        collind[j] <- sysanal.collind(sen.scaled[,ind[j,]])
        for ( k in 1:i )
        {
          par.set[j,k] <- names.par[ind[j,k]]
        }
      }
      if ( nrow(par.set) > 1 )
      {
        ind.sorted <- order(collind)
        res[[paste("collind.",i,sep="")]] <- 
          data.frame(par.set[ind.sorted,],collind=collind[ind.sorted])
      }
      else
      {
        res[[paste("collind.",i,sep="")]] <- 
          data.frame(par.set,collind=collind)
      }
    }
  }
  
  # return results:
  return(res)
}

################################################################################

# sysanal.boxcox
# ==============

# purpose:
# Box-Cox transformation

# arguments:
# data:        data to be transformed
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# transformed data

sysanal.boxcox <- function(data,lambda1=1,lambda2=1)
{
  if ( lambda1 == 0 )
  {
    return(ifelse(data>-lambda2,log(data+lambda2),NA))
  }
  else
  {
    return(ifelse(data>=-lambda2,((data+lambda2)^lambda1 - 1)/lambda1,NA))
  }
}

################################################################################

# sysanal.boxcox.deriv
# ====================

# purpose:
# calculate derivative of Box-Cox transformation

# arguments:
# data:        data to be transformed
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# derivative of Box-Cox transformation

sysanal.boxcox.deriv <- function(data,lambda1=1,lambda2=1)
{
  return(ifelse(data>-lambda2,(data+lambda2)^(lambda1 - 1),NA))
}

################################################################################

# sysanal.boxcox.inv
# ==================

# purpose:
# inverse Box-Cox transformation

# arguments:
# data:        data to be back-transformed
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# back-transformed data

sysanal.boxcox.inv <- function(data,lambda1=1,lambda2=1)
{
  if ( lambda1 == 0 )
  {
    return(exp(data)-lambda2)
  }
  else
  {
    return(ifelse(lambda1*data>-1,(lambda1*data+1)^(1/lambda1)-lambda2,
                  -lambda2))
  }
}

################################################################################

# sysanal.logsinh
# ===============

# purpose:
# log-sinh transformation

# arguments:
# data:        data to be transformed
# lambda1:     value of parameter lambda1 (trans. to const. var.) of log-sinh transformation
# lambda2:     value of parameter lambda2 (offset) of log-sinh transformation

# output:
# transformed data

sysanal.logsinh <- function(data,lambda1=1,lambda2=0)
{
  return(lambda1*log(sinh((lambda2+data)/lambda1)))
}

################################################################################

# sysanal.logsinh.deriv
# =====================

# purpose:
# calculate derivative of log-sinh transformation

# arguments:
# data:        data to be transformed
# lambda1:     value of parameter lambda1 (trans. to const. var.) of log-sinh transformation
# lambda2:     value of parameter lambda2 (offset) of log-sinh transformation

# output:
# derivative of log-sinh transformation

sysanal.logsinh.deriv <- function(data,lambda1=1,lambda2=0)
{
  return(1/tanh((lambda2+data)/lambda1))
}

################################################################################

# sysanal.logsinh.inv
# ===================

# purpose:
# inverse log-sinh transformation

# arguments:
# data:        data to be back-transformed
# lambda1:     value of parameter lambda1 (trans. to const. var.) of log-sinh transformation
# lambda2:     value of parameter lambda2 (offset) of log-sinh transformation

# output:
# back-transformed data

sysanal.logsinh.inv <- function(data,lambda1=1,lambda2=0)
{
  return(lambda1*log(exp(data/lambda1)+sqrt(exp(2*data/lambda1)+1)) - lambda2)
}

################################################################################

# sysanal.trans.to.interval
# =========================

# purpose:
# transforms the real axis to an interval (default: unit interval)

# arguments:
# x:       data to be transformed
# min:     minimum of the interval (default: 0)
# max:     maximum of the interval (default: 1)

# output:
# transformed data

sysanal.trans.to.interval <- function(x,min=0,max=1)
{
  y <- 0.5*(min+max) + (max-min)/pi*atan(x)
  return(y)
}

################################################################################

# sysanal.trans.from.interval
# ===========================

# purpose:
# transforms an interval (default: unit interval) to the real axis

# arguments:
# y:       data to be transformed
# min:     minimum of the interval (default: 0)
# max:     maximum of the interval (default: 1)

# output:
# transformed data

sysanal.trans.from.interval <- function(y,min=0,max=1)
{
  x <- tan(0.5*pi*(2*y-max-min)/(max-min))
  return(x)
}

################################################################################

sysanal.trans.par.normal.tovec <- function(mean,sd,cor,trans=T,max.cor=0.5)
{
  n <- length(mean)
  
  par <- rep(NA,n*(n+3)/2)
  par[1:n] <- mean
  par[(n+1):(2*n)] <- sd
  if ( n > 1 )
  {
    k <- 2*n
    for ( i in 1:(n-1) ) 
    {
      par[k+1:(n-i)] <- cor[i,(i+1):n]
      k <- k + n - i
    }
  }
  names(par) <- c(names(mean),names(mean),rep("cor",n*(n-1)/2))
  
  if ( trans )
  {
    par[(n+1):(2*n)] <- log(par[(n+1):(2*n)])
    par[(2*n+1):(n*(n+3)/2)] <- 
      sysanal.trans.from.interval(par[(2*n+1):(n*(n+3)/2)],
                                  min=-max.cor,max=max.cor)
  }
  
  return(par)
}

sysanal.trans.par.normal.fromvec <- function(par,trans=T,max.cor=0.5)
{
  n <- (-3+sqrt(9+8*length(par)))/2
  
  if ( length(par) != n*(n+3)/2 )
  {
    cat("sysanal.trans.par.normal.fromvec:",
        "illegal length of parameter vector:",length(par),"\n")
    mean <- NA
    sd   <- NA
    cor  <- NA
  }
  else
  {
    if ( trans )
    {
      par[(n+1):(2*n)] <- exp(par[(n+1):(2*n)])
      par[(2*n+1):(n*(n+3)/2)] <- 
        sysanal.trans.to.interval(par[(2*n+1):(n*(n+3)/2)],
                                  min=-max.cor,max=max.cor)
    }
    
    mean <- par[1:n]
    sd   <- par[(n+1):(2*n)]
    cor  <- diag(rep(1,n),nrow=n)
    if ( n > 1 )
    {
      k <- 2*n
      for ( i in 1:(n-1) )
      {
        cor[i,(i+1):n] <- par[k+1:(n-i)]
        cor[(i+1):n,i] <- par[k+1:(n-i)]
        k <- k + n - i
      }
    }
    names(mean)   <- names(par)[1:n]
    names(sd)     <- names(par)[1:n]
    rownames(cor) <- names(par)[1:n]
    colnames(cor) <- names(par)[1:n]
  }
  
  return(list(mean=mean,sd=sd,cor=cor))
}            

################################################################################

# sysanal.gnls
# ============

# purpose:
# calculate the generalized least squares parameter estimates of a nonlinear model

# arguments:
# model.name:  name of the function representing the model
#              (note that this model requires the parameters to be specified
#              explicitly in the function headers)
# y.obs:       observed data corresponding to model output
# var.inv:     inverse variance-covariance matrix of the multivariate normal distribution
#              that characterizes the error model
# par:         named parameter vector with initial values; 
#              passed to the model as separate arguments
# x:           list of named inputs passed to the model
# ...          further optional arguments passed to nls

# output:
# list of:
# model.name:  name of the function representing the model
#              (cf. sysanal.model)
# par:         named parameter vector with parameter estimates 
# y.obs:       observed data corresponding to model output
# y.det:       deterministic model results corresponding to the estimated parameters
# resid:       residuals
# var.inv:     inverse variance-covariance matrix of the multivariate normal distribution
#              that characterizes the error model
# x:           input object passed as first argument to the model
# res.nls:     results from application of the nls function

sysanal.gnls <- function(model.name,y.obs,var.inv,par,x=list(),...)
{
  A <- chol(var.inv)
  add.args <- ""
  if ( names(x) > 0 )
  {
    add.args <- paste(",",paste(names(x),collapse=","),sep="")
  }
  model.call <- paste(model.name,"(",paste(names(par),collapse=","),
                      add.args,")",sep="")
  data <- x
  data$y.obs.trans <- A%*%y.obs
  res.nls <- nls(as.formula(paste("y.obs.trans ~ A%*%",model.call)),
                 data=data,
                 start=par,...)
  
  coef  <- coef(res.nls)
  args <- as.list(coef)
  if ( length(x) > 0 ) args <- c(args,x)
  y.det <- do.call(model.name,args=args)  # reevaluate to avoid need for
  resid <- y.obs - y.det                  # back transformation
  
  return(list(model.name = model.name,
              par        = coef,
              y.obs      = y.obs,
              y.det      = y.det,
              resid      = resid,
              var.inv    = var.inv,
              x          = x,
              res.nls    = res.nls))
}

################################################################################

# sysanal.gnls.diag
# =================

# purpose:
# calculate regression diagnostics for generalized least squares

# arguments:
# res.gnls:    output from sysanal.gnls
# par.inc:     increments of parameters used to approximate the derivatives

# output:
# list of:
# par.inc:     parameter increments used to calculate V
# V:           matrix of partial derivatives of model results with respect to
#              parameters at the estimated parameter values
# sd.par:      standard errors of the parameters assuming known error variances
# var.par:     estimated variance-covariance matrix of the parameters 
#              assuming known error variances
# corr.par:    estimated correlation matrix of the parameters
# sd.rel:      estimated correction factor in standard deviations of the error model
# sd.par.rel:  standard errors of the parameters when estimating a common factor
#              in error variances
# var.par.rel: estimated variance-covariance matrix of the parameters 
#              when estimating a common factor in error variances

sysanal.gnls.diag <- function(res.gnls,par.inc)
{
  # sensitivites:
  # -------------
  
  V <- sysanal.sens.loc.explicitpars(par        = res.gnls$par,
                                     model.name = res.gnls$model.name,
                                     par.inc    = par.inc,
                                     x          = res.gnls$x)
  
  # calculate variance-covariance matrix of the estimator for given error variance:
  # -------------------------------------------------------------------------------
  
  var.par  <- solve( t(V) %*% res.gnls$var %*% V )
  rownames(var.par) = names(res.gnls$par)
  colnames(var.par) = names(res.gnls$par)
  
  sd.par   <- sqrt(diag(var.par))
  names(sd.par) <- names(res.gnls$par)
  
  corr.par <- (1/sd.par) %*% t(1/sd.par) * var.par
  rownames(corr.par) = names(res.gnls$par)
  colnames(corr.par) = names(res.gnls$par)
  
  # variance-covariance matrix of the estimator for estimated error variance:
  # -------------------------------------------------------------------------
  
  sd.rel <- sqrt(as.numeric(
    t(res.gnls$y.obs-res.gnls$y.det) %*% res.gnls$var.inv %*% (res.gnls$y.obs-res.gnls$y.det) /
      (length(res.gnls$y.obs)-length(res.gnls$par))
  ))
  var.par.rel <- var.par*sd.rel^2
  sd.par.rel  <- sd.par*sd.rel
  
  return(c(res.gnls,
           list(par.inc     = par.inc,
                V           = V,
                sd.par      = sd.par,
                var.par     = var.par,
                corr.par    = corr.par,
                sd.rel      = sd.rel,
                sd.par.rel  = sd.par.rel,
                var.par.rel = var.par.rel
           )))
}

################################################################################

# sysanal.gnls.test
# =================

# purpose:
# calculate test statistics for generalized least squares

# arguments:
# res.gnls:    output from sysanal.gnls

# output:
# list of:
# V:           matrix of partial derivatives of model results with respect to
#              parameters at the estimated parameter values
# sd.par:      standard errors of the parameters assuming known error variances
# var.par:     estimated variance-covariance matrix of the parameters 
#              assuming known error variances
# corr.par:    estimated correlation matrix of the parameters
# sd.rel:      estimated correction factor in standard deviations of the error model
# sd.par.rel:  standard errors of the parameters
#              when estimating a common factor in error variances
# var.par.rel: estimated variance-covariance matrix of the parameters 
#              when estimating a common factor in error variances

sysanal.gnls.test <- function(res.gnls,V,par,y.det)
{
  # chi2 test; to be compared with 
  # qchisq(1-alpha,length(y.obs)):
  # ------------------------------
  
  chi2 <- t(res.gnls$y.obs-y.det) %*% res.gnls$var.inv %*% (res.gnls$y.obs-y.det)
  
  # exact F test; to be compared with 
  # qf(1-alpha,length(par),length(y.obs)-length(par)):
  # --------------------------------------------------
  
  F.exact <- t(res.gnls$y.obs-y.det) %*% res.gnls$var.inv %*% V %*% 
    solve(t(V)%*%res.gnls$var.inv%*%V) %*% 
    t(V) %*% res.gnls$var.inv %*% (res.gnls$y.obs-y.det) /
    t(res.gnls$y.obs-y.det) %*% 
    ( res.gnls$var.inv -
        res.gnls$var.inv %*% V %*% 
        solve(t(V)%*%res.gnls$var.inv%*%V) %*% 
        t(V) %*% res.gnls$var.inv ) %*% 
    (res.gnls$y.obs-y.det) *
    (length(res.gnls$y.obs)-length(res.gnls$par))/length(res.gnls$par)
  
  # F test for linearized model; to be compared with 
  # qf(1-alpha,length(par),length(y.obs)-length(par)):
  # --------------------------------------------------
  
  F.lin <- t(par-res.gnls$par) %*% t(V) %*% res.gnls$var.inv %*% V  %*% (par-res.gnls$par) /
    ( t(res.gnls$y.obs-res.gnls$y.det) %*% res.gnls$var.inv %*% (res.gnls$y.obs-res.gnls$y.det) ) *
    (length(res.gnls$y.obs)-length(res.gnls$par))/length(res.gnls$par)
  
  return(list(chi2        = chi2,
              F.exact     = F.exact,
              F.lin       = F.lin))
}

################################################################################

# sysanal.gnls.predict
# ====================

# purpose:
# calculate predictions based on generalized least squares regression results

# arguments:
# diag.gnls:   output from sysanal.gnls.diag
# newx:        new model input x at which confidence intervals are to be calculated
# level:       probability level of conficence intervals (default 0.95)
# var:         error variance-covariance matrix at new model input (only variances used)

# output:
# list of:
# x:           x values at which model results are calculated
# y.pred:      predicted model results
# confint:     confidence intervals of deterministic model results
# confint.rel: confidence intervals of deterministic model results
#              when estimating a common factor in error variances
# predint:     prediction intervals (only if error variance.covariance matrix is provided)
# predint.rel: prediction intervals (only if error variance.covariance matrix is provided)
#              when estimating a common factor in error variances

sysanal.gnls.predict <- function(diag.gnls,newx,level=0.95,var=NA)
{
  # prediction and sensitivites:
  # ----------------------------
  
  args <- as.list(diag.gnls$par)
  if ( length(newx) > 0 ) args <- c(args,newx)
  res.par <- do.call(diag.gnls$model.name,args=args)
  newV <- sysanal.sens.loc.explicitpars(par        = diag.gnls$par,
                                        model.name = diag.gnls$model.name,
                                        par.inc    = diag.gnls$par.inc,
                                        x          = newx)
  
  # variance-covariance matrix and standard deviaitons of model results:
  # --------------------------------------------------------------------
  
  var.y    <- newV %*% diag.gnls$var.par %*% t(newV)
  sd.y     <- sqrt(diag(var.y))
  sd.y.rel <- diag.gnls$sd.rel*sd.y
  
  # calculate confidence intervals:
  # -------------------------------
  
  df          <- length(diag.gnls$y.obs)-length(diag.gnls$par)
  confint     <- sysanal.confint(est=res.par,sd=sd.y,df=NA,level=level)
  confint.rel <- sysanal.confint(est=res.par,sd=sd.y.rel,df=df,level=level)
  
  res <- list(x           = newx,
              y.pred      = res.par,
              confint     = confint,
              confint.rel = confint.rel)
  
  if ( is.matrix(var) )
  {
    sd.y.err     <- sqrt(sd.y^2+diag(var))
    sd.y.rel.err <- sqrt(sd.y.rel^2+diag(var)*diag.gnls$sd.rel^2)
    predint     <- sysanal.confint(est=res.par,sd=sd.y.err,df=NA,level=level)
    predint.rel <- sysanal.confint(est=res.par,sd=sd.y.rel.err,df=df,level=level)
    res <- c(res,list(predint=predint,predint.rel=predint.rel))
  }
  
  return(res)
}

################################################################################

# sysanal.confint
# ===============

# purpose:
# calculate confidence intervals based on the Student t distribution
# or on the normal distribution (if df=NA)

# arguments:
# est:         vector of point estimates
# sd:          vector of estimated standard deviations
# df:          degrees of freedom
# level:       probability level of confidence intervals (default 0.95)

# output:
# table of lower and upper bounds of confidence intervals

sysanal.confint <- function(est,sd,df,level)
{
  alpha <- 1- level
  confint <- matrix(NA,nrow=length(est),ncol=2)
  colnames(confint) <- c(paste(100*alpha/2,"%",sep=""),paste(100*(1-alpha/2),"%",sep=""))
  rownames(confint) <- names(est)
  if ( is.na(df) )
  {
    fact <- qnorm(1-alpha/2)
  }
  else
  {
    fact <- qt(1-alpha/2,df)
  }
  confint[,1] <- est - fact*sd
  confint[,2] <- est + fact*sd
  return(confint)
}

################################################################################

# sysanal.resid.diag
# ==================

# purpose:
# plot residual diagnostics plots

# arguments:
# obs:         vector of observations
# calc:        vector of calculated results

# output:
# diagnostic plots

sysanal.resid.diag <- function(obs,calc,header="")
{
  # calculate range of measurements and normalized residuals:
  # ---------------------------------------------------------
  
  obs.min        <- min(obs)
  obs.max        <- max(obs)
  calc.min       <- min(calc)
  calc.max       <- max(calc)
  resid.norm     <- (obs-calc)/sd(obs-calc)
  resid.norm.abs <- abs(resid.norm)
  resid.max      <- max(abs(resid.norm.abs))
  resid.lim      <- 1.1*resid.max*c(-1,1)
  marker         <- 19
  
  # divide plot area into four panels:
  # ----------------------------------
  
  par.def <- par(no.readonly=TRUE)
  par(mfrow=c(2,2),xaxs="i",yaxs="i",mar=c(4.5,4,3,1),oma=c(0,0,2,0))
  
  # plot sequence of residuals:
  # ---------------------------
  
  plot(resid.norm,main="Sequence of Residuals",
       ylim=resid.lim,pch=marker,cex=0.8,
       xlab="Data Points",ylab="Normalized Residuals")
  lines(c(0,length(resid.norm)),c(0,0))
  
  # plot residuals as function of predicted values:
  # -----------------------------------------------
  
  plot(calc,resid.norm,main="Residuals vs. Predicted",
       xlim=c(calc.min,calc.max),ylim=resid.lim,
       xlab="Predicted",ylab="Normalized Residuals",pch=marker,cex=0.8)
  lines(c(calc.min,calc.max),c(0,0))
  res.lm <- lm(resid.norm.abs ~ calc)
  x.new <- c(calc.min,calc.max)
  y.new <- predict(res.lm,newdata=data.frame(calc=x.new))
  lines(x.new, y.new)
  lines(x.new,-y.new)
  
  # plot histogram of residuals:
  # ----------------------------
  
  hist(resid.norm,freq=FALSE,main="Hist. of Residuals",
       xlim=resid.lim,
       xlab="Normalized Residuals",ylab="Density")
  lines(seq(-3,3,by=0.1),dnorm(seq(-3,3,by=0.1)))
  lines(resid.lim,c(0,0))
  
  # normal quantile plot:
  # ---------------------
  
  lim <- max(resid.max,qnorm(1-0.5/length(obs))+0.1)
  qqnorm(resid.norm,main="Sample vs. Normal Quant.",
         xlab="Normal Quantiles",ylab="Sample Quantiles",
         pch=marker,cex=0.8,
         xlim=1.1*lim*c(-1,1),ylim=1.1*lim*c(-1,1))
  lines(c(-10,10),c(-10,10))
  
  # reset plot attributes:
  # ----------------------
  
  mtext(header,side=3,outer=T,adj=0.5,cex=1.2)
  
  par(par.def)
}

################################################################################

# sysanal.resid.diag.boxcox
# =========================

# purpose:
# plot residual diagnostics plots for given Box-Cox transformation parameters

# arguments:
# obs:         vector of observations
# calc:        vector of calculated results
# lambda1:     value of parameter lambda1 (exponent) of Box-Cox transformation
# lambda2:     value of parameter lambda2 (offset) of Box-Cox transformation

# output:
# diagnostic plots

sysanal.resid.diag.boxcox <- function(obs,calc,lambda1=1,lambda2=1)
{
  sysanal.resid.diag(sysanal.boxcox(obs,lambda1,lambda2),
                     sysanal.boxcox(calc,lambda1,lambda2))
}

################################################################################

# sysanal.hessian
# ===============

sysanal.hessian <- function(fn,par,par.inc=NA,...)
{
  if ( is.na(par.inc[1]) ) par.inc <- 0.01*par
  n <- length(par)
  h <- matrix(NA,nrow=n,ncol=n)
  colnames(h) <- names(par)
  rownames(h) <- names(par)
  trace <- matrix(NA,nrow=2*n^2+1,ncol=n+1)
  colnames(trace) <- c(names(par),"fn")
  minimum <- TRUE
  maximum <- TRUE
  par.0 <- par
  res.0 <- fn(par.0,...)
  num.eval <- 1
  trace[1,] <- c(par.0,res.0)
  row.trace <- 1
  for ( i in 1:n )
  {
    for ( j in 1:i )
    {
      if ( i==j )
      {
        par.inc.tmp <- par.inc
        counter <- 0
        while (TRUE)
        {
          par.1 <- par.0
          par.1[i] <- par.1[i] + par.inc.tmp[i]
          res.1 <- fn(par.1,...)
          num.eval <- num.eval+1
          if ( !is.na(res.1) )
          { 
            if ( res.1 > res.0 ) maximum <- FALSE
            if ( res.1 < res.0 ) minimum <- FALSE
            par.2 <- par.0
            par.2[i] <- par.2[i] - par.inc.tmp[i]
            res.2 <- fn(par.2,...)
            num.eval <- num.eval+1
            if ( !is.na(res.2) )
            {
              if ( res.2 > res.0 ) maximum <- FALSE
              if ( res.2 < res.0 ) minimum <- FALSE
              h[i,i] <- (res.1 - 2*res.0 + res.2)/par.inc.tmp[i]^2
              trace[row.trace+1,] <- c(par.1,res.1)
              trace[row.trace+2,] <- c(par.2,res.2)
              row.trace <- row.trace + 2
              break
            }
          }
          counter <- counter + 1
          if ( counter > 5 ) stop("sysanal.hessian: unable to calculate hessian")
          par.inc.tmp <- 0.5*par.inc.tmp
        }
      }
      else
      {
        par.inc.tmp <- par.inc
        counter <- 0
        while (TRUE)
        {
          par.1 <- par.0
          par.1[i] <- par.1[i] + par.inc.tmp[i]
          par.1[j] <- par.1[j] + par.inc.tmp[j]
          res.1 <- fn(par.1,...)
          num.eval <- num.eval + 1
          if ( !is.na(res.1) )
          {
            if ( res.1 > res.0 ) maximum <- FALSE
            if ( res.1 < res.0 ) minimum <- FALSE
            par.2 <- par.0
            par.2[i] <- par.2[i] + par.inc.tmp[i]
            par.2[j] <- par.2[j] - par.inc.tmp[j]
            res.2 <- fn(par.2,...)
            num.eval <- num.eval + 1
            if ( !is.na(res.2) )
            {
              if ( res.2 > res.0 ) maximum <- FALSE
              if ( res.2 < res.0 ) minimum <- FALSE
              par.3 <- par.0
              par.3[i] <- par.3[i] - par.inc.tmp[i]
              par.3[j] <- par.3[j] + par.inc.tmp[j]
              res.3 <- fn(par.3,...)
              num.eval <- num.eval + 1
              if ( !is.na(res.3) )
              {
                if ( res.3 > res.0 ) maximum <- FALSE
                if ( res.3 < res.0 ) minimum <- FALSE
                par.4 <- par.0
                par.4[i] <- par.4[i] - par.inc.tmp[i]
                par.4[j] <- par.4[j] - par.inc.tmp[j]
                res.4 <- fn(par.4,...)
                num.eval <- num.eval + 1
                if ( !is.na(res.4) )
                {
                  if ( res.4 > res.0 ) maximum <- FALSE
                  if ( res.4 < res.0 ) minimum <- FALSE
                  h[i,j] <- (res.1 - res.2 - res.3 + res.4)/
                    (4*par.inc.tmp[i]*par.inc.tmp[j])
                  h[j,i] <- h[i,j]
                  trace[row.trace+1,] <- c(par.1,res.1)
                  trace[row.trace+2,] <- c(par.2,res.2)
                  trace[row.trace+3,] <- c(par.3,res.3)
                  trace[row.trace+4,] <- c(par.4,res.4)
                  row.trace <- row.trace + 4
                  break
                }
              }
            }
          }
          counter <- counter + 1
          if ( counter > 5 ) stop("sysanal.hessian: unable to calculate hessian")
          par.inc.tmp <- 0.5*par.inc.tmp
        }
      } 
    }
  }
  return(list(h        = h,
              minimum  = minimum,
              maximum  = maximum,
              num.eval = num.eval,
              trace    = trace))
}      

##############################################################################
#                                                                            #
# sysanal.smooth                                                             #
# --------------                                                             #
#                                                                            #
# Function for smoothing data points and estimating the derivative of the    #
# smoothed curve by local quadratic or optionally local linear regression.   #
# Local regression is implemented by using a Gaussian distribution of        #
# of weights centered a the point at which the smoothed curve has to be      #
# evaluated.                                                                 #
# The smoothing parameter is the standard deviation of the Gaussian weights. #
#                                                                            #
# Peter Reichert 05.09.2008 , last modification 01.05.2009                   #
#                                                                            #
##############################################################################

# Call:      sysanal.smooth(x,y,sigma,newx=NA,fac.extrap=1,method="quadratic")
# -----
#
# Input:
# ------
#
# x          vector of x-coordinates of data points to be smoothed
# y          vector of y-coordinates of data points to be smoothed
#            (x and y must be of the same length)
# sigma      standard deviation of Gaussian distribution used as weights for
#            local quadratic or local linear regression
# newx       optional vector of x-coordinates at which smoothed results and
#            derivatives are to be calculated (if not specified, results
#            are provided at the same locations as there is data available)
# fac.extrap calculate smoothed value only if data is available within
#            fac.extrap*sigma (default is 1)
# method     "quadratic" (default) indicates local quadratic regression, 
#            otherwise local linear regression is used
#
# Output:
# -------
#
# Data frame of
# x          x-coordinates at which smoothed results and derivatives are 
#            available
# y          smoothed results at the locations x
# ydot       derivatives of smoothed results at the locations x

sysanal.smooth <- function(x,y,sigma,newx=NA,fac.extrap=1,method="quadratic")
{
  # calculate x and y vectors for available data:
  ind <- !is.na(y)
  x.data <- x[ind]
  y.data <- y[ind]
  
  # check consistency of input:
  if ( length(x.data) < 1 )
  {
    stop("*** error in sysanal.smooth: no data available ***")
  }
  if ( length(x.data) != length(y.data) ) 
  {
    stop("*** error in sysanal.smooth: length of x and y different ***")
  }
  if ( ! sigma > 0 ) 
  {
    stop("*** error in sysanal.smooth: sigma is not positive ***")
  }
  
  # select x values for output:
  if ( is.na(newx[1]) ) newx <- x
  
  # set up ysmooth and ydot vectors:
  n <- length(newx)
  ysmooth <- rep(NA,n)
  ydot    <- rep(NA,n)
  
  # calclate smoothed values and derivatives:
  for ( i in 1:n )
  {
    # get indices of data within +/- 2*sigma:
    ind.extrap <- x.data >= newx[i]-fac.extrap*sigma & 
      x.data <= newx[i]+fac.extrap*sigma
    num.extrap <- sum(ifelse(ind.extrap,1,0))
    
    # calculate smoothed value only if data is available within +/- 2*sigma:
    if ( num.extrap > 0 ) 
    {
      # still use data within a 5 times larger interval
      # to calculate the smoothed value:
      fac.use <- 4*max(1,fac.extrap)
      ind.use <- x.data >= newx[i]-fac.use*sigma & 
        x.data <= newx[i]+fac.use*sigma
      x1  <- x.data[ind.use]-newx[i]
      x2  <- (x.data[ind.use]-newx[i])^2
      num.use <- sum(ifelse(ind.use,1,0))
      if ( num.use == 1 )  # use value
      {
        ysmooth[i] <- y.data[ind.use][1]
        ydot[i]    <- 0
      }
      else
      {
        if ( num.use == 2 )  # use weighted mean
        {
          weights <- dnorm(x.data[ind.use],mean=newx[i],sd=sigma)
          weights <- weights/sum(weights)
          ysmooth[i] <- weights[1]*y.data[ind.use][1] + 
            weights[2]*y.data[ind.use][2]
          if ( x.data[ind.use][2] != x.data[ind.use][1] )
          {
            ydot[i]    <- (y.data[ind.use][2]-y.data[ind.use][1])/
              (x.data[ind.use][2]-x.data[ind.use][1])
          }
        }
        else
        {
          if ( method != "quadratic" | num.use == 3 ) # use local linear
          {                                           # regression
            res.lm     <- lm(y.data[ind.use] ~ x1,
                             weights=dnorm(x.data[ind.use],
                                           mean=newx[i],sd=sigma))
            ysmooth[i] <- coef(res.lm)[1]
            ydot[i]    <- coef(res.lm)[2]
          }
          else  # use local quadratic regression
          {
            res.lm     <- lm(y.data[ind.use] ~ x1 + x2,
                             weights=dnorm(x.data[ind.use],
                                           mean=newx[i],sd=sigma))
            ysmooth[i] <- coef(res.lm)[1]
            ydot[i]    <- coef(res.lm)[2]
          }
        }
      }
    }
  }
  
  # return data frame:
  return(data.frame(x=newx,y=ysmooth,ydot=ydot))
}


##############################################################################
#                                                                            #
# sysanal.smooth_fun                                                         #
# ------------------                                                         #
#                                                                            #
# Function for smoothing piecewise linear functions.                         #
# The smoothing parameter is the standard deviation of the Gaussian weights. #
#                                                                            #
# Peter Reichert 16.12.2008 , last modification 01.05.2009                   #
#                                                                            #
##############################################################################

# Call:      sysanal.smooth_fun(data,z,sigma,newz=NA,newx=NA,fac.extrap=1,
#                               method="quadratic")
# -----
#
# Input:
# ------
#
# data       list of matrices specifying piecewise linear functions:
#            for each function, the independent variable x must be provided 
#            in the first column, the dependent variable y in the second
# z          vector of z-coordinates corresponding to the functions
# sigma      standard deviation of Gaussian distribution used as weights for
#            local quadratic regression in z-coordinates (see sysanal.smooth)
# newz       optional vector of z-coordinates at which smoothed functions
#            are to be calculated (if not specified, results are provided
#            at the same locations as there is data available)
# newx       optional vector of x-coordinates  at which smoothed function 
#            values are to be calculated
# fac.extrap calculate smoothed value only if data is available within
#            fac.extrap*sigma (default is 1)
# method     "quadratic" (default) indicates local quadratic regression, 
#            otherwise local linear regression is used
#
# Output:
# -------
#
# List of data frames with columns
# x          x-coordinates at which smoothed results are available
# y          smoothed results at the locations x

sysanal.smooth_fun <- function(data,z,sigma,newz=NA,newx=NA,fac.extrap=1,
                               method="quadratic")
{
  # check consistency of input:
  n <- length(data)
  if ( length(z) != n )
  { 
    stop("error in sysanal.smooth_fun: not same number of locations as functions")
  }
  
  # select z and x values for output:
  if ( is.na(newz[1]) ) newz <- z
  if ( is.na(newx[1]) )
  {
    x <- numeric(0)
    for ( j in 1:n ) x <- c(x,data[[j]][,1])
    newx <- sort(unique(x))
  }
  
  # interpolate input functions to selected x values:    
  y <- matrix(nrow=length(newx),ncol=n)
  for ( j in 1:n )
  {
    y[,j] <- approx(x=data[[j]][,1],y=data[[j]][,2],xout=newx)$y
  }
  
  # set up result data structure:
  res <- list()
  for ( j in 1:length(newz) )
  {
    res[[j]] <- data.frame(x=newx,y=rep(NA,length(newx)))
  }
  names(res) <- newz
  
  # calculate results:
  for ( i in 1:length(newx) )
  {
    res.smooth <- sysanal.smooth(x=z,y=y[i,],sigma=sigma,newx=newz,
                                 fac.extrap=fac.extrap,method=method)$y
    for ( j in 1:length(newz) )
    {
      res[[j]]$y[i] <- res.smooth[j] 
    }     
  }
  
  # return results:
  return(res)
}

############################################################################

# ======================================================== #
# Numerical Integration of Ordinary Differential Equations #
# ======================================================== #


# This library contains a didactical implementation of a numerical integrator
# of a system of ordinary differential equations.
# Note that this implementation is intended to demonstrate how such 
# techniques can be implemented in R. It does not represent the state of the 
# art of numerical integration of such differential equations.


# created and maintained by 
# Peter Reichert
# EAWAG
# Duebendorf
# Switzerland
# reichert@eawag.ch


# First version: Dec.  22, 2002
# Last revision: April 09, 2006


# Overview of functions
# =====================

# sysanal.ode:    numerical integration of deterministic ordinary differential 
#                 equations


# =========================================================================== #


# numerical integration of ordinary differential equations
# ========================================================

sysanal.ode <- function(rhs,x.ini,par,t.out,dt.max=NA,algorithm="euler",...)
{
  # -----------------------------------------------------------------------
  # This solver for a system of ordindary differential equations
  # was written for didactical purposes. It does not represent 
  # the state of the art in numerical integration techniques
  # but should rather demonstrate how the simplest integration
  # technique can be implemented for relatively simple use.
  # Check the package "odesolve" available from http://www.r-project.org
  # for professional solvers for ordinary differential equations.
  #
  # Arguments:
  # rhs:       function returning the right hand side of the 
  #            system of differential equations as a function
  #            of the arguments x (state variables), t (time),
  #            and par (model parameters).
  # x.ini:     start values of the state variables.
  # par:       model parameters (transferred to rhs).
  # t.out:     set of points in time at which the solution 
  #            should be calculated; the solution at the first
  #            value in t is set to x.ini, the solution at subsequent
  #            values is calculated by numerical integration.
  # dt.max:    maximum time step; if the difference between points
  #            in time at which output has to be provided, t, is
  #            larger than dt.max, then this output interval is 
  #            divided into step of at most dt.max internally to
  #            improve the accuracy of the solution at the next 
  #            output point.
  # algorithm: right now, the only options are "euler" (explicit first order
  #            Euler algorithm) and "euler.2.order" (explicit second order
  #            Euler algorithm).
  #
  # Return Value:
  # matrix with results for state variables (columns) at all output
  # time steps t.out (rows).
  #
  #                                        Peter Reichert    Dec.  22, 2002
  #                                        last modification April 09, 2006
  # -----------------------------------------------------------------------
  
  # determine number of equations and number of time steps:
  num.eq <- length(x.ini)
  steps <- length(t.out)
  
  # define and initialize result matrix:
  x <- matrix(nrow=steps,ncol=num.eq)
  colnames(x) <- names(x.ini)
  rownames(x) <- t.out
  x[1,] <- x.ini
  
  # perform integration:
  if ( algorithm == "euler" )
  {
    for ( i in 2:steps )
    {
      if ( is.na(dt.max) || dt.max >= t.out[i]-t.out[i-1] )
      {
        # output interval <= dt.max (or dt.max not available):
        # perform a single step to the next output time:
        x[i,] <- x[i-1,] + (t.out[i]-t.out[i-1])*rhs(x[i-1,],t.out[i-1],par,...)
      }
      else
      {
        # output interval > dt.max:
        # perform multiple steps of maximum size dt.max:
        x[i,] <- x[i-1,]
        steps.internal <- ceiling((t.out[i]-t.out[i-1])/dt.max)
        dt <- (t.out[i]-t.out[i-1])/steps.internal
        for ( j in 1:steps.internal )
        {
          t.current <- t.out[i-1] + (j-1)*dt
          x[i,] <- x[i,] + dt*rhs(x[i,],t.current,par,...)
        }
      }
    }
  }
  else
  {
    if ( algorithm == "euler.2.order" )
    {
      for ( i in 2:steps )
      {
        if ( is.na(dt.max) || dt.max >= t.out[i]-t.out[i-1] )
        {
          # output interval <= dt.max (or dt.max not available):
          # perform a single step to the next output time:
          t.mid <- 0.5*(t.out[i-1]+t.out[i])
          x.mid <- x[i-1,] + 0.5*(t.out[i]-t.out[i-1])*rhs(x[i-1,],t.out[i-1],par,...)
          x[i,] <- x[i-1,] + (t.out[i]-t.out[i-1])*rhs(x.mid,t.mid,par,...)
        }
        else
        {
          # output interval > dt.max:
          # perform multiple steps of maximum size dt.max:
          x[i,] <- x[i-1,]
          steps.internal <- ceiling((t.out[i]-t.out[i-1])/dt.max)
          dt <- (t.out[i]-t.out[i-1])/steps.internal
          for ( j in 1:steps.internal )
          {
            t.current <- t.out[i-1] + (j-1)*dt
            t.mid <- t.current + 0.5*dt
            x.mid <- x[i,] + 0.5*dt*rhs(x[i,],t.current,par,...)
            x[i,] <- x[i,] + dt*rhs(x.mid,t.mid,par,...)
          }
        }
      }
    }
    else
    {
      stop(paste("sysanal.ode: algorithm \"",
                 algorithm,
                 "\" not implemented",
                 sep=""))
    }
  }
  
  # return result matrix: 
  return(x)
}


# =========================================================================== #


sysanal.Var.B.M.L <- function(psi,L,i=NA,j=NA)
{
  L.decoded <- L
  L.encoded <- L
  if ( is.vector(L) )
  {
    L.decoded <- sysanal.decode(L)
  }
  else
  {
    L.encoded <- rownames(L)
  }
  n <- length(L.encoded)
  
  i.local <- i
  if ( is.na(i.local[1]) ) i.local <- 1:n
  j.local <- j
  if ( is.na(j.local[1]) ) j.local <- 1:n
  Var <- matrix(0,nrow=length(i.local),ncol=length(j.local))
  vars <- unique(L.decoded$var)
  beta <- 1/psi["corrlen"]^2
  if ( is.na(beta) ) beta <- 0
  for ( var in vars )
  {
    ind <- which(L.decoded$var==var)
    ind.i <- match(ind,i.local)
    ind.i <- ind.i[!is.na(ind.i)]
    ind.j <- match(ind,j.local)
    ind.j <- ind.j[!is.na(ind.j)]
    if ( length(ind.i)>0 & length(ind.j)>0 )
    {
      name.sd.B <- paste("sd.B",var,sep="_")
      var.B <- psi[name.sd.B]^2
      if ( is.na(var.B) ) var.B <- 0
      dist <- rep(1,length(ind.i)) %*% t(L.decoded$val[j.local[ind.j]]) -
        L.decoded$val[i.local[ind.i]] %*% t(rep(1,length(ind.j)))
      
      Var[ind.i,ind.j] <- var.B*exp(-beta*dist^2)
    }
  }
  rownames(Var) <- L.encoded[i.local]
  colnames(Var) <- L.encoded[j.local]
  return(Var)
}


sysanal.sd.Eps.L <- function(xi,L)
{
  L.decoded <- L
  L.encoded <- L
  if ( is.vector(L) )
  {
    L.decoded <- sysanal.decode(L)
  }
  else
  {
    L.encoded <- rownames(L)
  }
  n <- length(L.encoded)
  
  sd <- rep(NA,n)
  vars <- unique(L.decoded$var)
  for ( var in vars )
  {
    name.sd.Eps <- paste("sd.Eps",var,sep="_")
    sd.Eps <- xi[name.sd.Eps]
    if ( is.na(sd.Eps) ) stop(paste("*** parameter",
                                    name.sd.Eps,
                                    "not found in sd.Eps.L"))
    ind <- which(L.decoded$var==var)
    sd[ind] <- rep(sd.Eps,length(ind))
  }
  names(sd) <- L.encoded   
  return(sd)
}


sysanal.loglikeli.bias <- function(par,model,L,y.obs,Var.B,sd.Eps,
                                   lambda1=1,lambda2=1,...)
{
  # decode layout definition: 
  
  L.decoded <- L
  L.encoded <- L
  if ( is.vector(L) )
  {
    L.decoded <- sysanal.decode(L)
  }
  else
  {
    L.encoded <- rownames(L)
  }
  
  # calculate results of deterministic model:
  
  y.calc        <- model(par,L.decoded,...)
  
  # transform results and observations:
  
  y.calc.trans  <- sysanal.boxcox(y.calc,lambda1,lambda2)
  y.obs.trans   <- sysanal.boxcox(y.obs,lambda1,lambda2)
  boxcox.deriv  <- sysanal.boxcox.deriv(y.obs,lambda1,lambda2)
  
  # evaluate likelihood function:
  
  Sigma.B       <- Var.B(par,L.decoded)
  Sigma.Eps     <- diag(sd.Eps(par,L.decoded)^2)
  
  Sum.Sigma     <- Sigma.B + Sigma.Eps
  
  Sum.Sigma.inv <- solve(Sum.Sigma)
  
  log.det.Sum.Sigma <- determinant(Sum.Sigma,logarithm=TRUE)
  if ( log.det.Sum.Sigma$sign < 0 ) stop("determinant Sigma.Eps+Sigma.B < 0")
  
  loglikeli <- - 0.5 * length(L.encoded) * log(2*pi) -
    0.5 * as.numeric(log.det.Sum.Sigma$modulus) -
    0.5 * t(y.obs.trans-y.calc.trans) %*% (Sum.Sigma.inv ) %*% 
    (y.obs.trans-y.calc.trans) +
    sum(log(abs(boxcox.deriv)))
  
  return(loglikeli)
}


sysanal.loglikeli.bias.blockdesign <- function(par,
                                               model,
                                               L,
                                               y.obs,
                                               Var.B,
                                               sd.Eps,
                                               lambda1 = 1,
                                               lambda2 = 1,
                                               ...)
{
  # decode layout definition:
  
  L.decoded <- L
  L.encoded <- L
  if ( is.vector(L) )
  {
    L.decoded <- sysanal.decode(L)
  }
  else
  {
    L.encoded <- rownames(L)
  }
  n <- length(L.encoded)
  
  # calculate results of deterministic model:
  
  y.calc <- model(par,L.decoded,...)
  
  # transform results and observations:
  
  y.calc.trans  <- sysanal.boxcox(y.calc,lambda1,lambda2)
  y.obs.trans   <- sysanal.boxcox(y.obs,lambda1,lambda2)
  boxcox.deriv  <- sysanal.boxcox.deriv(y.obs,lambda1,lambda2)
  
  # evaluate likelihood function:
  
  loglikeli <- 0
  vars <- unique(L.decoded$var)
  for ( var in vars )
  {
    ind <- which(L.decoded$var==var)
    
    Sigma.B   <- Var.B(par,L.decoded[ind,])
    Sigma.Eps <- diag(sd.Eps(par,L.decoded[ind,])^2)
    
    Sum.Sigma <- Sigma.Eps + Sigma.B
    
    Sum.Sigma.inv <- solve(Sum.Sigma)
    
    log.det.Sum.Sigma <- determinant(Sum.Sigma,logarithm=TRUE)
    if ( log.det.Sum.Sigma$sign < 0 ) stop("determinant Sigma.Eps+Sigma.B < 0")
    
    loglikeli <- loglikeli -
      0.5 * length(ind) * log(2*pi) -
      0.5 * as.numeric(log.det.Sum.Sigma$modulus) -
      0.5 * t(y.obs.trans[ind]-y.calc.trans[ind]) %*%
      Sum.Sigma.inv %*%
      (y.obs.trans[ind]-y.calc.trans[ind]) +
      sum(log(abs(boxcox.deriv[ind])))
  }
  
  return(loglikeli)
}


# define prediction functions:
# ----------------------------

sysanal.predict.bias.cond <- function(par,model,L1,y.obs,
                                      Var.B,sd.Eps,
                                      L2=NA,y.calc=NA,
                                      lambda1=1,lambda2=1,
                                      ...)
{
  # decode likelihood definitions:
  
  L1.decoded <- L1
  L1.encoded <- L1
  if ( is.vector(L1) )
  {
    L1.decoded <- sysanal.decode(L1)
  }
  else
  {
    L1.encoded <- rownames(L1)
  }
  n1 <- length(L1.encoded)
  
  L.decoded <- L1.decoded
  L.encoded <- L1.encoded
  n2 <- 0
  n <- n1
  
  L2.available <- TRUE
  if ( is.vector(L2) ) { if ( is.na(L2[1]) ) L2.available <- FALSE }   
  if ( L2.available )
  { 
    L2.decoded <- L2
    L2.encoded <- L2
    if ( is.vector(L2) )
    {
      L2.decoded <- sysanal.decode(L2)
    }
    else
    {
      L2.encoded <- rownames(L2)
    }
    n2 <- length(L2.encoded)
    
    L.decoded <- rbind(L1.decoded,L2.decoded)
    L.encoded <- c(L1.encoded,L2.encoded)
    n <- n1 + n2
  }
  
  # calculate results of deterministic model:
  
  if ( is.na(y.calc[1]) )
  {   
    y.calc    <- model(par,L.decoded,...)
  }
  else
  {
    if ( length(y.calc) != n )
    {
      cat("*** y.calc is not of correct length:",length(y.calc),
          "instead of",n,"\n")
      return(NA)
    }
  }
  y.calc.L1 <- y.calc[1:n1]
  y.calc.L2 <- NA
  if ( n2 > 0 )
  {
    y.calc.L2    <- y.calc[n1+(1:n2)]
  }
  
  # transform results and observations:
  
  y.calc.L1.trans  <- sysanal.boxcox(y.calc.L1,lambda1,lambda2)
  y.calc.L2.trans  <- sysanal.boxcox(y.calc.L2,lambda1,lambda2)
  y.obs.trans      <- sysanal.boxcox(y.obs,lambda1,lambda2)
  
  # calculate predictions for layout 1:   
  
  Sigma.B.L1   <- Var.B(par,L1.decoded)
  Sigma.Eps.L1 <- diag(sd.Eps(par,L1.decoded)^2)
  
  Sum.Sigma.L1 <- Sigma.B.L1 + Sigma.Eps.L1
  
  Sum.Sigma.L1.inv <- solve(Sum.Sigma.L1)
  
  B.var.L1  <- diag(Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% Sigma.Eps.L1)
  B.mean.L1 <- as.numeric(Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% (y.obs.trans-y.calc.L1.trans))
  names(B.mean.L1) <- L1.encoded
  names(B.var.L1)  <- L1.encoded
  
  Y.mean.L1 <- y.calc.L1.trans + 
    Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% ( y.obs.trans - y.calc.L1.trans )
  names(Y.mean.L1) <- L1.encoded
  Y.var.L1  <- diag(Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% Sigma.Eps.L1) +
    diag(Sigma.Eps.L1) 
  names(Y.var.L1)  <- L1.encoded
  
  # calculate predictions for layout 2:
  
  Y.mean.L2 <- rep(NA,max(1,n2))
  Y.var.L2  <- rep(NA,max(1,n2))
  B.mean.L2 <- rep(NA,max(1,n2))
  B.var.L2  <- rep(NA,max(1,n2))
  Sigma.B.L2L1 <- matrix(NA,nrow=n2,ncol=n1,dimnames=list(L2.encoded,L1.encoded))
  if ( n2 > 0 )
  {
    for ( i in 1:n2 )
    {
      v <- as.vector(Var.B(par,L.decoded,i=1:n1,j=n1+i))
      Sigma.B.L2L1[i,] <- v
      Sigma.B.L2.i.i   <- as.numeric(Var.B(par,L2.decoded,i,i))
      
      Y.mean.L2[i] <- y.calc.L2.trans[i] + t(v) %*% Sum.Sigma.L1.inv %*% 
        ( y.obs.trans - y.calc.L1.trans )
      Y.var.L2[i]  <- Sigma.B.L2.i.i + 
        sd.Eps(par,L2.decoded[i,])^2 - 
        t(v) %*% Sum.Sigma.L1.inv %*% v 
      B.mean.L2[i] <- t(v) %*% Sum.Sigma.L1.inv %*% 
        ( y.obs.trans - y.calc.L1.trans )
      B.var.L2[i]  <- Sigma.B.L2.i.i - 
        t(v) %*% Sum.Sigma.L1.inv %*% v 
    }
    names(Y.mean.L2) <- L2.encoded
    names(Y.var.L2)  <- L2.encoded
    names(B.mean.L2) <- L2.encoded
    names(B.var.L2)  <- L2.encoded
  }
  
  return(list(y.calc.L1        = y.calc.L1.trans,
              B.mean.L1        = B.mean.L1,
              B.var.L1         = B.var.L1,
              Y.mean.L1        = Y.mean.L1,
              Y.var.L1         = Y.var.L1,
              y.calc.L2        = y.calc.L2.trans,
              B.mean.L2        = B.mean.L2,
              B.var.L2         = B.var.L2,
              Y.mean.L2        = Y.mean.L2,
              Y.var.L2         = Y.var.L2,
              Sum.Sigma.L1.inv = Sum.Sigma.L1.inv,
              Sigma.B.L1       = Sigma.B.L1,
              Sigma.B.L2L1     = Sigma.B.L2L1))
}


sysanal.predict.bias.cond.blockdesign <- function(par,model,L1,y.obs,
                                                  Var.B,sd.Eps,
                                                  L2=NA,y.calc=NA,
                                                  lambda1=1,lambda2=1,
                                                  ...)
{
  # decode likelihood definitions:
  
  L1.decoded <- L1
  L1.encoded <- L1
  if ( is.vector(L1) )
  {
    L1.decoded <- sysanal.decode(L1)
  }
  else
  {
    L1.encoded <- rownames(L1)
  }
  n1 <- length(L1.encoded)
  
  L.decoded <- L1.decoded
  L.encoded <- L1.encoded
  n2 <- 0
  n <- n1
  
  L2.available <- TRUE
  if ( is.vector(L2) ) { if ( is.na(L2[1]) ) L2.available <- FALSE }   
  if ( L2.available )
  { 
    L2.decoded <- L2
    L2.encoded <- L2
    if ( is.vector(L2) )
    {
      L2.decoded <- sysanal.decode(L2)
    }
    else
    {
      L2.encoded <- rownames(L2)
    }
    n2 <- length(L2.encoded)
    
    L.decoded <- rbind(L1.decoded,L2.decoded)
    L.encoded <- c(L1.encoded,L2.encoded)
    n <- n1 + n2
  }
  
  # calculate results of deterministic model:
  
  if ( is.na(y.calc[1]) )
  {   
    y.calc    <- model(par,L.decoded,...)
  }
  else
  {
    if ( length(y.calc) != n )
    {
      cat("*** y.calc is not of correct length:",length(y.calc),
          "instead of",n,"\n")
      return(NA)
    }
  }
  y.calc.L1 <- y.calc[1:n1]
  y.calc.L2 <- NA
  if ( n2 > 0 )
  {
    y.calc.L2    <- y.calc[n1+(1:n2)]
  }
  
  # transform results and observations:
  
  y.calc.L1.trans  <- sysanal.boxcox(y.calc.L1,lambda1,lambda2)
  y.calc.L2.trans  <- sysanal.boxcox(y.calc.L2,lambda1,lambda2)
  y.obs.trans      <- sysanal.boxcox(y.obs,lambda1,lambda2)
  
  # set up arrays and calculate results:
  
  B.mean.L1 <- rep(NA,n1)
  B.var.L1  <- rep(NA,n1)
  Y.mean.L1 <- rep(NA,n1)
  Y.var.L1  <- rep(NA,n1)
  B.mean.L2 <- rep(NA,max(1,n2))
  B.var.L2  <- rep(NA,max(1,n2))
  Y.mean.L2 <- rep(NA,max(1,n2))
  Y.var.L2  <- rep(NA,max(1,n2))
  list.Sum.Sigma.L1.inv <- list()
  list.Sigma.B.L1       <- list()
  list.Sigma.B.L2L1     <- list()
  vars <- unique(L.decoded$var)
  for ( var in vars )
  {
    ind1 <- which(L1.decoded$var==var)
    
    Sigma.B.L1   <- Var.B(par,L1.decoded[ind1,])
    sigma.Eps.L1 <- sd.Eps(par,L1.decoded[ind1,])
    Sigma.Eps.L1 <- diag(sigma.Eps.L1^2)
    
    Sum.Sigma.L1 <- Sigma.B.L1 + Sigma.Eps.L1
    
    Sum.Sigma.L1.inv <- solve(Sum.Sigma.L1)
    
    list.Sigma.B.L1[[var]]       <- Sigma.B.L1
    list.Sum.Sigma.L1.inv[[var]] <- Sum.Sigma.L1.inv
    
    B.var.L1[ind1]  <- diag(Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% Sigma.Eps.L1)
    B.mean.L1[ind1] <- Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% 
      (y.obs.trans[ind1]-y.calc.L1.trans[ind1])
    
    # calculate predictions for layout 1:
    
    Y.mean.L1[ind1] <- y.calc.L1.trans[ind1] + 
      Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% 
      ( y.obs.trans[ind1] - y.calc.L1.trans[ind1] )
    Y.var.L1[ind1]  <- diag(Sigma.B.L1 %*% Sum.Sigma.L1.inv %*% Sigma.Eps.L1) +
      diag(Sigma.Eps.L1)
    
    # calculate predictions for layout 2: 
    
    if ( n2 > 0 )
    {
      ind2 <- which(L2.decoded$var==var)
      if ( length(ind2) > 0 )
      {
        list.Sigma.B.L2L1[[var]] <- matrix(nrow=length(ind2),ncol=length(ind1))
        
        for ( i in 1:length(ind2) )
        {
          v <- as.vector(Var.B(par,L.decoded,i=ind1,j=n1+ind2[i]))
          list.Sigma.B.L2L1[[var]][i,] <- v
          Sigma.B.L2.i.i   <- as.numeric(Var.B(par,L2.decoded,ind2[i],ind2[i]))
          
          Y.mean.L2[ind2[i]] <- y.calc.L2.trans[ind2[i]] +  
            t(v) %*% Sum.Sigma.L1.inv %*% 
            ( y.obs.trans[ind1] - y.calc.L1.trans[ind1] )
          Y.var.L2[ind2[i]]  <- Sigma.B.L2.i.i + 
            sd.Eps(par,L2.decoded[ind2[i],])^2 - 
            t(v) %*% Sum.Sigma.L1.inv %*% v 
          
          B.mean.L2[ind2[i]] <- t(v) %*% Sum.Sigma.L1.inv %*% 
            ( y.obs.trans[ind1] - y.calc.L1.trans[ind1] )
          B.var.L2[ind2[i]]  <- Sigma.B.L2.i.i - 
            t(v) %*% Sum.Sigma.L1.inv %*% v 
        }
      }
    }
  }
  names(B.mean.L1) <- L1.encoded
  names(B.var.L1)  <- L1.encoded
  names(Y.mean.L1) <- L1.encoded
  names(Y.var.L1)  <- L1.encoded
  names(B.mean.L2) <- L2.encoded
  names(B.var.L2)  <- L2.encoded
  names(Y.mean.L2) <- L2.encoded
  names(Y.var.L2)  <- L2.encoded
  
  return(list(y.calc.L1        = y.calc.L1.trans,
              B.mean.L1        = B.mean.L1,
              B.var.L1         = B.var.L1,
              Y.mean.L1        = Y.mean.L1,
              Y.var.L1         = Y.var.L1,
              y.calc.L2        = y.calc.L2.trans,
              B.mean.L2        = B.mean.L2,
              B.var.L2         = B.var.L2,
              Y.mean.L2        = Y.mean.L2,
              Y.var.L2         = Y.var.L2,
              Sum.Sigma.L1.inv = list.Sum.Sigma.L1.inv,
              Sigma.B.L1       = list.Sigma.B.L1,
              Sigma.B.L2L1     = list.Sigma.B.L2L1))
}


sysanal.predict.bias <- function(parsamp.L1,model,L1,y.obs,
                                 predict.bias.cond,
                                 probs=c(0.025,0.05,0.25,0.5,0.75,0.95,0.975),
                                 L2=NA,y.calc.samp=NA,
                                 lambda1=1,lambda2=1,
                                 ...)
{
  # decode likelihood definitions:
  
  L1.decoded <- L1
  L1.encoded <- L1
  if ( is.vector(L1) )
  {
    L1.decoded <- sysanal.decode(L1)
  }
  else
  {
    L1.encoded <- rownames(L1)
  }
  n1 <- length(L1.encoded)
  
  L.decoded <- L1.decoded
  L.encoded <- L1.encoded
  n2 <- 0
  n <- n1
  
  L2.available <- TRUE
  if ( is.vector(L2) ) { if ( is.na(L2[1]) ) L2.available <- FALSE }   
  if ( L2.available )
  { 
    L2.decoded <- L2
    L2.encoded <- L2
    if ( is.vector(L2) )
    {
      L2.decoded <- sysanal.decode(L2)
    }
    else
    {
      L2.encoded <- rownames(L2)
    }
    n2 <- length(L2.encoded)
    
    L.decoded <- rbind(L1.decoded,L2.decoded)
    L.encoded <- c(L1.encoded,L2.encoded)
    n <- n1 + n2
  }
  
  # initialize result arrays:
  
  y.L1.samp <- matrix(nrow=nrow(parsamp.L1),ncol=n1)
  colnames(y.L1.samp) <- L1.encoded
  B.L1.samp <- y.L1.samp
  Y.L1.samp <- y.L1.samp
  neg.var.B.L1 <- rep(0,n1)
  names(neg.var.B.L1) <- L1.encoded
  neg.var.Y.L1 <- neg.var.B.L1
  
  y.L2.samp    <- NA
  B.L2.samp    <- NA
  Y.L2.samp    <- NA
  neg.var.B.L2 <- NA
  neg.var.Y.L2 <- NA
  if ( n2 > 0 )
  {
    y.L2.samp <- matrix(nrow=nrow(parsamp.L1),ncol=n2)
    colnames(y.L2.samp) <- L2.encoded
    B.L2.samp <- y.L2.samp
    Y.L2.samp <- y.L2.samp
    neg.var.B.L2 <- rep(0,n2)
    names(neg.var.B.L2) <- L2.encoded
    neg.var.Y.L2 <- neg.var.B.L2
  }
  
  # calculate samples:
  
  par.old <- rep(NA,ncol(parsamp.L1))
  num.eval <- 0
  for ( j in 1:nrow(parsamp.L1) )
  {
    par <- parsamp.L1[j,]
    y.calc <- NA
    if ( length(dim(y.calc.samp)) == 2 ) y.calc <- y.calc.samp[j,]
    
    if ( j==1 | sum((par-par.old)^2) != 0 ) # do not reevaluate if parameter
    {                                       # values stayed the same
      par.old <- par
      num.eval <- num.eval + 1      
      res <- predict.bias.cond(par     = par,
                               model   = model,
                               L1      = L1.decoded,
                               y.obs   = y.obs,
                               L2      = L2,
                               y.calc  = y.calc,
                               lambda1 = lambda1,
                               lambda2 = lambda2,
                               ...)
    }
    y.L1.samp[j,] <- res$y.calc.L1
    for ( i in 1:n1 )
    {
      var <- res$Y.var.L1[i]
      if ( var < 0 ) 
      { 
        cat("* Warning: negative variance of Y:",var,"at",L1.encoded[i],"\n") 
        var <- 0; neg.var.Y.L1[i] <- neg.var.Y.L1[i]+1
      }
      Y.L1.samp[j,i] <- rnorm(1,res$Y.mean.L1[i],sqrt(var))
    }
    if ( ! is.na(res$B.mean.L1[1]) )
    {
      for ( i in 1:n1 )
      {   
        var <- res$B.var.L1[i]
        if ( var < 0 ) 
        { 
          cat("* Warning: negative variance of B:",var,"at",L1.encoded[i],"\n") 
          var <- 0; neg.var.B.L1[i] <- neg.var.B.L1[i]+1 
        }
        B.L1.samp[j,i] <- rnorm(1,res$B.mean.L1[i],sqrt(var))   
      }
    }
    if ( n2 > 0 )
    {
      y.L2.samp[j,] <- res$y.calc.L2
      for ( i in 1:n2 )
      {
        var <- res$Y.var.L2[i]
        if ( var < 0 ) 
        { 
          cat("* Warning: negative variance of Y:",var,"at",L2.encoded[i],"\n") 
          var <- 0; neg.var.Y.L2[i] <- neg.var.Y.L2[i]+1 
        }
        Y.L2.samp[j,i] <- rnorm(1,res$Y.mean.L2[i],sqrt(var))
      }
      if ( ! is.na(res$B.mean.L2[1]) )
      {
        for ( i in 1:n2 )
        {   
          var <- res$B.var.L2[i]
          if ( var < 0 ) 
          { 
            cat("* Warning: negative variance of B:",var,"at",L2.encoded[i],"\n") 
            var <- 0; neg.var.B.L2[i] <- neg.var.B.L2[i]+1 
          }
          B.L2.samp[j,i] <- rnorm(1,res$B.mean.L2[i],sqrt(var))
        }
      }      
    }
  }
  neg.var.B.L1 <- neg.var.B.L1/n1
  neg.var.Y.L1 <- neg.var.Y.L1/n1
  if ( n2 > 0 ) 
  {
    neg.var.B.L2 <- neg.var.B.L2/n2
    neg.var.Y.L2 <- neg.var.Y.L2/n2
  }
  
  # derive quantiles:
  
  y.L1.quant <- matrix(nrow=length(probs),ncol=n1)
  colnames(y.L1.quant) <- L1.encoded
  rownames(y.L1.quant) <- probs
  B.L1.quant <- y.L1.quant
  yplusB.L1.quant <- y.L1.quant
  Y.L1.quant <- y.L1.quant
  for ( i in 1:n1 )
  {
    y.L1.quant[,i]      <- quantile(y.L1.samp[,i],probs=probs,na.rm=TRUE)
    B.L1.quant[,i]      <- quantile(B.L1.samp[,i],probs=probs,na.rm=TRUE)
    yplusB.L1.quant[,i] <- quantile(y.L1.samp[,i]+B.L1.samp[,i],
                                    probs=probs,na.rm=TRUE)
    Y.L1.quant[,i]      <- quantile(Y.L1.samp[,i],probs=probs,na.rm=TRUE)
  }
  
  y.L2.quant <- NA
  B.L2.quant <- NA
  yplusB.L2.quant <- NA
  Y.L2.quant <- NA
  if ( n2 > 0 )
  {
    y.L2.quant <- matrix(nrow=length(probs),ncol=n2)
    colnames(y.L2.quant) <- L2
    rownames(y.L2.quant) <- probs
    B.L2.quant <- y.L2.quant
    yplusB.L2.quant <- y.L2.quant
    Y.L2.quant <- y.L2.quant
    for ( i in 1:n2 )
    {
      y.L2.quant[,i]      <- quantile(y.L2.samp[,i],probs=probs,na.rm=TRUE)
      B.L2.quant[,i]      <- quantile(B.L2.samp[,i],probs=probs,na.rm=TRUE)
      yplusB.L2.quant[,i] <- quantile(y.L2.samp[,i]+B.L2.samp[,i],
                                      probs=probs,na.rm=TRUE)
      Y.L2.quant[,i]      <- quantile(Y.L2.samp[,i],probs=probs,na.rm=TRUE)
    }
  }
  
  return(list(y.L1.samp       = y.L1.samp,
              B.L1.samp       = B.L1.samp,
              Y.L1.samp       = Y.L1.samp,
              y.L2.samp       = y.L2.samp,
              B.L2.samp       = B.L2.samp,
              Y.L2.samp       = Y.L2.samp,
              y.L1.quant      = y.L1.quant,
              B.L1.quant      = B.L1.quant,
              yplusB.L1.quant = yplusB.L1.quant,
              Y.L1.quant      = Y.L1.quant,
              y.L2.quant      = y.L2.quant,
              B.L2.quant      = B.L2.quant,
              yplusB.L2.quant = yplusB.L2.quant,
              Y.L2.quant      = Y.L2.quant,
              neg.var.B.L1    = neg.var.B.L1,
              neg.var.Y.L1    = neg.var.Y.L1,
              neg.var.B.L2    = neg.var.B.L2,
              neg.var.Y.L2    = neg.var.Y.L2))
}



sysanal.predict.bias.lin <- function(par.mean.L1,
                                     par.var.L1,
                                     model,
                                     L1,
                                     y.obs,
                                     predict.bias.cond,
                                     par.inc = 0.01*par.mean.L1,
                                     L2      = NA,
                                     y.calc  = NA,
                                     V.y     = NA,
                                     lambda1 = 1,
                                     lambda2 = 1,
                                     logpar  = F,
                                     ...)
{
  # decode likelihood definitions:
  
  L1.decoded <- L1
  L1.encoded <- L1
  if ( is.vector(L1) )
  {
    L1.decoded <- sysanal.decode(L1)
  }
  else
  {
    L1.encoded <- rownames(L1)
  }
  n1 <- length(L1.encoded)
  
  L.decoded <- L1.decoded
  L.encoded <- L1.encoded
  n2 <- 0
  n <- n1
  
  L2.available <- TRUE
  if ( is.vector(L2) ) { if ( is.na(L2[1]) ) L2.available <- FALSE }   
  if ( L2.available )
  { 
    L2.decoded <- L2
    L2.encoded <- L2
    if ( is.vector(L2) )
    {
      L2.decoded <- sysanal.decode(L2)
    }
    else
    {
      L2.encoded <- rownames(L2)
    }
    n2 <- length(L2.encoded)
    
    L.decoded <- rbind(L1.decoded,L2.decoded)
    L.encoded <- c(L1.encoded,L2.encoded)
    n <- n1 + n2
  }
  
  # calculate transformed data:
  
  y.obs.trans <- sysanal.boxcox(y.obs,lambda1,lambda2)
  
  # calculate solution and conditional predictions:
  
  par <- par.mean.L1
  par.local <- par
  if ( logpar ) par.local <- exp(par)
  res.par <- predict.bias.cond(par     = par.local,
                               model   = model,
                               L1      = L1.decoded,
                               y.obs   = y.obs,
                               L2      = L2.decoded,
                               y.calc  = y.calc,
                               lambda1 = lambda1,
                               lambda2 = lambda2,
                               ...)
  boxcox.deriv <- sysanal.boxcox.deriv(sysanal.boxcox.inv(c(res.par$y.calc.L1,
                                                            res.par$y.calc.L2),
                                                          lambda1,lambda2),
                                       lambda1,lambda2)
  
  # calculate jacobian of bias and of model plus bias:
  
  V.yplusB <- matrix(NA,nrow=n,ncol=length(par))
  rownames(V.yplusB) <- L.encoded
  colnames(V.yplusB) <- names(par)
  V.y.test <- V.yplusB
  V.B      <- V.yplusB
  V.Eps    <- matrix(NA,nrow=n1,ncol=length(par))
  rownames(V.Eps) <- L1.encoded
  colnames(V.Eps) <- names(par)
  if ( !is.matrix(V.y) & !is.data.frame(V.y) ) # calculate numerical deriv 
  {                                            # of all expected val.
    for( j in 1:length(par) )
    {
      par.j <- par
      par.j[j] <- par.j[j] + par.inc[j]
      if ( logpar ) par.j <- exp(par.j)
      res.par.j <- predict.bias.cond(par     = par.j,
                                     model   = model,
                                     L1      = L1.decoded,
                                     y.obs   = y.obs,
                                     L2      = L2,
                                     y.calc  = NA,
                                     lambda1 = lambda1,
                                     lambda2 = lambda2,
                                     ...)
      
      V.y.test[,j] <- ( c(res.par.j$y.calc.L1,res.par.j$y.calc.L2) - 
                          c(res.par$y.calc.L1,res.par$y.calc.L2) ) /
        par.inc[j]
      V.yplusB[,j] <- ( c(res.par.j$y.calc.L1 + res.par.j$B.mean.L1,
                          res.par.j$y.calc.L2 + res.par.j$B.mean.L2) - 
                          c(res.par$y.calc.L1   + res.par$B.mean.L1,
                            res.par$y.calc.L2   + res.par$B.mean.L2) ) /
        par.inc[j]
      V.B[,j]      <- ( c(res.par.j$B.mean.L1,res.par.j$B.mean.L2) - 
                          c(res.par$B.mean.L1,res.par$B.mean.L2) ) /
        par.inc[j]
      V.Eps[,j]    <- ( (y.obs.trans - res.par.j$y.calc.L1 - res.par.j$B.mean.L1) -
                          (y.obs.trans - res.par$y.calc.L1 - res.par$B.mean.L1) ) /
        par.inc[j]
    }
    #write.table(V.y.test,"V.y.calculated.dat",row.names=T,col.names=NA,sep="\t")
  }
  else
  {
    # extend Y.y to other parameters, assume no sensitivity
    
    if ( nrow(V.y) != n ) stop("incorrect number of rows in V.y")
    V.y.local <- matrix(0,nrow=n,ncol=length(par))
    colnames(V.y.local) = names(par)
    rownames(V.y.local) = L.encoded
    for ( j in colnames(V.y) )  # copy entries from V.y
    {
      V.y.local[,j] <- V.y[,j]*boxcox.deriv
    }
    #write.table(V.y.local,"V.y.derived.dat",row.names=T,col.names=NA,sep="\t")
    
    # calculate derivatives of matrix factors:
    
    if ( ! is.list(res.par$Sum.Sigma.L1.inv) )
    {
      A <- res.par$Sigma.B.L1   %*% res.par$Sum.Sigma.L1.inv
      if ( n2 > 0 ) B <- res.par$Sigma.B.L2L1 %*% res.par$Sum.Sigma.L1.inv
      for( j in 1:length(par) )
      {
        par.j <- par
        par.j[j] <- par.j[j] + par.inc[j]
        if ( logpar ) par.j <- exp(par.j)
        res.par.j <- predict.bias.cond(par     = par.j,
                                       model   = model,
                                       L1      = L1.decoded,
                                       y.obs   = y.obs,
                                       L2      = L2,
                                       y.calc  = c(res.par$y.calc.L1,res.par$y.calc.L2),
                                       lambda1 = lambda1,
                                       lambda2 = lambda2,
                                       ...)
        # note that y.calc is not correct for these parameter values.
        # this makes all results that depend on y.calc incorrect.
        # this is no problem as we will only use results that do not depend on 
        # y.calc
        
        dA <- ( res.par.j$Sigma.B.L1   %*% res.par.j$Sum.Sigma.L1.inv - A ) /
          par.inc[j] 
        a1 <- dA %*% ( y.obs.trans - res.par$y.calc.L1 )
        b1 <- A %*% ( - V.y.local[1:n1,j] )            
        V.yplusB[1:n1,j] <- V.y.local[1:n1,j] + a1 + b1
        V.B[1:n1,j]      <- a1 + b1
        V.Eps[,j]        <- - V.yplusB[1:n1,j]
        if ( n2 > 0 )
        {      
          dB <- ( res.par.j$Sigma.B.L2L1 %*% res.par.j$Sum.Sigma.L1.inv - B ) /
            par.inc[j]
          a2 <- dB %*% ( y.obs.trans - res.par$y.calc.L1 )
          b2 <- B %*% ( - V.y.local[1:n1,j] )             
          V.yplusB[n1+(1:n2),j] <- V.y.local[n1+(1:n2),j] + a2 + b2
          V.B[n1+(1:n2),j]      <- a2 + b2
        }
      }
    }
    else   # blockdesign
    {
      for( j in 1:length(par) )
      {
        par.j <- par
        par.j[j] <- par.j[j] + par.inc[j]
        if ( logpar ) par.j <- exp(par.j)
        res.par.j <- predict.bias.cond(par     = par.j,
                                       model   = model,
                                       L1      = L1.decoded,
                                       y.obs   = y.obs,
                                       L2      = L2,
                                       y.calc  = c(res.par$y.calc.L1,res.par$y.calc.L2),
                                       lambda1 = lambda1,
                                       lambda2 = lambda2,
                                       ...)
        # note that y.calc is not correct for these parameter values.
        # this makes all results that depend on y.calc incorrect.
        # this is no problem as we will only use results that do not depend on 
        # y.calc
        
        vars <- unique(L.decoded$var)
        for ( var in vars )
        {
          ind1 <- which(L1.decoded$var==var)
          A <- res.par$Sigma.B.L1[[var]]   %*% res.par$Sum.Sigma.L1.inv[[var]]
          dA <- ( res.par.j$Sigma.B.L1[[var]]   
                  %*% res.par.j$Sum.Sigma.L1.inv[[var]] - A ) /
            par.inc[j] 
          a1 <- dA %*% ( y.obs.trans[ind1] - res.par$y.calc.L1[ind1] )
          b1 <- A %*% ( - V.y.local[ind1,j] )
          V.yplusB[ind1,j] <- V.y.local[ind1,j] + a1 + b1
          V.B[ind1,j]      <- a1 + b1
          V.Eps[ind1,j]    <- - V.yplusB[ind1,j]
          ind2 <- which(L2.decoded$var==var)
          if ( length(ind2) > 0 )
          {
            B <- res.par$Sigma.B.L2L1[[var]] %*% res.par$Sum.Sigma.L1.inv[[var]]
            dB <- ( res.par.j$Sigma.B.L2L1[[var]] 
                    %*% res.par.j$Sum.Sigma.L1.inv[[var]] - B ) /
              par.inc[j]
            a2 <- dB %*% ( y.obs.trans[ind1] - res.par$y.calc.L1[ind1] )
            b2 <- B %*% ( - V.y.local[ind1,j] )
            V.yplusB[n1+ind2,j]  <- V.y.local[n1+ind2,j] + a2 + b2
            V.B[n1+ind2,j]       <- a2 + b2
          }
        }
      }
    }
  }
  
  # calculate results:
  
  y.cond      <- c(res.par$y.calc.L1,res.par$y.calc.L2)
  B.mean.cond <- c(res.par$B.mean.L1,res.par$B.mean.L2)
  Y.mean.cond <- c(res.par$Y.mean.L1,res.par$Y.mean.L2)
  B.var.cond  <- c(res.par$B.var.L1, res.par$B.var.L2)
  Y.var.cond  <- c(res.par$Y.var.L1, res.par$Y.var.L2)
  
  y.mean      <- y.cond
  B.mean      <- B.mean.cond
  B.var       <- diag( V.B %*% par.var.L1 %*% t(V.B) ) + B.var.cond
  Eps.mean    <- y.obs.trans - res.par$y.calc.L1 - res.par$B.mean.L1
  Eps.var     <- diag( V.Eps %*% par.var.L1 %*% t(V.Eps) ) + res.par$B.var.L1
  yplusB.mean <- y.cond + B.mean.cond
  yplusB.var  <- diag( V.yplusB %*% par.var.L1 %*% t(V.yplusB) ) + B.var.cond
  Y.mean      <- yplusB.mean
  Y.var       <- diag( V.yplusB %*% par.var.L1 %*% t(V.yplusB) ) + Y.var.cond
  
  names(y.mean)      <- L.encoded
  names(B.mean)      <- L.encoded
  names(B.var)       <- L.encoded
  names(Eps.mean)    <- L1.encoded
  names(Eps.var)     <- L1.encoded
  names(yplusB.mean) <- L.encoded
  names(yplusB.var)  <- L.encoded
  names(Y.mean)      <- L.encoded
  names(Y.var)       <- L.encoded
  
  return(list(y.mean      = y.mean,
              B.mean      = B.mean,
              B.var       = B.var,
              Eps.mean    = Eps.mean,
              Eps.var     = Eps.var,
              yplusB.mean = yplusB.mean,
              yplusB.var  = yplusB.var,
              Y.mean      = Y.mean,
              Y.var       = Y.var))
}


sysanal.intpol.irregular <- function(x,
                                     x.data,
                                     y.data,               
                                     method=c("neighbor","invdist"),
                                     scales.x=NA,
                                     fact=1)
{
  #   interpolation on an irregular grid
  #
  #   x:        vector of values of independent variable at which to interpolate;
  #             xout can also be a matrix, interpolation is then for each row
  #   x.data:   matrix of independent independent variables; each row of x.data
  #             belongs to one input specification (x.data can be a vector if there
  #             is only one input dimension)
  #   y.data:   vector of values of the dependent variable; one value for each
  #             row of x.data
  #   method:   "neighbor" does nearest neighbor interpolation (piecewise 
  #             constant function of value at nearest data point)
  #             "invdist2" weights the provided dependent values according to 
  #             the inverse of the distance to xout squared
  #   scales.x: values of x and of x.data are scaled by these scales before 
  #             interpolation; default is all scales equal to 1
  #
  #   --------------------------------------
  
  
  # call of each row if x is a matrix:
  
  if ( is.vector(x.data) )
  {
    x.data <- matrix(x.data,ncol=1)
    x      <- matrix(x,ncol=1)
  }
  if ( is.matrix(x) )
  {
    y <- rep(NA,nrow(x))
    for ( i in 1:nrow(x) )
    {
      y[i] <- sysanal.intpol.irregular(as.vector(x[i,]),x.data,y.data,method,scales.x,fact)
    }
    return(y)
  }
  
  # check input:
  
  n.data <- length(y.data)
  dim.x  <- length(x)
  y <- NA
  if ( nrow(x.data)!=n.data | ncol(x.data)!=dim.x )
  {
    print("incorrect dimensions in sysanal.interpol.irregular")
    return(y)
  }
  if ( is.na(scales.x[1]) ) scales.x <- rep(1,dim.x)
  
  # calculate distances to provided points:
  
  x.mat <- matrix(rep(x,n.data),nrow=n.data,byrow=TRUE)
  dist <- sqrt(apply(((x.mat-x.data)%*%diag(1/scales.x))^2,1,sum))
  
  ind.min <- which.min(dist)
  if ( method[1] == "neighbor" )
  {
    y <- y.data[ind.min]
  }
  else
  {
    if ( method[1] == "invdist" )
    {
      if ( dist[ind.min] == 0)
      {
        y <- y.data[ind.min]
      }
      else
      {
        w <- 1/dist^(dim.x*fact)
        w <- w/sum(w)
        y <- sum(w*y.data)
      }
    }
    else
    {
      print(paste("unknown method:",method[1]))
    }
  }
  
  return(y)
}



# didactical emulator (preliminary version, added May 30, 2012):
# ====================

emulate <- function(x, ...) UseMethod("emulate")

sysanal.emulator.create <- function(inp.design,
                                    res.design,
                                    par,
                                    sd,
                                    lambda,
                                    alpha    = 2,
                                    sd.smooth = 0,
                                    pri.mean = sysanal.emulator.pri.mean.lin,
                                    pri.var  = sysanal.emulator.pri.var.normal)
{
  #sysanal.package("fpc")
  
  n.design <- length(res.design)
  dim.inp  <- length(lambda)
  if ( is.vector(inp.design) ) 
  {
    inp.design <- matrix(inp.design,nrow=n.design)
  }
  if ( nrow(inp.design) != n.design ) stop("number of outputs must be equal to number of inputs")
  if ( ncol(inp.design) != dim.inp )  stop("dimension of input must be equal to length of lambda")
  
  emu <- list()
  emu$inp.design  <- inp.design
  emu$res.design  <- res.design
  emu$par         <- par
  emu$sd          <- sd
  emu$lambda      <- lambda
  emu$alpha       <- alpha
  emu$sd.smooth   <- sd.smooth
  emu$pri.mean    <- pri.mean
  emu$pri.var     <- pri.var
  emu$pri.var.num <- pri.var(inp=inp.design,sd=sd,lambda=lambda,alpha=alpha,sd.smooth=sd.smooth)
  emu$delta       <- res.design - pri.mean(inp=inp.design,par=par)
  emu$inv         <- solve(emu$pri.var.num)
  #emu$inv         <- solvecov(emu$pri.var.num)$inv
  emu$v           <- emu$inv %*% emu$delta
  class(emu) <- "emulator.GASP"
  
  return(emu)
}


sysanal.emulator.pri.mean.lin <- function(inp,par)
{
  if ( is.vector(inp) ) inp <- matrix(inp,ncol=length(par)-1)
  par <- as.numeric(par)
  
  n.inp   <- nrow(inp)
  dim.inp <- ncol(inp)
  dim.par <- length(par)
  
  if ( dim.par != dim.inp+1 ) stop("number of parameters must be dim.inp+1")
  
  res <- par[1] + inp %*% par[-1]
  rownames(res) <- 1:n.inp
  colnames(res) <- "y"
  
  return(res)
}


sysanal.emulator.pri.var.normal <- function(inp,sd,lambda,alpha=2,sd.smooth=0,rows=NA,cols=NA)
{
  if ( is.vector(inp) ) inp <- matrix(inp,ncol=length(lambda))
  n.inp   <- nrow(inp)
  dim.inp <- ncol(inp)
  if ( length(lambda) != dim.inp ) stop("length of lambda must be dim.inp")
  if ( is.na(rows[1]) ) rows <- 1:n.inp
  if ( is.na(cols[1]) ) cols <- 1:n.inp
  var <- matrix(0,nrow=length(rows),ncol=length(cols))
  
  s <- matrix(0,nrow=length(rows),ncol=length(cols))
  for ( k in 1:dim.inp )
  {
    s <- s + (abs(inp[rows,k]%o%rep(1,length(cols))-rep(1,length(rows))%o%inp[cols,k])/
                lambda[k])^alpha
  }
  var <- sd*sd*exp(-s)
  var <- var+diag(sd.smooth*sd.smooth,nrow=n.inp)[rows,cols]
  
  colnames(var) <- cols
  rownames(var) <- rows
  return(var)
}


emulate.emulator.GASP <- function(emulator,inp)
{
  if ( is.vector(inp) ) inp <- matrix(inp,ncol=length(emulator$lambda))
  n.inp <- nrow(inp)
  n.design <- nrow(emulator$inp.design)
  y   <- rep(NA,n.inp)
  var <- rep(NA,n.inp)
  for ( i in 1:n.inp )
  {
    pri.mean <- emulator$pri.mean(inp      = matrix(inp[i,],nrow=1),par=emulator$par)
    k        <- emulator$pri.var(inp       = rbind(emulator$inp.design,matrix(inp[i,],nrow=1)),
                                 sd        = emulator$sd,
                                 lambda    = emulator$lambda,
                                 alpha     = emulator$alpha,
                                 sd.smooth = emulator$sd.smooth,
                                 rows      = n.design+1,
                                 cols      = 1:n.design)
    K        <- emulator$pri.var(inp       = matrix(inp[i,],nrow=1),
                                 sd        = emulator$sd,
                                 lambda    = emulator$lambda,
                                 alpha     = emulator$alpha,
                                 sd.smooth = emulator$sd.smooth)
    y[i]   <- pri.mean + k %*% emulator$v
    var[i] <- K - k %*% emulator$inv %*% t(k)
  }
  return(list(inp=inp,y=y,var=var))
}




