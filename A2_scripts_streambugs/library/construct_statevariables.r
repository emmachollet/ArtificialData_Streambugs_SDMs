
construct.statevariables <- function(POM,Algae,Invertebrates,Reaches,Habitats)
{
  # construct state variables for all reaches
  if(length(Reaches)!=length(Habitats)) stop("No of Reaches",length(Reaches),"not equal to",
                                            "No of Habitats", length(Habitats))
  ynames <- character(0)
  for ( i in 1:length(Reaches) )
  {
    ynames <- c(ynames,
                paste(Reaches[i],Habitats[i],
                      c(paste(POM,"POM",sep="_"),
                        paste(Algae,"Algae",sep="_"),
                        paste(Invertebrates,"Invertebrates",sep="_")),
                      sep="_"))
  }
  if ( length(ynames) != length(Reaches)*(length(POM)+length(Algae)+length(Invertebrates)) ) stop("problem constructing state variables")
  
  return(ynames)
}

# construct.ini.pars <- function(y.names,y.ini,method="fixed")
# 
# # set initial values of state variables # units?   
#   
# {
#   y.ini <- rep(y.ini,length(ynames))
#       
#   if(method=="runif")
#   {
#      y.ini <- runif(length(ynames),min=0.5*y.ini, max=1.5*y.ini)
#   }
#   
#   names(y.ini) <- ynames
#   
#   return(y.ini)
# }