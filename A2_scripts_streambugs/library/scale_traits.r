################################################################
#
# Convert trait data in affinity scores
# =====================================
#
# creation:      06.12.2013
#
################################################################

# This function adds affinity scores of the "indifferent" category to the others and normalizes the 
# affinity scores between 0 and 1.

max.trait.score <- function(x)
{
  if(sum(is.na(x)) == length(x)){
    return(NA)
  } else {
    return(max(x, na.rm=TRUE))
  }
}
  
scale.traits <- function(db.traits)
{
  db.traits.temp <- as.data.frame(matrix(NA, nrow=nrow(db.traits), ncol=ncol(db.traits), dimnames=dimnames(db.traits)))
  db.traits.temp[,"Taxon"] <- as.character(db.traits[,"Taxon"])
  
  ind.rm <- c()
  
    ## Divide the trait scores by the maximum score attributed
    ## to a given taxon for a given trait
    ## --------------------------------------------------------
      
      # 'Indiferrent' categories
      cind <- c(grep("_ind", colnames(db.traits)),     
               grep("_eut", colnames(db.traits)),     # we have to document this for including other databases
               grep("_eury", colnames(db.traits)))
      cind2 <- grep("Taxon", colnames(db.traits))
  
  if(sum(cind) > 0 & length(cind) > 1)
  {
    warning("More than one trait category correspond to an indeterminate category")
  } else {
    if(sum(cind) > 0 & length(cind) == 1)
    {

      db.traits.temp[is.na(db.traits[,cind]),- cind2] =                                # just copy rows with NA in the "indifferent" category:
        db.traits[is.na(db.traits[,cind]),- cind2]
  
      db.traits.temp[!is.na(db.traits[,cind]),-c(cind, cind2)] <-                      # add affinity scores of indifferent to other categories
        replace(x      = db.traits[!is.na(db.traits[,cind]),-c(cind, cind2)],          # replace NA with 0
                list   = is.na(db.traits[!is.na(db.traits[,cind]), -c(cind, cind2)]),
                values = 0) + 
        replace(x      = db.traits[!is.na(db.traits[,cind]),cind],                     # replace NA with 0 
                list   = is.na(db.traits[!is.na(db.traits[,cind]),cind]),              
                values = 0)

    } else {
      db.traits.temp[,-cind2] <- db.traits[,-cind2]
    }
    
    if(sum(cind) > 0 & length(cind) == 1)                                             # remove indifferent category
    {
      db.traits.temp <- db.traits.temp[, - cind] 
    }
    
    maxi <- apply(db.traits.temp[,- cind2], 1, max.trait.score)                       # determine max. trait score                                
    
    db.traits.temp[,- cind2] <- db.traits.temp[,- cind2] / maxi                       # normalize scores between 0 and 1
    
  }  
  return(db.traits.temp)
}
