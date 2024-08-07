derive.par.env.traits <- function(file.db.temp=NA,
                                  file.db.current=NA,
                                  file.db.sapro=NA,
                                  file.db.pH=NA,
                                  file.db.spear=NA)
{
  par.env.trait.global <- c()
  
  if(!is.na(file.db.temp))
  {
    if(grepl("freshwaterecology", file.db.temp))
    {
      par.env.trait.global <- c(par.env.trait.global,
                                "tempmaxKval_class1"=3+273.15,  # K
                                "tempmaxKval_class2"=8+273.15,
                                "tempmaxKval_class3"=14+273.15,
                                "tempmaxKval_class4"=22+273.15)
    }
    
    if(grepl("tachet", file.db.temp))
    {
      par.env.trait.global <- c(par.env.trait.global,
                                "tempmaxKval_class1"=7.5+273.15,  # K
                                "tempmaxKval_class2"=22.5+273.15)
    }
  }
  
  if(!is.na(file.db.current))
  {
    if(grepl("freshwaterecology", file.db.current))
    {
      warning("Current velocity trait from freshwaterecology data base is not implemented yet")
    }
    
    if(grepl("tachet", file.db.current))
    {
      par.env.trait.global <- c(par.env.trait.global,
                                "currentmsval_class1" = 0,             # m/s
                                "currentmsval_class2" = 0.125,
                                "currentmsval_class3" = 0.375,
                                "currentmsval_class4" = 0.625)
    }
  }
  
  if(!is.na(file.db.sapro))
  {
    if(grepl("freshwaterecology", file.db.sapro))
    {
      par.env.trait.global <- c(par.env.trait.global,
                                "saprowqclassval_class0" = 0,        # water quality class 
                                "saprowqclassval_class1" = 1,        # water quality class 
                                "saprowqclassval_class2" = 2, 
                                "saprowqclassval_class3" = 3, 
                                "saprowqclassval_class4" = 4)
    }
    
    if(grepl("tachet", file.db.sapro))
    {
      par.env.trait.global <- c(par.env.trait.global,
                                "saprowqclassval_class0" = 0,        # water quality class 
                                "saprowqclassval_class1" = 1,        # water quality class 
                                "saprowqclassval_class2" = 2, 
                                "saprowqclassval_class3" = 3, 
                                "saprowqclassval_class4" = 4)
    }
  }
  
  if(!is.na(file.db.pH))
  {
    if(grepl("freshwaterecology", file.db.pH))
    {
      par.env.trait.global <- c(par.env.trait.global,
                                "pHval_class1" = 3.5,        # pH
                                "pHval_class2" = 10.5)
    }
    
    if(grepl("tachet", file.db.pH))
    {
      par.env.trait.global <- c(par.env.trait.global,
                                "pHval_class1" = 2,        # pH 
                                "pHval_class2" = 4.25, 
                                "pHval_class3" = 4.75, 
                                "pHval_class4" = 5.25,
                                "pHval_class5" = 5.75,
                                "pHval_class6" = 6.25)
    }
  }
  
  if(!is.na(file.db.spear))
  {
    par.env.trait.global <- c(par.env.trait.global,
                              "orgmicropollTUval_class1"  = -4,       # Toxic unit
                              "orgmicropollTUval_class2"  =  0)
  }
  
  return(par.env.trait.global)
  
}