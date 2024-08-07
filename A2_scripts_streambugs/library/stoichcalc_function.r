if ( !require(stoichcalc) )
{
   install.packages("stoichcalc")
   library(stoichcalc)
}


calc.stoich <- function( par, returns="parout",verbose=FALSE ) #returns=nu
{
    # par = par.stoich.default
    
    # define composition vectors of substances an organisms:
    # ------------------------------------------------------

    NH4    <- c(H      = 4*1/14,              # gH/gNH4-N
                N      = 1,                   # gN/gNH4-N
                charge = 1/14)                # chargeunits/gNH4-N
  # NO3    <- c(O      = 3*16/14,             # gO/gNO3-N
  #             N      = 1,                   # gN/gNO3-N
  #             charge = 1/14)                # chargeunits/gNO3-N
    HPO4   <- c(O      = 4*16/31,             # gO/gHPO4-P
                H      = 1*1/31,              # gH/gHPO4-P
                P      = 1,                   # gP/gHPO4-P
                charge = -2/31)               # chargeunits/gHPO4-P
    HCO3   <- c(C      = 1,                   # gC/gHCO3-C
                O      = 3*16/12,             # gO/gHCO3-C
                H      = 1*1/12,              # gH/gHCO3-C
                charge = -1/12)               # chargeunits/gHCO3-C
    O2     <- c(O      = 1)                   # gO/gO2-O
    H      <- c(H      = 1,                   # gH/molH
                charge = 1)                   # chargeunits/molH
    H2O    <- c(O      = 1*12,                # gO/molH2O
                H      = 2*1)                 # gH/molH2O
    Algae    <- c(C    = par$Algae_aC,        # gC/gAlgae
                  O      = par$Algae_aO,      # gO/gAlgae
                  H      = par$Algae_aH,      # gH/gAlgae
                  N      = par$Algae_aN,      # gN/gAlgae
                  P      = par$Algae_aP       # gP/gAlgae
                 )
    filamentousAlgae <- crustyAlgae <- Algae
    Invertebrates   <- c(C      = par$Invertebrates_aC,      # gC/gInvertebrates
                         O      = par$Invertebrates_aO,      # gO/gInvertebrates
                         H      = par$Invertebrates_aH,      # gH/gInvertebrates
                         N      = par$Invertebrates_aN,      # gN/gInvertebrates
                         P      = par$Invertebrates_aP       # gP/gInvertebrates
                        )
    FPOM   <- c(C      = par$FPOM_aC,     # gC/gPOM
                O      = par$FPOM_aO,     # gO/gPOM
                H      = par$FPOM_aH,     # gH/gPOM
                N      = par$FPOM_aN,     # gN/gPOM
                P      = par$FPOM_aP      # gP/gPOM
               )
   
    if(length(par$CPOM_aC) > 0)
    {
      CPOM <- CPOMp  <- c(C      = par$CPOM_aC,     # gC/gPOM
                          O      = par$CPOM_aO,     # gO/gPOM
                          H      = par$CPOM_aH,     # gH/gPOM
                          N      = par$CPOM_aN,     # gN/gPOM
                          P      = par$CPOM_aP      # gP/gPOM
      )
      
    } else {
      CPOM <- CPOMp   <- FPOM
    }

    CPOMa   <- Invertebrates
    SusPOM  <- FPOM
    PREY    <- Invertebrates
    
    # define list of composition vectors:
    # -----------------------------------

    subst.comp <- list(NH4  = NH4,
                     # NO3  = NO3,
                       HPO4 = HPO4,
                       HCO3 = HCO3,
                       O2   = O2,
                       H    = H,
                       H2O  = H2O,
                       Algae = Algae,
                       filamentousAlgae  = filamentousAlgae,
                       crustyAlgae = crustyAlgae,
                       Invertebrates = Invertebrates,
                       FPOM = FPOM,
                       SusPOM = SusPOM,
                       PREY = PREY,
                       CPOMp = CPOMp,
                       CPOMa = CPOMa,
                       CPOM = CPOM)

    # compile and print composition matrix:
    # -------------------------------------

    alpha <- calc.comp.matrix(subst.comp,verbose=verbose)

   # print(alpha)
   # nu.basis <- calc.stoich.basis(alpha)
   # print(nu.basis)

    #test
   #  nu.basis %*% t(alpha)

    # Derivation of Process Stoichiometry:
    # ====================================

    # growth of scraper (general case: filamentous 
    # and crusty algae with the same elemental composition)
    # ------------------------------------------------------

    # substances/organisms relevant for growth of consumers:

    subst.gro.CONS.scra <- c("NH4","HPO4","HCO3","O2","H","H2O",#"DOM",
                             "Algae","Invertebrates","FPOM")

    # define parameters for constraints:

    par$Y.CONS.scra <- min(
       1,
       (par$Algae_EC-par$Scra_fe*par$POM_EC)/par$Invertebrates_EC,
       (alpha["N","Algae"]-par$Scra_fe*alpha["N","FPOM"])/alpha["N","Invertebrates"],
       (alpha["P","Algae"]-par$Scra_fe*alpha["P","FPOM"])/alpha["P","Invertebrates"])
             

    # process definition

    const.gro.CONS.scra <- list(c("Invertebrates" = 1,"Algae" = par$Y.CONS.scra),
                                c("FPOM"  = 1,"Algae" = par$Scra_fe) )

    nu.gro.CONS.scra    <- calc.stoich.coef (alpha       = alpha,
                                             name        = "gro.scra",
                                             subst       = subst.gro.CONS.scra,
                                             subst.norm  = "Invertebrates",
                                             nu.norm     = 1,
                                             constraints = const.gro.CONS.scra,
                                             verbose=verbose)
  #  print(nu.gro.CONS.scra)

    
    # growth of shredder (general case)
    # ---------------------------------
    
    # substances/organisms relevant for growth of consumers:
    
    subst.gro.CONS.shred <- c("NH4","HPO4","HCO3","O2","H","H2O",
                              "Invertebrates","CPOM","FPOM")    #,"DOM"
    
    # define parameters for constraints:
    
    par$Y.CONS.shred <- min(1,
                            (par$POM_EC-par$Shred_fe*par$POM_EC)/par$Invertebrates_EC,
                            (alpha["N","CPOM"]-par$Shred_fe*alpha["N","FPOM"])/alpha["N","Invertebrates"],
                            (alpha["P","CPOM"]-par$Shred_fe*alpha["P","FPOM"])/alpha["P","Invertebrates"])
    
    
    # process definition
    
    const.gro.CONS.shred <- list(c("Invertebrates" = 1,"CPOM" = par$Y.CONS.shred),
                                 c("FPOM" = 1,"CPOM" = par$Shred_fe) )
    
    nu.gro.CONS.shred  <- calc.stoich.coef (alpha       = alpha,
                                            name        = "gro.shred",
                                            subst       = subst.gro.CONS.shred,
                                            subst.norm  = "Invertebrates",
                                            nu.norm     = 1,
                                            constraints = const.gro.CONS.shred,
                                            verbose=verbose)
    #  print(nu.gro.CONS.shred)
    
    
    # growth of shredder feeding on dead plants
    # ------------------------------------------

    # substances/organisms relevant for growth of consumers:

    subst.gro.CONS.shredp <- c("NH4","HPO4","HCO3","O2","H","H2O",
                              "Invertebrates","CPOMp","FPOM")    #,"DOM"

    # define parameters for constraints:

    par$Y.CONS.shredp <- min(1,
       (par$POM_EC-par$Shred_fe*par$POM_EC)/par$Invertebrates_EC,
       (alpha["N","CPOMp"]-par$Shred_fe*alpha["N","FPOM"])/alpha["N","Invertebrates"],
       (alpha["P","CPOMp"]-par$Shred_fe*alpha["P","FPOM"])/alpha["P","Invertebrates"])
  
       
    # process definition

    const.gro.CONS.shredp <- list(c("Invertebrates" = 1,"CPOMp" = par$Y.CONS.shredp),
                                 c("FPOM" = 1,"CPOMp" = par$Shred_fe) )

    nu.gro.CONS.shredp  <- calc.stoich.coef (alpha       = alpha,
                                            name        = "gro.shredp",
                                            subst       = subst.gro.CONS.shredp,
                                            subst.norm  = "Invertebrates",
                                            nu.norm     = 1,
                                            constraints = const.gro.CONS.shredp,
                                            verbose=verbose)
  #  print(nu.gro.CONS.shredp)

    # growth of shredder feeding on dead animals
    # ------------------------------------------
    
    # substances/organisms relevant for growth of consumers:
    
    subst.gro.CONS.shreda <- c("NH4","HPO4","HCO3","O2","H","H2O",
                               "Invertebrates","CPOMa","FPOM")    #,"DOM"
    
    # define parameters for constraints:
    
    par$Y.CONS.shreda <- min(1,
                             (par$Invertebrates_EC-par$Shred_fe*par$POM_EC)/par$Invertebrates_EC,
                             (alpha["N","CPOMa"]-par$Shred_fe*alpha["N","FPOM"])/alpha["N","Invertebrates"],
                             (alpha["P","CPOMa"]-par$Shred_fe*alpha["P","FPOM"])/alpha["P","Invertebrates"])
    
    
    # process definition
    
    const.gro.CONS.shreda <- list(c("Invertebrates" = 1,"CPOMa" = par$Y.CONS.shreda),
                                  c("FPOM" = 1,"CPOMa" = par$Shred_fe) )
    
    nu.gro.CONS.shreda  <- calc.stoich.coef (alpha       = alpha,
                                             name        = "gro.shreda",
                                             subst       = subst.gro.CONS.shreda,
                                             subst.norm  = "Invertebrates",
                                             nu.norm     = 1,
                                             constraints = const.gro.CONS.shreda,
                                             verbose=verbose)
    #  print(nu.gro.CONS.shreda)
    
    # growth of shredder feeding on living macrophytes
    # -------------------------------------------------
    
    # substances/organisms relevant for growth of consumers:
    
    subst.gro.CONS.shredm <- c("NH4","HPO4","HCO3","O2","H","H2O",
                               "Invertebrates","filamentousAlgae","FPOM")    #,"DOM"
    
    # define parameters for constraints:
    
    par$Y.CONS.shredm <- min(1,
                             (par$filamentousAlgae_EC-par$Shred_fe*par$POM_EC)/par$Invertebrates_EC,
                             (alpha["N","filamentousAlgae"]-par$Shred_fe*alpha["N","FPOM"])/alpha["N","Invertebrates"],
                             (alpha["P","filamentousAlgae"]-par$Shred_fe*alpha["P","FPOM"])/alpha["P","Invertebrates"])
    
    
    # process definition
    
    const.gro.CONS.shredm <- list(c("Invertebrates" = 1,"filamentousAlgae" = par$Y.CONS.shredm),
                                  c("FPOM" = 1,"filamentousAlgae" = par$Shred_fe) )
    
    nu.gro.CONS.shredm  <- calc.stoich.coef (alpha       = alpha,
                                             name        = "gro.shredm",
                                             subst       = subst.gro.CONS.shredm,
                                             subst.norm  = "Invertebrates",
                                             nu.norm     = 1,
                                             constraints = const.gro.CONS.shredm,
                                             verbose=verbose)
    #  print(nu.gro.CONS.shredm)
    
    # growth of filterer feeding on FPOM
    # ----------------------------------

    # substances/organisms relevant for growth of consumers:

    subst.gro.CONS.filt <- c("NH4","HPO4","HCO3","O2","H","H2O",
                              "Invertebrates","SusPOM","FPOM")    #,"DOM"

    # define parameters for constraints:

    par$Y.CONS.filt <- min(1,
       (par$POM_EC-par$Filt_fe*par$POM_EC)/par$Invertebrates_EC,
       (alpha["N","SusPOM"]-par$Filt_fe*alpha["N","FPOM"])/alpha["N","Invertebrates"],
       (alpha["P","SusPOM"]-par$Filt_fe*alpha["P","FPOM"])/alpha["P","Invertebrates"])


    # process definition

    const.gro.CONS.filt <- list(c("Invertebrates" = 1,"SusPOM" = par$Y.CONS.filt),
                                c("FPOM" = 1,"SusPOM" = par$Filt_fe) )

    nu.gro.CONS.filt  <- calc.stoich.coef (alpha       = alpha,
                                           name        = "gro.filt",
                                           subst       = subst.gro.CONS.filt,
                                           subst.norm  = "Invertebrates",
                                           nu.norm     = 1,
                                           constraints = const.gro.CONS.filt,
                                           verbose=verbose)
  #  print(nu.gro.CONS.filt)

    # growth of collectors
    # --------------------

    # substances/organisms relevant for growth of consumers:

    subst.gro.CONS.coll <- c("NH4","HPO4","HCO3","O2","H","H2O","Invertebrates","FPOM")

    # define parameters for constraints:
    par$Y.CONS.coll <- min(1,
       (par$POM_EC-par$Coll_fe*par$POM_EC)/par$Invertebrates_EC,
       (alpha["N","FPOM"]-par$Coll_fe*alpha["N","FPOM"])/alpha["N","Invertebrates"],
       (alpha["P","FPOM"]-par$Coll_fe*alpha["P","FPOM"])/alpha["P","Invertebrates"])


    # process definition

    const.gro.CONS.coll <- list(c("Invertebrates" = -(par$Coll_fe-1),
                                  "FPOM" = par$Y.CONS.coll))

    nu.gro.CONS.coll    <- calc.stoich.coef (alpha       = alpha,
                                             name        = "gro.coll",
                                             subst       = subst.gro.CONS.coll,
                                             subst.norm  = "Invertebrates",
                                             nu.norm     = 1,
                                             constraints = const.gro.CONS.coll,
                                             verbose=verbose)
 #   print(nu.gro.CONS.coll)

    # growth of predators
    # --------------------

    # substances/organisms relevant for growth of consumers:

    subst.gro.CONS.pred <- c("NH4","HPO4","HCO3","O2","H","H2O","Invertebrates","FPOM","PREY")

    # define parameters for constraints:

    par$Y.CONS.pred <- min(1,
          (par$Invertebrates_EC-par$Pred_fe*par$POM_EC)/par$Invertebrates_EC,                            
          (alpha["N","Invertebrates"]-par$Pred_fe*alpha["N","FPOM"])/alpha["N","Invertebrates"],
          (alpha["P","Invertebrates"]-par$Pred_fe*alpha["P","FPOM"])/alpha["P","Invertebrates"])                       

    # process definition

    const.gro.CONS.pred <- list(c("Invertebrates" = 1,"PREY" = par$Y.CONS.pred),
                                c("FPOM" = 1,"PREY" = par$Pred_fe))


    nu.gro.CONS.pred    <- calc.stoich.coef (alpha       = alpha,
                                        name        = "gro.pred",
                                        subst       = subst.gro.CONS.pred,
                                        subst.norm  = "Invertebrates",
                                        nu.norm     = 1,
                                        constraints = const.gro.CONS.pred,
                                        verbose=verbose)
 #   print(nu.gro.CONS.pred)

    # growth of piercers feeding on living macrophytes
    # ------------------------------------------------
    
    # substances/organisms relevant for growth of consumers:
    
    subst.gro.CONS.piercm <- c("NH4","HPO4","HCO3","O2","H","H2O",
                               "Invertebrates","filamentousAlgae","FPOM")    #,"DOM"
    
    # define parameters for constraints:
    
    par$Y.CONS.piercm <- min(1,
                             (par$Algae_EC-par$Pierc_fe*par$POM_EC)/par$Invertebrates_EC,
                             (alpha["N","filamentousAlgae"]-par$Pierc_fe*alpha["N","FPOM"])/alpha["N","Invertebrates"],
                             (alpha["P","filamentousAlgae"]-par$Pierc_fe*alpha["P","FPOM"])/alpha["P","Invertebrates"])
    
    
    # process definition
    
    const.gro.CONS.piercm <- list(c("Invertebrates" = 1,"filamentousAlgae" = par$Y.CONS.piercm),
                                  c("FPOM" = 1,"filamentousAlgae" = par$Pierc_fe) )
    
    nu.gro.CONS.piercm  <- calc.stoich.coef (alpha       = alpha,
                                             name        = "gro.piercm",
                                             subst       = subst.gro.CONS.piercm,
                                             subst.norm  = "Invertebrates",
                                             nu.norm     = 1,
                                             constraints = const.gro.CONS.piercm,
                                             verbose=verbose)
    #  print(nu.gro.CONS.piercm)
    
    # growth of piercers feeding on animals
    # -------------------------------------
    
    # substances/organisms relevant for growth of consumers:
    
    subst.gro.CONS.pierca <- c("NH4","HPO4","HCO3","O2","H","H2O",
                               "Invertebrates", "PREY", "FPOM")    #,"DOM"
    
    # define parameters for constraints:
    
    par$Y.CONS.pierca <- min(1,
                             (par$Invertebrates_EC-par$Pierc_fe*par$POM_EC)/par$Invertebrates_EC,
                             (alpha["N","PREY"]-par$Pierc_fe*alpha["N","FPOM"])/alpha["N","Invertebrates"],
                             (alpha["P","PREY"]-par$Pierc_fe*alpha["P","FPOM"])/alpha["P","Invertebrates"])
    
    
    # process definition
    
    const.gro.CONS.pierca <- list(c("Invertebrates" = 1,"PREY" = par$Y.CONS.pierca),
                                  c("FPOM" = 1,"PREY" = par$Pierc_fe) )
    
    nu.gro.CONS.pierca  <- calc.stoich.coef (alpha       = alpha,
                                             name        = "gro.pierca",
                                             subst       = subst.gro.CONS.pierca,
                                             subst.norm  = "Invertebrates",
                                             nu.norm     = 1,
                                             constraints = const.gro.CONS.pierca,
                                             verbose=verbose)
    #  print(nu.gro.CONS.pierca)
    
    # death of algae (crusty and filamentous)
    # ---------------------------------------

      # substances/organisms relevant for death of periphyton:
      
      subst.death.alg <- c("NH4","HPO4","HCO3","O2","H","H2O","Algae","FPOM")

    # define parameters for constraints:

     par$Y.death.alg <- min(1,
                              par$Algae_EC/par$POM_EC,
                              par$Algae_aN /par$FPOM_aN,
                              par$Algae_aP /par$FPOM_aP)

    # process definition

       const.death.alg <- list(c("FPOM" = 1,"Algae" = par$Y.death.alg))
      
       nu.death.alg    <- calc.stoich.coef (alpha       = alpha,
                                            name        = "death.alg",
                                            subst       = subst.death.alg,
                                            subst.norm  = "Algae",
                                            nu.norm     = -1,
                                            constraints = const.death.alg,
                                            verbose=verbose)
       #print(nu.death.alg)
       if ( round(nu.death.alg[,"O2"],digits=10) < 0 ) cat("\n !!! warning: death.alg consumes oxygen! nu O2: ",
                                             round(nu.death.alg[,"O2"],digits=3),"\n\n  ")    
      
    
    # death of consumers
    # ------------------    
         
    ## Direct transformation of dead invertebrates in FPOM
    ## ---------------------------------------------------
    
    # substances/organisms relevant for death of consumers:
    
    subst.death1.cons <- c("NH4","HPO4","HCO3","O2","H","H2O","Invertebrates","FPOM")
    
    # define parameters for constraints:
    
    par$Y.death1.cons <- min(1,
                            par$Invertebrates_EC/par$POM_EC,
                            par$Invertebrates_aN /par$FPOM_aN,
                            par$Invertebrates_aP /par$FPOM_aP)
    
    # process definition
    
    const.death1.cons <- list(c("FPOM" = 1,"Invertebrates" = par$Y.death1.cons))
    
    nu.death1.cons    <- calc.stoich.coef (alpha       = alpha,
                                          name        = "death1.cons",
                                          subst       = subst.death1.cons,
                                          subst.norm  = "Invertebrates",
                                          nu.norm     = -1,
                                          constraints = const.death1.cons,
                                          verbose=verbose)
    #Achtung - Sauerstoffverbrauch!
    if ( round(nu.death1.cons[,"O2"],digits=10) < 0 ) cat("\n !!! warning: death1.cons consumes oxygen! nu O2: ",
                                                         round(nu.death1.cons[,"O2"],digits=3),"\n\n  ")  
    #print(nu.death1.cons)
    
    ## Transformation of dead invertebrates in CPOMa
    ## ---------------------------------------------
    
    # substances/organisms relevant for death of consumers:
    
    subst.death2.cons <- c("NH4","HPO4","HCO3","O2","H","H2O","Invertebrates","CPOMa")
    
    # define parameters for constraints:
    
    par$Y.death2.cons <- min(1,
                             par$Invertebrates_EC/par$Invertebrates_EC,
                             par$Invertebrates_aN /par$Invertebrates_aN,
                             par$Invertebrates_aP /par$Invertebrates_aP)
    
    # process definition
    
    const.death2.cons <- list(c("CPOMa" = 1,"Invertebrates" = par$Y.death2.cons))
    
    nu.death2.cons    <- calc.stoich.coef (alpha       = alpha,
                                           name        = "death2.cons",
                                           subst       = subst.death2.cons,
                                           subst.norm  = "Invertebrates",
                                           nu.norm     = -1,
                                           constraints = const.death2.cons,
                                           verbose=verbose)
    #Achtung - Sauerstoffverbrauch!
    if ( round(nu.death2.cons[,"O2"],digits=10) < 0 ) cat("\n !!! warning: death2.cons consumes oxygen! nu O2: ",
                                                          round(nu.death2.cons[,"O2"],digits=3),"\n\n  ")  
    #print(nu.death2.cons)
    
       
       
    # combine stoichiometric matrix nu:
    # =================================

    nu <- rbind(nu.gro.CONS.scra,
                nu.gro.CONS.shred,
                nu.gro.CONS.shredp,
                nu.gro.CONS.shreda,
                nu.gro.CONS.shredm,
                nu.gro.CONS.filt,
                nu.gro.CONS.coll,
                nu.gro.CONS.pred,
                nu.gro.CONS.piercm,
                nu.gro.CONS.pierca,
                nu.death.alg,
                nu.death1.cons,
                nu.death2.cons)
    
 
    if(returns=="parout")
      
    {  
    
     # combine parameter vector as output: 
     # ===================================

     parout <- c(
     
                 "Cons_Invertebrates_Scra_Scra"  = nu["gro.scra","Invertebrates"],
                 "Cons_Invertebrates_Scra_Algae" = nu["gro.scra","Algae"],
                 "Cons_Invertebrates_Scra_FPOM"  = nu["gro.scra","FPOM"],
                       
                 "Cons_Invertebrates_Shred_Shred" = nu["gro.shred","Invertebrates"],
                 "Cons_Invertebrates_Shred_CPOM"  = nu["gro.shred","CPOM"],
                 "Cons_Invertebrates_Shred_FPOM"  = nu["gro.shred","FPOM"],
                 
                 "Cons_Invertebrates_Shredp_Shredp" = nu["gro.shredp","Invertebrates"],
                 "Cons_Invertebrates_Shredp_CPOMp"  = nu["gro.shredp","CPOMp"],
                 "Cons_Invertebrates_Shredp_FPOM"   = nu["gro.shredp","FPOM"],
                 
                 "Cons_Invertebrates_Shredm_Shredm" = nu["gro.shredm","Invertebrates"],
                 "Cons_Invertebrates_Shredm_filamentousAlgae"  = nu["gro.shredm","filamentousAlgae"],
                 "Cons_Invertebrates_Shredm_FPOM"   = nu["gro.shredm","FPOM"],
                 
                 "Cons_Invertebrates_Shreda_Shreda" = nu["gro.shreda","Invertebrates"],
                 "Cons_Invertebrates_Shreda_CPOMa"  = nu["gro.shreda","CPOMa"],
                 "Cons_Invertebrates_Shreda_FPOM"   = nu["gro.shreda","FPOM"],
                 
                 "Cons_Invertebrates_Filt_Filt"   = nu["gro.filt","Invertebrates"],
                 "Cons_Invertebrates_Filt_SusPOM" = nu["gro.filt","SusPOM"],
                 "Cons_Invertebrates_Filt_FPOM"   = nu["gro.filt","FPOM"],

                 "Cons_Invertebrates_Coll_Coll"   = nu["gro.coll","Invertebrates"],
                 "Cons_Invertebrates_Coll_FPOM"   = nu["gro.coll","FPOM"],
 
                 "Cons_Invertebrates_Pred_Pred"   = nu["gro.pred","Invertebrates"],
                 "Cons_Invertebrates_Pred_Prey"   = nu["gro.pred","PREY"],
                 "Cons_Invertebrates_Pred_FPOM"   = nu["gro.pred","FPOM"],
     
                 "Cons_Invertebrates_Piercm_Piercm" = nu["gro.piercm","Invertebrates"],
                 "Cons_Invertebrates_Piercm_filamentousAlgae"  = nu["gro.piercm","filamentousAlgae"],
                 "Cons_Invertebrates_Piercm_FPOM"   = nu["gro.piercm","FPOM"],
                 
                 "Cons_Invertebrates_Pierca_Pierca" = nu["gro.pierca","Invertebrates"],
                 "Cons_Invertebrates_Pierca_Prey"   = nu["gro.pierca","PREY"],
                 "Cons_Invertebrates_Pierca_FPOM"   = nu["gro.pierca","FPOM"],
                 
                 "Death_Algae_Algae" = nu["death.alg","Algae"],
                 "Death_Algae_FPOM"  = nu["death.alg","FPOM"],  
                 "Death_Algae_O2"    = nu["death.alg","O2"],
                                 
                 "Death1_Invertebrates_Invertebrates" = nu["death1.cons","Invertebrates"],     
                 "Death1_Invertebrates_FPOM"         = nu["death1.cons","FPOM"],
                 "Death1_Invertebrates_O2"            = nu["death1.cons","O2"],
                 
                 "Death2_Invertebrates_Invertebrates" = nu["death2.cons","Invertebrates"],     
                 "Death2_Invertebrates_CPOMa"         = nu["death2.cons","CPOMa"],
                 "Death2_Invertebrates_O2"            = nu["death2.cons","O2"]
                 
                 )
    
     return(parout)
    } else
    
    {
     if(returns=="nu")
     {
       return(nu)
     } else warning("returns argument ",returns," unknown\n")
    }  
       
}

par.stoich.default <- list(  Invertebrates_aC = 0.36,      # gC/gInvertebrates
                             Invertebrates_aO = 0.50,      # gO/gInvertebrates
                             Invertebrates_aH = 0.07,      # gH/gInvertebrates
                             Invertebrates_aN = 0.06,      # gN/gInvertebrates
                             Invertebrates_aP = 0.01,      # gP/gInvertebrates
                             Algae_aC         = 0.36,      # gC/gAlgae
                             Algae_aO         = 0.50,      # gO/gAlgae
                             Algae_aH         = 0.07,      # gH/gAlgae
                             Algae_aN         = 0.06,      # gN/gAlgae
                             Algae_aP         = 0.01,      # gP/gAlgae
                             FPOM_aC          = 0.36,      # gC/gPOM
                             FPOM_aO          = 0.50,      # gO/gPOM
                             FPOM_aH          = 0.07,      # gH/gPOM
                             FPOM_aN          = 0.06,      # gN/gPOM
                             FPOM_aP          = 0.01,      # gP/gPOM
                             Invertebrates_EC = 22000,     # J/gInvertebrates
                             Algae_EC         = 20000,     # J/gAlgae
                             POM_EC           = 18000,     # J/gPOM
                             Scra_fe          = 0.1,       # gDM/gDM
                             Coll_fe          = 0.1,       # gDM/gDM
                             Filt_fe          = 0.1,       # gDM/gDM
                             Shred_fe         = 0.6,       # gDM/gDM
                             Pred_fe          = 0.2,       # gDM/gDM 
                             Pierc_fe         = 0.1        # gDM/gDM
                          )
                            
                            
# par.stoich.default$POM_EC      = 24000
# parout <- calc.stoich( par=par.stoich.default,returns="parout")
# nu <- calc.stoich( par=par.stoich.default,returns="nu")
# testwarn <- calc.stoich( par=par.stoich.default,returns="Schmarrn" )