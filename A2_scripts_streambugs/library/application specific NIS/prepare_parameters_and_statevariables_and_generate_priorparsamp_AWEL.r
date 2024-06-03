#source("input/model_inputs_catchment_AWEL_ABC.r")

#observed.abund <- read.table("input/observed_taxa_AWEL.dat",sep="\t",header=TRUE)


# set up state variables: 
# -----------------------

Invertebrates     <- c(t(Invertebrates))
n.Invertebrates   <- length(Invertebrates)


CPOM.dif <- FALSE

if(CPOM.dif == TRUE)
{
  POM               <- c("FPOM", "CPOMa", "CPOMp")
}
if(CPOM.dif == FALSE)
{
  POM               <- c("FPOM","CPOM")
}
n.POM             <- length(POM)

n.Algae           <- 2
Algae             <- paste("Alga",1:n.Algae,sep="")

Algae <- c("crustyAlgae","filamentousAlgae")

y.names <- construct.statevariables(POM,Algae,Invertebrates,Reaches,Habitats)

y.names <- decode.statevarnames(y.names)


# construct matrix with observationss

observed.matrix <- construct.observation.matrix(data=observed.abund,taxa=y.names$taxa )


# set-up list with parameter distributions:
# =========================================

# environmental input parameters
# ------------------------------

# reach and habitat specific parameters - environmental inputs
# w L T I0 fshade P N
par.env <- construct.envpars(Reaches,Habitats,env.data)

# choose exponential distribution for shading at M08 (0%)
#par.env[["M08_hab1_fshade"]] <- list(dist="Exponential",dist.par= 0.05)

# Invertebrates trait parameters
# ------------------------------


#file.db.food      = "databases/db_traits_fuzzy.dat"   # if they should be used than Cedrics stoich.cal assign_par_stoich_CM
file.db.feeding    = "databases/db_traits_tachet.dat"  #XXX
file.db.feeding    = "databases/db_freshwaterecology_20131206.dat"  #XXX
file.db.microhab   = "databases/db_traits_tachet.dat"  # AP
file.db.current    = "databases/db_traits_tachet.dat"
file.db.temp       = "databases/db_freshwaterecology_20131206.dat"
file.db.sapro      = "databases/db_freshwaterecology_20131206.dat"
#file.db.pH        = "databases/db_freshwaterecology_20131206.dat"
file.db.spear      = "databases/db_traits_rspear_Eurasia_20130114.dat"


# mid values of classes of environmental conditions used for traits
par.env.global <- derive.par.env.traits(file.db.temp=file.db.temp,
                                        file.db.current=file.db.current,
                                        file.db.sapro=file.db.sapro,
                                        #file.db.pH,
                                        file.db.spear=file.db.spear)  



#par.invtraits <- construct.invpars.traits(Invertebrates,file.db.bodymass="databases/db_bodymass.dat")

par.invtraits <- construct.invpars.traits(Invertebrates,
                                          file.db.feeding=file.db.feeding,
                                          file.db.microhab=file.db.microhab,
                                          file.db.current=file.db.current,
                                          file.db.temp=file.db.temp,
                                          file.db.sapro=file.db.sapro,
                                          file.db.spear=file.db.spear)

# other parameters:
# -----------------

# read in definitions for parameter distributions from file
par <- sysanal.read.distdef("input/parameter_input_Pfsee.dat")

# automatic generated parmaters for Organisms:


for ( i in c(Invertebrates,"Algae") )
{
  par[[paste(i,"_fgrotax",sep="")]]   <- par[["Organism_fgrotax"]]    # -
  par[[paste(i,"_fbasaltax",sep="")]] <- par[["Organism_fbasaltax"]]  # -
}

par.feedtype <- par.invtraits[grepl("feedingtype",names(par.invtraits))]

for ( i in c(Invertebrates) )
{
  par[[paste(i,"_hdens",sep="")]] <- par[["Invertebrates_hdens"]]   
  
  # test if pure filter feeders:
  
  if( sum(par.feedtype[paste(i,"feedingtype_aff",sep="_")]+par.feedtype[paste(i,"feedingtype_pff",sep="_")])>0 & 
        sum(par.feedtype[paste(i,"feedingtype_shr",sep="_")]+par.feedtype[paste(i,"feedingtype_min",sep="_")],
            par.feedtype[paste(i,"feedingtype_xyl",sep="_")]+par.feedtype[paste(i,"feedingtype_gra",sep="_")],
            par.feedtype[paste(i,"feedingtype_gat",sep="_")]+par.feedtype[paste(i,"feedingtype_pre",sep="_")])==0  
  )
  {
    par[[paste(i,"_Kfood",sep="")]] <- par[["Filt_Kfood"]] 
    cat(i ,"is pure filterfeeder, Kfood is set to",par[["Filt_Kfood"]],"\n" )
  } else {par[[paste(i,"_Kfood",sep="")]] <- par[["Invertebrates_Kfood"]]  }
}


#remove parameters that are no longer needed
ind_rm <- which(names(par)=="Organism_fgrotax"|names(par)=="Organism_fbasaltax"|
                  names(par)=="Invertebrates_Kfood"|names(par)=="Invertebrates_hdens"|
                  names(par)=="Filt_Kfood"
)
par <- par[-ind_rm]

# fixed stoichiometric parameters:

par.stoich.fix <- numeric(0)

# generate stoichiometric parameters that are fixed (mineralization, respiration,production):

# POM:
taxa.POM <- unique(y.names$y.taxa[y.names$y.groups=="POM"])
# mineralization:
par.stoich.fix[paste("Miner",taxa.POM,taxa.POM,sep="_")] <- -1

# Algae:
taxa.Algae <- unique(y.names$y.taxa[y.names$y.groups=="Algae"])
# respiration:
par.stoich.fix[paste("Resp",taxa.Algae,taxa.Algae,sep="_")]  <- -1
# production:
par.stoich.fix[paste("Prod",taxa.Algae,taxa.Algae,sep="_")]  <- 1

# Invertebrates:
taxa.Invertebrates <- unique(y.names$y.taxa[y.names$y.groups=="Invertebrates"])
# respiration:
par.stoich.fix[paste("Resp",taxa.Invertebrates,taxa.Invertebrates,sep="_")]  <- -1


# combine parameters

par.unc <- c(par,
             par.env,
             par.env.global,
             par.invtraits,
             par.stoich.fix)


divide.pars <- divide.parameters(par.unc=par.unc)


# parameters without uncertainty:

par.prior.delta <- divide.pars$prior.delta
names.par.prior.delta <- names(par.prior.delta)
par.prior.delta <- as.numeric(par.prior.delta)
names(par.prior.delta) <- names.par.prior.delta

write.table(par.prior.delta,paste("output/par_prior_delta_cerinp_",name.run,".dat",sep=""),
            sep="\t",col.names=FALSE,row.names=names(par.prior.delta))

# parameters with uncertainty

prior.par.defs <- divide.pars$prior.par.defs
prior.par.defs.untrans <- prior.par.defs

# # transform logarithmic priors to normals:
# 
# prior.par.defs <- transform.priors(prior.par.defs)

# parameters with uncertainty

par.samp <- generate.par.samp.matrix(n.samp=n.samp,par.unc=prior.par.defs)
# # convert strings to numerical codes:
# par.samp <- encode.envpars(par.samp)
# class(par.samp) <- "numeric"

write.table(t(par.samp),paste("output/prior_parsamp_cerinp_",name.run,".dat",sep=""),
            sep="\t",col.names=TRUE,row.names=FALSE)
