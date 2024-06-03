# set-up list with parameter distributions:
# =========================================

# environmental input parameters
# ------------------------------

# reach and habitat specific parameters - environmental inputs
# w L T I0 fshade P N
par.env <- construct.envpars(Reaches,Habitats,env.data)


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

write.table(par.prior.delta,paste(outputfolder,"/par_prior_delta_cerinp_",name.run,".dat",sep=""),
            sep="\t",col.names=FALSE,row.names=names(par.prior.delta))

# parameters with uncertainty

prior.par.defs <- divide.pars$prior.par.defs
prior.par.defs.untrans <- prior.par.defs

# # transform logarithmic priors to normals:
# 
# prior.par.defs <- transform.priors(prior.par.defs)

# parameters with uncertainty

cat("n.samp:",n.samp)
par.samp <- generate.par.samp.matrix(n.samp=n.samp,par.unc=prior.par.defs,verbose=TRUE)
# # convert strings to numerical codes:
# par.samp <- encode.envpars(par.samp)
# class(par.samp) <- "numeric"

write.table(t(par.samp),paste(outputfolder,"/prior_parsamp_cerinp_",name.run,".dat",sep=""),
            sep="\t",col.names=TRUE,row.names=FALSE)
