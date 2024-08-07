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
#Algae             <- paste("Alga",1:n.Algae,sep="")
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

write.table(par.invtraits,paste("output/par_invtraits_",name.run,".dat",sep=""),sep="\t",col.names=F,row.names=T)

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

par.prior.delta <- divide.pars$prior.delta
names.par.prior.delta <- names(par.prior.delta)
par.prior.delta <- as.numeric(par.prior.delta)
names(par.prior.delta) <- names.par.prior.delta

prior.par.defs <- divide.pars$prior.par.defs
prior.par.defs.untrans <- prior.par.defs

# transform logarithmic priors to normals:

prior.par.defs <- transform.priors(prior.par.defs)


# generate par samp to get a parameter sample at the mean of the marginals

parx <- generate.par.samp.matrix(n.samp=1,par.unc=par.unc)
# convert strings to numerical codes:

class(parx) <- "numeric"

# convert T in degC to K and CSusPOM to DSusPOM by multiplying with Filt_scope 
# if those parameters are uncertain:

par0 <- parx[,1]
#par0 <- convert.TdegC.CSusPOM(as.matrix(par0))  #rbind(as.matrix(par.prior.delta))
# Achtung, par0 enthält auch die par.prior.delta!

par.D <- as.vector(par0)
names(par.D) <- names(par0)

# parameters for likelihood function
# ----------------------------------

D.drift <- 0.05
p.obs   <- 0.5
p.abs   <- 0.1

# set presence/absence threshold depending on abundance
# D.crit <- calc.D.crit(y.names=y.names,par=par.D,par.prior.delta=NULL,observed.abund=observed.abund)
# write.table(t(D.crit),"D_crit_new.dat",sep="\t",col.names=T,row.names=F)


# par.stoich.out <- calc.stoich(par=as.list(par.D),returns="parout")
# par.stoich.taxa <- assign.par.stoich(par.invtraits,par.stoich.out,y.names)
# par.fix <- c(par.D,par.stoich.taxa)
# 
# pdf(paste("output/complete_foodwebs",name.run,"feedtypesFWB.pdf",sep=""),width=9,height=5,onefile=T)
# plot.foodweb(y.names,pars=par.fix,cex=1.1,font=2,title="complete foodweb",ncrit=8,
#              lcrit=8,lwd=3,bg=NA,lcol=colors()[555]) #,texts=F,pointcol=T)
# dev.off()



# transform parameters at the prior mean (=par.D) 
# -----------------------------------------------
par.unc.prior.mean <- rep(NA,length(prior.par.defs))
names(par.unc.prior.mean) <- names(prior.par.defs)

for(i in 1:length(prior.par.defs))
{
  if(grepl("transformed",names(prior.par.defs)[i] ))
  {
    ind <- which(names(par.D)==substr(names(prior.par.defs)[i],1,nchar(names(prior.par.defs)[i])-12))
    par.unc.prior.mean[i] <- log(par.D[ind])
  } else {
    ind <- which(names(par.D)==names(prior.par.defs)[i])
    par.unc.prior.mean[i] <- par.D[ind]
  }
}

# set-up inputs:
# --------------

inp <- NA