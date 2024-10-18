## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## --- Apply Streambugs to create artificial data ---
##
## --- May 2024 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PRELIMINARIES ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sys.setenv(LANG="EN")
seed <- 13
set.seed(seed) # Always set seed to a lucky number
getwd()        # Show working directory. It needs to be the location of 'main.r'
rm(list=ls())  # Free work space
graphics.off() # Clean graphics display

## Load libraries ####
# ~~~~~~~~~~~~~~~

if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") }               # to sort, join, merge data
if ( !require("tidyr") ) { install.packages("tidyr"); library("tidyr") }               # to sort, join, merge data
# if ( !require("tidyverse") ) { install.packages("tidyverse"); library("tidyverse") }               # to sort, join, merge data
if ( !require("ggplot2") ) { install.packages("ggplot2"); library("ggplot2") }         # to do nice plots
# if ( !require("plotly") ) { install.packages("plotly"); library("plotly") }            # to do nice "interactive" plots
# if ( !require("sf") ) { install.packages("sf"); library("sf") }                        # to read layers for plotting maps
if ( !require("pracma") ) { install.packages("pracma"); library("pracma") }            # to do polynomial interpolation (of ecological traits)                # to read layers for plotting maps
if ( !require("stringr") ) { install.packages("stringr"); library("stringr") }                        # to handle character vectors
if ( !require("parallel") ) { install.packages("parallel"); library("parallel") }                        # to run things on parallel on the server

## Directory and file definitions ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir.inputs      <- "../A1_inputs_streambugs/"
dir.outputs     <- "../A3_outputs_streambugs/"
dir.utilities   <- "../00_utilities/"

file.midat.all    <- "All_2339samples_134taxa_13envfact_ProcessedDataset.csv" # output of: PreprocessSwissInvertebrateEnvironmentalModellingData
file.midat.bdm    <- "BDM_950samples_294taxa_13envfact_ProcessedDataset.csv" # output of: PreprocessSwissInvertebrateEnvironmentalModellingData
file.taxonomy     <- "invertebrates_taxonomy_2023-08-03.dat"  # output of: PreprocessSwissInvertebrateEnvironmentalModellingData

# file.par.update   <- "Parameters_update_Vermeiren_Chollet.dat"
file.par.update   <- "T1d_BDM_CH_Tvsappestsubst_maxpost_trait_pars_2023-12-15.dat"
file.par.correc   <- "correction_preference_traits.csv"
file.synth.points <- "Synthetic_environmental_data_2024-03-19.dat" # output of: PreprocessSwissInvertebrateEnvironmentalModellingData
file.selected.taxa      <- "selected_taxa_analysis.csv"

source(paste0(dir.utilities, "utilities_global.r"))
source("utilities_streambugs.r")
source("functions_run_streambugs.r")
source("library/application specific NIS/load_library.r")

# prepare inputs for plotting results on swiss maps
# map.inputs        <- map.inputs(directory = paste0(dir.utilities,"swiss.map.gdb"))

## Analysis options ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

server <- T

# run simulations to produce ICE
ICE <- F # if False, simulations run for whole dataset
no.sites.ice <- 50
# custom.range <- c(5:10)
# no.steps.ice <- length(custom.range)
no.steps.ice <- 40
env.fact.ice <- "tempmaxC"

# global options
run.C            <- T # run C version of the rhs of Streambugs
return.res.add   <- F # calculate additional output (rates, limitation factors) when running Streambugs
plot.biomass     <- F # plot time evolution of biomass (heavier computation for dataframes and plots if many sites and many taxa)
no.class.new     <- 30 # number of new classes/points for preference trait interpolation

# has to be defined globally but is used in fct to construct parameters
inp <- NA # set NA for time-dependent input and parameters

# set-up output times

# # sequenced output times of different timesteps length
# n.years1       <- 2 # number of years
# n.days.appart1 <- 1/10  # number of days from one time step to another (can be smaller than one day, e.g., 0.5)
# tout1          <- c(seq(0, n.years1, by = n.days.appart1/365)) # set up time vector of [n.years*365/n.days.appart] steps
# n.years2       <- 10 # number of years
# n.days.appart2 <- 1  # number of days from one time step to another (can be smaller than one day, e.g., 0.5)
# tout2          <- c(seq(n.years1, n.years2, by = n.days.appart2/365)) # set up time vector of [n.years*365/n.days.appart] steps
# tout <- c(tout1, tout2)
# n.years <- n.years2

# tout <- c(#seq(0,5,by=0.01),
#           seq(0,3,by=1/10/365),
#           seq(3.0005,10,by=1/2/365),
#           seq(10.003,n.years,by=1/365)) # t steadystate:170-200

# continuous output times of same timesteps length
n.years       <- 100 # number of years
n.days.appart <- 1  # number of days from one time step to another (can be smaller than one day, e.g., 0.5)
tout          <- c(seq(0, n.years, by = n.days.appart/365)) # set up time vector of [n.years*365/n.days.appart] steps
# tout <- c(0, 0.001, 0.002)

# select method for ode solver
method <- "rk4"
# method <- "lsoda"

n.time.steps  <- length(tout) # number of time steps simulated
cat("Simulating", n.years, "years", # "every", n.days.appart, "day,", 
    " total:", n.time.steps, "time steps")

# catchment scale variable by which we select the taxa pool (ask Rosi Siber for different catchments delimitation)
catch.variable <- "Watershed" 
# catch.variable <- "Watershed_Reg"

# tune parameters to get stronger presence/absence response to environmental factors, necessary for our application
# if Flag True, then Value is applied to the tuned parameter
par.adjust <- list("curve.curr.temp"     = list("Flag" = F, "Value" = 0),  # assign Value (e.g., -10) to curve parameter of limitation factor of current and temperature, to get stronger (pres/abs) taxa
                   "curve.orgmic.sapro"  = list("Flag" = F, "Value" = 10),   # assign Value (e.g., -10) to curve parameter of limitation factor of orgmicropoll and saproby, to get stronger (pres/abs) taxa
                   "interc.orgmic.sapro" = list("Flag" = F, "Value" = 4),    # assign Value (e.g., -10) to curve parameter of limitation factor of orgmicropoll and saproby, to get stronger (pres/abs) taxa
                   "sens.anal.init.cond" = list("Flag" = F, "Value" = 0.5),  # do a sensitivity analysis of results to initial conditions (initial densities, original 1gDM) multiplied by Value
                   "hdens"               = list("Flag" = F, "Value" = 0.05), # multiply hdens by small number (e.g., 0.05) to get strong self-inhibition, main impact on abundance and not pres/abs 
                   "thresh.indiv"        = list("Flag" = F, "Value" = 1),    # define presence data points above specific threshold Value (e.g., 1) of number of individuals
                   "update.traits"       = list("Flag" = T, "Value" = ".hybrid"))

correct.shade        <- 0.6 # factor to correct shading (fshade) for now, otherwise algae is too limited

# catchment selection for simulation (can be specific or all to analyze taxa pool and results)
select.catch <- "all"
# select.catch <- "Aare"
# select.catch <- "Doubs"
# select.catch <- "Limmat"
# select.catch <- "Reuss"
# select.catch <- "RheinabBS"
# select.catch <- "Rhone"
# select.catch <- "RhoneLeman"
# select.catch <- "Ticino"
# select.catch <- "Rhein"

# select specific site(s) for debugging and analysis
select.sites   <- "random"
# select.sites <- "CH114BEAa"
# select.sites <- "SynthPoint2533Ti"
# select.sites <- "SynthPoint6969Rh"
# select.sites <- c("SynthPoint88Rh", "SynthPoint21Rh") # Rhone
# select.sites   <- "CH076ZGRe" # e.g., site with problem in Reuss (when running all reaches at once, 10 sites, 3 years)

# select number of sites per catchment for simulation
if(ICE){
  n.sites.per.catch    <- no.sites.ice * no.steps.ice
} else {
  # n.sites.per.catch <- "all"
  n.sites.per.catch <- 2
}
sites.selection      <- list("n.sites" = n.sites.per.catch, 
                             "select.sites" = select.sites)

# specify how we select the taxa pool
# metric: based on prevalence or number of presence points in the catchment
# threshold: threshold above which taxa is selected
taxa.selection <- c("metric" = "Prevalence", 
                    "threshold" = 0.1)
# taxa.selection <- c("metric" = "NbPresPoints", 
#                     "threshold" = 7)

# select taxonomic levels to filter taxa pool
select.taxonomy = c("species", "genus", "family") 

# select colors for present/absent plots
vect.col.pres.abs <- c( "Present" = "#00BFC4", "Absent" = "#F8766D")
scales::show_col(vect.col.pres.abs)

# define run name for all output files
name.par.adjust <- ""
# for (i in 1:length(par.adjust)) {
#   name.par.adjust <- paste0(name.par.adjust, ifelse(par.adjust[[i]][["Flag"]], paste0(names(par.adjust)[i], par.adjust[[i]][["Value"]], "_"), ""))
# }

name.run <- paste0(
  ifelse(ICE, paste0("ICE_", no.steps.ice, "Steps_"), ""),
  "PolyInterp", no.class.new, "_",
  ifelse(run.C, "runC_", "runR_"),
  name.par.adjust,
  # method, "Met_",
  # "CorrectPar_",
  ifelse(return.res.add, "ResAdd_", ""),
  # "Prev", taxa.selection[2], "_",
  # "Shade", correct.shade, "_",
  # n.sites, "Sites_",
  n.years, "Yea_",
  n.time.steps, "Steps_" #,
)
print(name.run)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DATA PROCESSING ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# read data files
data.env.midat  <- read.csv(paste0(dir.inputs, file.midat.all))
data.env.synth  <- read.table(paste(dir.inputs,file.synth.points,sep=""), header=T, sep="\t", stringsAsFactors=F)
data.inv.bdm    <- read.csv(paste0(dir.inputs, file.midat.bdm))
data.taxonomy   <- read.table(paste(dir.inputs,file.taxonomy,sep=""),header=T,sep="\t", stringsAsFactors=F)
data.par.update <- read.table(paste(dir.inputs,file.par.update,sep=""),header=T,sep="\t", stringsAsFactors=F)
data.par.correc <- read.csv(paste0(dir.inputs, file.par.correc), sep = ";")
data.selected.taxa  <- read.csv(paste0(dir.utilities, file.selected.taxa), header = T, sep = ";", stringsAsFactors = F)

## Environmental data ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### preprocess dataset ####

# remove duplicated sites in data from midat
colnames(data.env.midat)[which(colnames(data.env.midat) == "SiteId")] <- "ReachID" # change column name of site id to match streambugs
data.env.midat <- data.env.midat[!duplicated(data.env.midat$ReachID), ] # remove duplicated sites

# bind env data from midat and synthetic points
colnames(data.env.synth)[which(colnames(data.env.synth) == "SiteId")] <- "ReachID" # change column name of site id to match streambugs
rind <- sample(1:dim(data.env.synth)[1], 2300) # subselect randomly some synthetic points otherwise too many of them
data.env.synth <- data.env.synth[rind,]
data.env.synth$MonitoringProgram <- "SyntheticPoint" # specify "fake" monitoring program
data.env.all <- bind_rows(data.env.midat, data.env.synth)
remove(data.env.synth) # remove entire synthetic sample points data from environement because slows down the system
data.env.all <- filter(data.env.all, tempmaxC > 4) # remove sites with temperature bellow 4 degrees
data.env.all$Watershed <- gsub("Doux", "Doubs", data.env.all$Watershed) # correct cute mistake


# clean ReachID names
data.env.all$ReachID <- gsub("CSCF_", "", data.env.all$ReachID) # remove some characters in ReachID names
data.env.all$ReachID <- gsub("_", "", data.env.all$ReachID) # remove some characters in ReachID names
data.env.all$ReachID <- gsub("-", "", data.env.all$ReachID) # remove some characters in ReachID names
data.env.all$ReachID <- gsub("\\.", "", data.env.all$ReachID) # remove some characters in ReachID names
data.env.all$ReachID <- paste0(data.env.all$ReachID, substr(data.env.all[,catch.variable], start = 1, stop = 2)) # add first two letters of catch at the end of ReachID

# correct some environmental factors
data.env.all$shade <- correct.shade*data.env.all$shade # correct shade to still have algae growing
data.env.all$w <- 1 # set width 1 m to have gDM/m as output

# clean column names of taxa
cind.taxa <- which(grepl("Occurrence.", colnames(data.env.all)))
colnames(data.env.all)[-cind.taxa] <- gsub("\\.", "_", colnames(data.env.all)[-cind.taxa]) # adapt column names otherwise its creates problems when running streambugs

### select inputs ####

# select environmental inputs
vect.info.inputs <- c("ReachID", "Habitat", "Reachname")
vect.info.plots <- c("X", "Y", catch.variable, "MonitoringProgram")
vect.env.inputs <- c("Temp_average",                    # used in metabolic rate, calculated as yearly mean temperature with catchment height and area
                     "L", "w", "I0", "C_SusPOM",        # fixed for all sites
                     "shade", "Lit_Inp" ,               # drive algae growth and food for shredders, calculated with fraction of forest (deciduous) on the river bank
                     "C_P", "C_N",                      # drive algae growth, calculated with fraction of different land use in the catchment
                     "currentms", "tempmaxC",           # used with invertebrate preference trait for limiting growth rate, calculated with catchment and river attributes (height, area, slope, discharge) 
                     "orgmicropollTU", "saprowqclass",  # used with invertebrate preference trait for increasing death rate, calculated with catchment land use attributes
                     "microhabaf_type1", "microhabaf_type2", "microhabaf_type3", "microhabaf_type4" # assigned with a "random" number (1/4 = 0.25), not needed for our application
)
names(vect.env.inputs) <- vect.env.inputs
vect.all.inputs <- c(vect.info.inputs, vect.info.plots, vect.env.inputs)
names(vect.all.inputs) <- vect.all.inputs

# subselect relevant direct and indirect factors for analysis
vect.dir.env.fact <- vect.env.inputs[c(10:13)]
vect.indir.env.fact <- vect.env.inputs[c(6:9)]
vect.analysis.env.fact <- c(vect.dir.env.fact, vect.indir.env.fact)
print(vect.analysis.env.fact)

# set Na, "random", or specific number for env. inputs not needed for this application
for (factor in vect.all.inputs[-which(vect.all.inputs %in% colnames(data.env.all))]) {
    if(grepl("type", factor)){
        data.env.all[,factor] <- 0.25 # xxxnis
    } else {
        data.env.all[,factor] <- "random" # set to "random" for any other empty class
    }
}
cind.env.inp <- which(colnames(data.env.all) %in% vect.env.inputs)
rind <- which(complete.cases(data.env.all[,cind.env.inp])) # rows without NA in environmental input
cind.all.inp <- which(colnames(data.env.all) %in% vect.all.inputs)
data.env.inputs <- data.env.all[rind,c(cind.all.inp, cind.taxa)]

### plot hist and maps ####
 
# # make dataframe with combination of direct factors for plotting later
# df.comb.dir.env.fact <- t(combn(vect.dir.env.fact,2))
# 
# # histograms
# if(select.catch == "all"){
#     plot.data <- data.env.inputs[, c("ReachID", "X", "Y", vect.analysis.env.fact)]
# } else {
#     temp.data.env <- data.env.inputs[which(data.env.inputs[,catch.variable] == select.catch), ]
#     plot.data <- temp.data.env[, c("ReachID", "X", "Y", vect.analysis.env.fact)]
# }
# length(vect.analysis.env.fact)
# par(mfrow=c(3,3))
# for (fact in vect.analysis.env.fact) {
#     hist(plot.data[,fact], main = fact, xlab = fact)
# }
# 
# # plot on swiss map
# file.name <- paste0("MapEnvFact_", dim(plot.data)[1], "Sites")
# if(!file.exists(paste0(dir.outputs, file.name, ".pdf"))){
#     for(fact in vect.analysis.env.fact){
#       fact <- vect.analysis.env.fact[1]
#       plot.data[,fact] <- round(as.numeric(plot.data[,fact]), digits = 2)
#     }
#   # !! problems with x/Y coordinates
#     list.plots <- maps.env.fact(map.inputs = map.inputs, 
#                                 vect.env.fact = vect.analysis.env.fact, vect.info = c("X", "Y", "ReachID"),
#                                 data = plot.data)
#     print.pdf.plots(list.plots = list.plots, width = 15, height = 8,
#                     dir.output = dir.outputs, file.name = file.name)
# }

## Invertebrate data ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get taxa for analysis (in global utilities)
selected.taxa.analysis <- data.selected.taxa$Occurrence.taxa
names(selected.taxa.analysis) <- gsub("Occurrence.", "", selected.taxa.analysis)

### select catchments ####

# clean column names taxa
cind.taxa <- which(grepl("Occurrence.", colnames(data.inv.bdm)))
colnames(data.inv.bdm)[-cind.taxa] <- gsub("\\.", "_", colnames(data.inv.bdm)[-cind.taxa]) # adapt column names
cind.taxa <- which(grepl("Occurrence", colnames(data.inv.bdm))) 
colnames(data.inv.bdm)[cind.taxa] <- gsub("_", "", colnames(data.inv.bdm)[cind.taxa]) # change species name to remove "_"

# select catchments
table.catch.all <- table(data.env.inputs[,catch.variable]) # analyse distribution of sites per catchment in the whole dataset
table.catch.all
vect.catch.all <- names(table.catch.all)
table.catch.bdm <- table(data.inv.bdm[,catch.variable]) # analyse distribution of sites per catchment in the BDM dataset
table.catch.bdm
vect.catch.select <- vect.catch.all[!is.na(match(vect.catch.all, names(table.catch.bdm)))] # remove catchments that don't have BDM sites
too.small.catch <- names(which(table.catch.all < 30)) # remove catchments that have too few sites
vect.catch.select <- vect.catch.select[-which(vect.catch.select %in% too.small.catch)] 
table.catch.all[vect.catch.select] # check final selection of catchments

# subselect sites in selected ctachments
data.env.inputs <- filter(data.env.inputs, data.env.inputs[,catch.variable] %in% vect.catch.select)

### produce ICE data ####

if(ICE){
  summary(as.factor(data.env.inputs$MonitoringProgram))
  # rind.bio.data <- which(data.env.inputs$MonitoringProgram != "SyntheticPoint")
  # rind.ice.site <- sample(rind.bio.data, no.sites.ice)
  # if(sites.selection["select.sites"] == "random"){
  #   if(sites.selection["n.sites"] != "all" & as.numeric(sites.selection["n.sites"]) < dim(env.data)[1]){
  #     env.data <- env.data[1:sites.selection["n.sites"],] # temporary take only few sites for trial
  #   }
  # } else {
  #   rind.sites.selected <- which(grepl(sites.selection["select.sites"], data.env.inputs$ReachID))
  #   if(sites.selection["n.sites"] != length(rind.sites.selected)){
  #     rind.sites.selected <- c(1:(n.sites - length(rind.sites.selected)), rind.sites.selected) 
  #   }
  #   env.data <- env.data[rind.sites.selected,] # temporary take only few sites for trial
  # }
  if(!("random" %in% sites.selection[["select.sites"]])){
    rind.select.site <- which(data.env.inputs$ReachID %in% sites.selection[["select.sites"]])
    if(no.sites.ice != length(rind.select.site)){
      rind.ice.site <- c(sample(1:nrow(data.env.inputs), no.sites.ice - length(rind.select.site)), rind.select.site)
    } else {
      rind.ice.site <- rind.select.site
    }
  } else {
    rind.ice.site <- sample(1:nrow(data.env.inputs), no.sites.ice)
  }
  data.base.ice <- data.env.inputs[rind.ice.site,]
  
  # get range of env.factor
  min.env.fact <- min(data.env.inputs[,env.fact.ice])
  max.env.fact <- max(data.env.inputs[,env.fact.ice])
  
  # create range of all env.factor
  if(exists("custom.range")){
    vect.env.fact.ice <- custom.range
  } else {
    range.size <- max.env.fact-min.env.fact
    step.size  <- range.size/(no.steps.ice-1)
    vect.env.fact.ice <- seq(from=min.env.fact,
                             to=max.env.fact,
                             by=step.size)
  }

  data.ICE <- data.frame()
  
  for (i in 1:no.sites.ice) {
    # i <- 1
    vect.select.site <- data.base.ice[i,]
    temp.data.ICE <- vect.select.site[rep(1, each = no.steps.ice),]
    temp.data.ICE[,env.fact.ice] <- vect.env.fact.ice
    temp.data.ICE$ICE_id <- seq(1, no.steps.ice, by = 1)
    temp.data.ICE$SiteID <- temp.data.ICE$ReachID
    temp.data.ICE$ReachID <- paste0(temp.data.ICE$SiteID, temp.data.ICE$ICE_id)
    data.ICE <- bind_rows(data.ICE, temp.data.ICE)
  }
  
  data.env.inputs.orig <- data.env.inputs
  data.env.inputs <- data.ICE
  vect.catch.select.orig <- vect.catch.select
  vect.catch.select <- unique(data.ICE$Watershed)
  
  file.name <- paste0(dir.outputs,"WideDataInputs_ICE_",
                      no.sites.ice,"Sites",".csv")
  write.table(data.base.ice, file = file.name, sep =";", row.names = F)
}

if(ICE){
  n.sites <- no.sites.ice
} else {
  n.sites <- ifelse(n.sites.per.catch == "all", 
                    dim(data.env.inputs)[1], 
                    length(vect.catch.select) * n.sites.per.catch) # compute nb total of sites for catchments selected
}
# vect.catch.select <- vect.catch.select[-which(vect.catch.select == "Doux")] # exemple on how to remove a catchment

# modify name run (for all files name) with selected catchments
if(select.catch == "all"){ 
    name.run <- paste0(length(vect.catch.select), "Catch_",
                       n.sites, "Sites_",
                       name.run)
} else {
    vect.catch.select.orig <- vect.catch.select
    vect.catch.select <- select.catch 
    name.run <- paste0(vect.catch.select, "_", 
                       n.sites.per.catch, "Sites_",
                       name.run)
}
names(vect.catch.select) <- vect.catch.select
print(name.run)

# construct dataset with taxonomic level and prevalence per catchment for each taxon
list.df.prev.pres <- get.df.prev.pres.taxonomy(data.inv = data.inv.bdm, data.taxonomy = data.taxonomy,
                                          catch.variable = catch.variable, vect.catch.select = vect.catch.select, dir.outputs = dir.outputs)

data.prevalence <- list.df.prev.pres$Prevalence
data.nb.pres <- list.df.prev.pres$NbPresPoints
all.taxa <- data.prevalence$Taxon

if(taxa.selection["metric"] == "Prevalence"){
    data.taxa.selection <- data.prevalence
} else if(taxa.selection["metric"] == "NbPresPoints"){
    data.taxa.selection <- data.nb.pres
}

### get parameters for all taxa ####

# rename colnames (trait preference) of updated parameters from Vermeiren et al., 2021
data.par.update <- data.par.update %>%
  rename(
    "Taxon" = X,
    "tempmaxtolval_class1" = vco,
    "tempmaxtolval_class2" = cod,
    "tempmaxtolval_class3" = mod,
    "tempmaxtolval_class4" = war,
    "currenttolval_class1" = st,
    "currenttolval_class2" = slow,
    "currenttolval_class3" = mod.1,
    "currenttolval_class4" = high,
    "saprotolval_class0"   = xeno,
    "saprotolval_class1"   = oligo,
    "saprotolval_class2"   = beta,
    "saprotolval_class3"   = alpha,
    "saprotolval_class4"   = poly
    )

list.par.update <- list("update" = data.par.update, "correc" = data.par.correc)

# extract trait scores and mass of all taxa
file.name <- paste0("df.preferences_", "PolyInterp", no.class.new, "_",
    ifelse(par.adjust[["round.inv.traits"]][["Flag"]], "round.inv.traits", 
           ifelse(par.adjust[["update.traits"]][["Flag"]], paste0("update", par.adjust[["update.traits"]][["Value"]]), "orig")), 
    "_")
file.rds.name <- paste0(file.name, "266Taxa.rds")

if(file.exists(paste0(dir.outputs, file.rds.name))){
    list.df.preferences <- readRDS(paste0(dir.outputs, file.rds.name))
    mass.all.inv <- list.df.preferences[["mass"]]
    vect.all.inv <- names(mass.all.inv)
} else {
    list.par.all.taxa <- construct.variables.par.catch(catch = taxa.selection["metric"], data.env.inputs = data.env.inputs, 
                                                       list.par.update = list.par.update, par.adjust = par.adjust, no.class.new = no.class.new,
                                                       data.taxa.selection = data.taxa.selection, taxa.selection = taxa.selection, selected.taxa.analysis = selected.taxa.analysis,
                                                       sites.selection = sites.selection, select.taxonomy = select.taxonomy,
                                                       catch.variable = catch.variable, name.run = name.run,
                                                       dir.inputs = dir.inputs, dir.outputs = dir.outputs)
    # list.par.all.taxa$par.fixc[which(grepl("Limoniidae", names(list.par.all.taxa$par.fixc)) & grepl("class", names(list.par.all.taxa$par.fixc)))]
    system.def <- streambugs.get.sys.def(y.names=list.par.all.taxa$y.names$y.names, par=list.par.all.taxa$par.fixc)
    vect.all.inv <- list.par.all.taxa$Invertebrates
    
    # write updated and extracted parameters in csv files for all taxa
    list.df.preferences <- write.df.preferences(system.def = system.def,
                                                Invertebrates = vect.all.inv,
                                                file.name = file.name,
                                                dir.outputs = dir.outputs)
    mass.all.taxa <- system.def$par.taxaprop.direct$parvals[,"M"] # mean individual biomass of taxon i
    for (taxon in list.par.all.taxa$y.names$taxa) {
        names(mass.all.taxa)[which(grepl(taxon, names(mass.all.taxa)))] <- taxon
    }
    mass.all.taxa <- mass.all.taxa[which(!duplicated(names(mass.all.taxa)))]
    mass.all.inv <- mass.all.taxa[-c(1:3)]
    list.df.preferences[["mass"]] <- mass.all.inv
    saveRDS(list.df.preferences, file = paste0(dir.outputs, 
                                               file.name, length(vect.all.inv), "Taxa", ".rds"))
}

# make long dataframe trait preference
exp.trans.curve <- c(-5, -10)
long.df.trait <- data.frame()
for (taxon in vect.all.inv) {
    # taxon <- vect.all.inv[1]
    # print(taxon)
    for (factor in vect.dir.env.fact) {
        # factor <- vect.dir.env.fact[1]
        # print(factor)
        if(grepl("temp", factor)){
            df.pref <- list.df.preferences$tempmaxtolval
        } else if(grepl("current", factor)){
            df.pref <- list.df.preferences$currenttolval
        } else if(grepl("orgmicro", factor)){
            df.pref <- list.df.preferences$orgmicropolltolval
        } else if(grepl("sapro", factor)){
            df.pref <- list.df.preferences$saprotolval
        }
        temp.df <- data.frame(Value = df.pref$Values)
        temp.df$Pref <- df.pref[, taxon]
        temp.df$Factor <- factor
        temp.df$Taxon <- taxon
        for (curve in exp.trans.curve) {
          # curve <- exp.trans.curve[1]
          temp.df[,paste0("ExpTrans_", curve)] <- exp.transform(temp.df$Pref, intercept = 0, curv = curve)
        }
        
        long.df.trait <- rbind(long.df.trait, temp.df)
    }
}

file.name <- paste0("PreferenceTraitAllTaxa_MixUpdateOrig_PolyInterp", no.class.new, ".pdf")
if(!file.exists(paste0(dir.outputs, file.name))){
  list.plots <- list()
  for(taxon in vect.all.inv){
    # taxon <- vect.all.inv[5]
    temp.df.trait <- filter(long.df.trait, Taxon == taxon)
    plot.data <- gather(temp.df.trait, key = "Transformation", value = "Score", -c("Value", "Factor", "Taxon"))
    p <- ggplot(plot.data, aes(x = Value, y = Score, color = Transformation))
    p <- p + geom_point()
    p <- p + geom_line()
    p <- p + facet_wrap(. ~ Factor, scales = "free_x")
    p <- p + scale_y_continuous(limits = c(0,1))
    p <- p + labs(title = taxon)
    # ggplotly(p)
    list.plots[[taxon]] <- p
  }
  print.pdf.plots(list.plots = list.plots, width = 5, height = 5, dir.output = dir.outputs, info.file.name = "", file.name = file.name,
                  png = F)
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SIMULATIONS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Run streambugs ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(name.run)
start.all <- proc.time()

if (server){
    
    # construct variables and parameters per catchment
    list.variables.par.catch <- mclapply(vect.catch.select, mc.cores = length(vect.catch.select), FUN = construct.variables.par.catch,
                                       data.env.inputs = data.env.inputs, 
                                       list.par.update = list.par.update, par.adjust = par.adjust, no.class.new = no.class.new,
                                       data.taxa.selection = data.taxa.selection, 
                                       taxa.selection = taxa.selection, 
                                       selected.taxa.analysis = selected.taxa.analysis,
                                       sites.selection = sites.selection, 
                                       select.taxonomy = select.taxonomy, 
                                       catch.variable = catch.variable, 
                                       plot.foodweb = T, name.run = name.run,
                                       dir.inputs = dir.inputs, dir.outputs = dir.outputs)
    
    # streambugs.debug = 1
    # # print metadata
    # sink(file = paste0(dir.outputs, name.run, "debug.txt"))
    
    # run streambugs
    list.results <- mclapply(list.variables.par.catch, mc.cores = length(vect.catch.select), FUN = run.streambugs.catch,
                           tout = tout, return.res.add = return.res.add, name.run = name.run, 
                           dir.output = dir.outputs, run.C = run.C,
                           write.plot.results = F,
                           method = method
                           # atol = 1e-06,
                           # rtol = 1e-06
    )# ,
    # hmax = 0.1)
} else {
    
    # construct variables and parameters per catchment
    list.variables.par.catch <- lapply(vect.catch.select, FUN = construct.variables.par.catch,
                                       data.env.inputs = data.env.inputs, 
                                       list.par.update = list.par.update, par.adjust = par.adjust, no.class.new = no.class.new,
                                       data.taxa.selection = data.taxa.selection, 
                                       taxa.selection = taxa.selection, 
                                       selected.taxa.analysis = selected.taxa.analysis,
                                       sites.selection = sites.selection, 
                                       select.taxonomy = select.taxonomy, 
                                       catch.variable = catch.variable, 
                                       plot.foodweb = T, name.run = name.run,
                                       dir.inputs = dir.inputs, dir.outputs = dir.outputs)
    
    # streambugs.debug = 1
    # # print metadata
    # sink(file = paste0(dir.outputs, name.run, "debug.txt"))

    # run streambugs
    list.results <- lapply(list.variables.par.catch, FUN = run.streambugs.catch,
                           tout = tout, return.res.add = return.res.add, name.run = name.run, 
                           dir.output = dir.outputs, run.C = run.C,
                           write.plot.results = F,
                           method = method
                           # atol = 1e-06,
                           # rtol = 1e-06
                           )# ,
                           # hmax = 0.1)
    
}

# sink()
duration.all <- proc.time() - start.all
cat(duration.all,"\n")

# print metadata
sink(file = paste0(dir.outputs, name.run, "metadata.txt"))
for (i in 1:length(list.results)) {
    # list.results$Aare$list.metadata.catch
    print(list.results[[i]][["list.metadata.catch"]])
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
}
cat("The entire simulation time is:", duration.all)
sink()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PLOTS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Time series biomass ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract vectors of all sites and taxa modeled
vect.sites <- unique(unlist(sapply(list.results, "[[", 1)["reaches",]))
vect.taxa <- unique(unlist(sapply(list.results, "[[", 1)["taxa",]))
vect.inv <- as.vector(unique(unlist(sapply(list.results, "[[", "Invertebrates"))))
vect.occ.inv <- paste("Occurrence.", vect.inv, sep = "")
names(vect.inv) <- vect.occ.inv

temp.vect.taxa <- names(selected.taxa.analysis)[names(selected.taxa.analysis) %in% vect.taxa]
# temp.vect.taxa <- vect.taxa[1:3]

# plot.biomass <- T
if(plot.biomass){
    
    # make list of long dataframe biomass for each taxon
    list.df.res.taxa <- list()
    
    for (taxon in temp.vect.taxa) {
        # taxon <- temp.vect.taxa[2]
        long.df.taxon.biomass <- data.frame()
        for (n in 1:length(list.results)) { # loop over catchments modeled
            # n <- 1
            temp.all.res <- list.results[[n]]
            temp.matrix.res <- as.matrix(temp.all.res$res)
            temp.sites <- temp.all.res$y.names$reaches
            if(taxon %in% temp.all.res$y.names$taxa){
                c.ind <- which(grepl(taxon, colnames(temp.matrix.res)))
                temp.df <- temp.matrix.res[,c(1, c.ind)]
                for (site in temp.sites) {
                    # site <- temp.sites[2]
                    c.ind2 <- which(grepl(paste0(site, "_"), colnames(temp.df)))
                    temp.df2 <- as.data.frame(temp.df[,c(1,c.ind2)])
                    if(any(is.na(temp.df2[,2]))){ # print warning if NA in results
                        cat("\nWARNING NA in results:\nTaxon:", taxon,
                            "\nSite:", site,
                            "\nTime step:",  min(which(is.na(temp.df2[,2]))), "\n")
                    }
                    colnames(temp.df2)[2] <- "Value"
                    temp.df2$Factor <- "Biomass"
                    temp.df2$Taxon <- taxon
                    temp.df2$Site <- site
                    temp.df2$FactorType <- "Biomass"
                    long.df.taxon.biomass <- rbind(long.df.taxon.biomass, temp.df2)
                }
            }
        }
        
        # long.df.taxon.biomass <- filter(long.df.taxon.biomass, Site %in% "CH070TGRh") # to plot one selected site only
        list.df.res.taxa[[taxon]] <- long.df.taxon.biomass
    }
    
    # Plot time series biomass
    
    # plot each site separately in a grid for each taxon 
    list.plots.sep.sites <- lapply(list.df.res.taxa, FUN = function(plot.data){
        # plot.data <- list.df.res.taxa[[1]]
        taxon <- unique(plot.data$Taxon)
        p <- ggplot(plot.data, aes(x = time, y = Value))
        p <- p + geom_line(size=1)
        p <- p + facet_wrap(. ~ Site)
        p <- p + labs(title = taxon) + ylab("Biomass [gDM/m2]")
        if(!grepl("POM", taxon) & !grepl("Algae", taxon)){
            p <- p + ylim(0,1)
        }    
        p <- p + theme_bw()
        p <- p + theme(legend.position="none")
        return(p)
    })
    # list.plots.sep.sites[[2]]
    
    
    # MODIFY BACK LATER
    # plot all sites together for each taxon 
    list.plots.all.sites <- lapply(list.df.res.taxa, FUN = function(plot.data){
        # plot.data <- list.df.res.taxa[[3]]
      # plot.data$Site <- as.numeric(gsub(select.sites, "", plot.data$Site))
        taxon <- unique(plot.data$Taxon)
        q <- ggplot(plot.data, aes(x = time, y = Value, color = Site))
        q <- q + geom_line(size = 1, alpha = 0.7)
        q <- q + labs(title = taxon) + ylab("Biomass [gDM/m2]")
        if(!grepl("POM", taxon) & !grepl("Algae", taxon)){ # fix y axis for invertebrates
            q <- q + ylim(0,1)
        }
        # q <- q + scale_colour_gradient(low = "blue", high = "red")
        # q <- q + xlim(9,10)
        q <- q + theme_bw()
        # if(length(unique(long.df.taxon.biomass$Site)) > 10){
            # q <- q + theme(legend.position="none")
        # }
        return(q)
    })
    # list.plots.all.sites[[2]]
    # ggplotly(list.plots.all.sites[[3]])
    
    # print pdf with plots
    # file.name <- paste0("_IndivBiomass_ggplot_POM_gridSites")
    file.name <- paste0("_IndivBiomass_ggplot_allTaxa_gridSites")
    print.pdf.plots(list.plots = list.plots.sep.sites, width = 15, height = 12, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name,
                    png = F)
    file.name <- paste0("IndivBiomass_ggplot_allTaxa_allSites")
    # file.name <- paste0("IndivBiomass_ggplot_POM_allSites")
    # paste0(name.run, file.name)
    # file.name <- paste0("IndivBiomass_ggplot_allTaxa_", select.sites)
    # file.name <- paste0("IndivBiomass_ggplot_allTaxa_", "CH070TGRh")
    
    print.pdf.plots(list.plots = list.plots.all.sites, width = 8, height = 6, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name,
                    png = F)
}


# Env data of modeled sites ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract env.data selected for the simulation
data.env.res <- data.frame()
for (n in 1:length(list.results)) {
    # n <- 1
    temp.env.data <- list.results[[n]][["env.data"]][,vect.all.inputs]
    data.env.res <- rbind(data.env.res, temp.env.data)
}
# data.env.res <- left_join(data.env.res, data.env.inputs[,c("ReachID", "X", "Y", "MonitoringProgram", catch.variable)], by = "ReachID")

# # hist distribution selected env. data
# par(mfrow=c(3,3))
# for (fact in vect.analysis.env.fact) {
#     hist(data.env.res[,fact], main = fact, xlab = fact)
# }

# plot selected env data multidmensional
plot.data <- data.env.res
plot.data$saprowqclass <- as.factor(plot.data$saprowqclass)
plot.data$ReachID <- gsub("CSCF", "", plot.data$ReachID)
p <- ggplot(plot.data, aes(x = tempmaxC, y = currentms, 
                           size = orgmicropollTU, shape = saprowqclass,
                           color = Watershed)) #, label = ReachID))
p <- p + geom_point(alpha = 0.6)
if(ICE){
  p <- p + geom_point(data = data.base.ice, shape = 13, color = "#7CAE00", size = 5)
}
# if(dim(data.env.res)[1] < 10){ # write SiteID next to plot
#     p <- p + geom_text(size = 3, vjust = 0, nudge_y = 0.05)
# }
p <- p + theme_bw()
# p

# file.name <- "EnvData_MultiDim"
# print.pdf.plots(list.plots = list(p), 
#                 width = 10, height = 8, 
#                 dir.output = dir.outputs, info.file.name = name.run, file.name = file.name,
#                 png = F)

q <- ggplot(plot.data, aes(x = tempmaxC, y = currentms, 
                           color = orgmicropollTU, shape = saprowqclass)) #, label = ReachID))
q <- q + geom_point(alpha = 0.6)
if(ICE){
  q <- q + geom_point(data = data.base.ice, shape = 13, color = "#7CAE00", size = 5)
}
q <- q + scale_color_continuous(trans='reverse')
q <- q + facet_wrap(. ~ Watershed)
q <- q + theme_bw()
q <- q + theme(strip.background = element_rect(fill = "white"))
# q

file.name <- "EnvData_MultiDim"
print.pdf.plots(list.plots = list(p,q), 
                width = 10, height = 8, 
                dir.output = dir.outputs, info.file.name = name.run, file.name = file.name,
                png = F)

# Additional results ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(return.res.add){
    list.factors.type <- list(
        "Biomass" = c("Biomass"),
        "Rates" = c("r_basal", "r_resp", "r_prod", "r_death", "r_cons_tot", "r_miner"),
        "DeathFactors" = c("fsapro", "forgmicropoll"),
        "AlgaeFactors" = c("flimI", "flimP", "flimN",
                           "fshade", "flimnutrients", "fselfshade"),
        "InhibFactors" = c("fmicrohab", "fcurrent",
                           "ftempmax",  "fselfinh", "ffoodlim"),
        "FoodFactors" = c("sum_food", "sum_food_pref"))
    fact.types <- names(list.factors.type)
    # mycolors <- sample(scales::hue_pal()(20), length(unlist(list.factors.type)), replace = T)
    # mycolors <- rainbow(length(unlist(list.factors.type)))
    mycolors <- c(
      "Biomass" = "darkgreen",
      "r_basal" = "#ff7f0e",
      "r_resp" = "#2ca02c",
      "r_prod" = "#d62728",
      "r_death" = "#9467bd",
      "r_cons_tot" = "#8c564b",
      "r_miner" = "#e377c2",
      "fsapro" = "#7f7f7f",
      "forgmicropoll" = "#bcbd22",
      "flimI" = "#17becf",
      "flimP" = "#1f77b4",
      "flimN" = "#ff7f0e",
      "fshade" = "#2ca02c",
      "flimnutrients" = "#d62728",
      "fselfshade" = "#9467bd",
      "fmicrohab" = "#8c564b",
      "fcurrent" = "#e377c2",
      "ftempmax" = "#7f7f7f",
      "fselfinh" = "#bcbd22",
      "ffoodlim" = "#17becf",
      "sum_food" = "#1f77b4",
      "sum_food_pref" = "#ff7f0e"
    )

    scales::show_col(mycolors)
    
    processed.add.res <- get.plot.data.add.res(catch.results = list.results[[1]],
                                             list.factors.type = list.factors.type)

    # # further analyses
    # test.add.res <- list.results$Aare$res.add
    # colnames.add.res <- colnames(test.add.res)
    # add.res.fpom <- test.add.res[,which(grepl("FPOM", colnames.add.res))]
    # colnames(add.res.fpom)
    # write.csv(add.res.fpom, file = paste0(dir.outputs, name.run, "add.res.FPOM.csv"))
    # system.def.test <- streambugs.get.sys.def(y.names=list.results$Aare$y.names$y.names, par=list.results$Aare$par.fixc)
    # sink(file = paste0(dir.outputs, name.run, "sys.def.txt"))
    # print(system.def.test)
    # sink()
    
    list.df.res.add.taxa <- processed.add.res$list.df.add.res.taxa
    plot.data.all.taxa.sites <- processed.add.res$plot.data
    
    # plot all rates and limitation factors
    temp.vect.taxa <- unique(plot.data.all.taxa.sites$Taxon)
    # temp.vect.sites <- unique(plot.data$Site)
    list.plots <- list()
    for (taxon in temp.vect.taxa) {
      # taxon <- temp.vect.taxa[6]
        plot.data <- plot.data.all.taxa.sites %>%
            filter(Taxon == taxon) # %>%
        # filter(Site != 107)
        p <- ggplot(plot.data, aes(x = time, y = Value, color = Factor))
        p <- p + geom_line(size=0.7, alpha = 0.7)
        p <- p + facet_grid(FactorType ~ Site, scales = "free")
        if(!grepl("Algae", taxon)){ p <- p + scale_color_manual(values = mycolors)}
        p <- p + theme_bw()
        p <- p + labs(title = taxon)
        # ggplotly(p)
        
        list.plots[[taxon]] <- p
    }
    file.name <- paste0("AddRes_", length(temp.vect.taxa), "Taxa")
    # paste0(name.run, file.name)
    print.pdf.plots(list.plots = list.plots, width = 23, height = 8, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
                    png = F)
}

# Merge results in dataframes ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# fix parameters used to calculate steady state
ratio.tail <- 1/10 # which fraction of time vector should we take (at the end) to calculate the steady state by making the average
threshold.orig <- 0.001 # standard deviation tolerance of the distribution of the "averages" steady state, to be considered at equilibrium
threshold.min <- 6.162e-05*0.01 # 1% change in individual of the smallest invertebrate
threshold.max <- 3.600e-02*0.01 # 1% change in individual of the biggest invertebrate
threshold.med <- 5.350e-04*0.01 # 1% change in individual of the median invertebrate
threshold.mean <- 2.852e-03*0.01
threshold.test <- 3.947e-03*0.01 # 1% change in individual of the mean of the size of all invertebrates
threshold.select <- threshold.test

# analyze results for 1 sites, 20 steps, differences in steady state parameters (tail and threshold)

# options()$nwarnings  ## check current state
# # [1] 50
# 
# oo <- options(nwarnings=10000)  ## set, store old options
# 
# options()$nwarnings  ## check
# # [1] 10000
# 
# options(oo)  ## restore old options
# 
# extract results
# n <- 2
# catch <- names(list.results)[1]
# temp.res <- as.matrix(list.results[[n]][["res"]])
# temp.parfix <- list.results[[n]][["par.fixc"]]
# temp.y.names <- list.results[[n]][["y.names"]]
# temp.inv <- list.results[[n]][["Invertebrates"]]
# temp.sites <- temp.y.names$reaches
# summary(mass.all.inv)
# 
# # analysis steady state and warnings
# res.ss <- calc.res.steadystate.update(res=temp.res, times=tout, par=temp.parfix, y.names=temp.y.names,
#                                            ratio.tail = ratio.tail, threshold = threshold.select, warn.pom = F)
# D.mean.mod <- res.ss$D.mean.mod
# print(res.ss$Warnings)
# 
# # print warnings steady state
# sink(file = paste0(dir.outputs, name.run,
#                    "warnings_calcSS_",
#                    ratio.tail, "tail_",
#                    catch, "_",
#                    threshold.select, "thresh",
#                    ".txt"))
# cat("Ratio tail:", ratio.tail,
#     "\nThreshold selected:", threshold.select,
#     "\nNumber of warnings:", length(res.ss$Warnings),
#     "\nWarnings steady state:\n")
# print(res.ss$Warnings)
# sink()
# 
# # ind.war <- grepl("Baetisrhodani", res.ss$Warnings)
# # res.ss$Warnings[ind.war]
# 
# # plot time series of abundance for different analysis
# rind.time.tail <- c(ceiling(length(tout)*(1-ratio.tail)):length(tout))
# rind.time.select <- rind.time.tail # tout[c(1:100)]
# # colfunc <- colorRampPalette(c("darkturquoise","gold1", "red1"))
# # col.lines <- c(colfunc(no.steps.ice))
# # scales::show_col(col.lines)
# list.plots <- list()
# for (taxon in temp.inv) {
#   print(taxon)
#   plot.data <- data.frame()
#   for(site in temp.sites){
#     # site <- temp.sites[1]
#     temp.plot.data <- as.data.frame(temp.res) %>%
#       slice(rind.time.select) %>%
#       select(time, which(grepl(taxon, colnames(temp.res)) & grepl(site, colnames(temp.res)))) %>%
#       rename(biomass = 2) %>%
#       mutate(site = site)
#     plot.data <- bind_rows(plot.data, temp.plot.data)
#   }
#   # taxon <- temp.inv[1]
#   # mass.taxon <- mass.all.inv[which(grepl(taxon, names(mass.all.inv)))]
#   # temperature.range <- round(data.env.res$tempmaxC, digits = 2)
#   # plot.data <- as.data.frame(temp.res[rind.time.tail, which(grepl(taxon, colnames(temp.res)))])
#   # plot.data <- as.data.frame(temp.res[which(temp.res[,"time"] %in% rind.time.select), which(grepl(taxon, colnames(temp.res)))])
# 
#   # colnames(plot.data) <- temperature.range # 1:20 # paste0(select.sites, "_", 1:20)
#   # plot.data$Time <- rind.time.tail
#   # plot.data$Time <- rind.time.select
# 
#   # plot.data <- gather(plot.data, key = Site, value = Biomass, -Time)
#   plot.data$abundance <- plot.data$biomass/mass.all.inv[taxon]
#   # plot.data$Site <- as.factor(plot.data$Site)
#   # plot.data$Site <- factor(plot.data$Site, levels = temperature.range)
#   # plot.data$Site <- factor(plot.data$Site, levels= temp)
# 
#   p <- ggplot(plot.data, aes(x = time, y = abundance, color = site), size = 3)
#   p <- p + geom_line()
#   # p <- p + scale_color_manual(values = col.lines)
#   p <- p + labs(title = taxon) #,
#                 # color = "Temperature")
#   # p
#   list.plots[[taxon]] <- p
# }
# # file.name <- paste0("TimeSeriesAbundance_", ratio.tail, "TimeTail")
# 
# file.name <- paste0("TimeSeriesAbundance_", ratio.tail, "TimeTail_", catch)
# print.pdf.plots(list.plots = list.plots, width = 8, height = 6, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name,
#                 png = F)


# # analyze change in results with different tout
# temp.res <- as.matrix(list.results[[1]][["res"]])
# site <- vect.sites[1]
# length(tout)
# print(name.run)
# res1.site <- as.data.frame(temp.res[,c(1, which(grepl(site, colnames(temp.res))))]) # 7301 tout
# res1.site$time <- round(res1.site$time, digits = 4)
# res2.site <- read.table(paste0(dir.outputs, "streambugs_results_", site,
#                               "_RheinabBS_5Sites_ICE_update.traits.hybrid_rk4Met_10Yea_3651Steps_.dat"),
#                         header = T, sep = "\t") # 731 tout
# res2.site$time <- round(res2.site$time, digits = 4)
# vect.time <- res2.site$time # choose the coarser time steps
# # vect.time.analysis <- vect.time[c((length(vect.time)/2):length(vect.time))]
# length(vect.time)
# vect.time.analysis <- vect.time[c(round(9*length(vect.time)/10):length(vect.time))]
# 
# list.plots <- list()
# for (taxon in vect.taxa) {
#   # taxon <- vect.taxa[1]
#   temp1 <- res1.site %>%
#     filter(time %in% vect.time.analysis) %>%
#     select(time, which(grepl(taxon,colnames(res1.site)))) %>%
#     mutate(method = "lsoda")
#   temp2 <- res2.site %>%
#     filter(time %in% vect.time.analysis) %>%
#     select(time, which(grepl(taxon,colnames(res2.site)))) %>%
#     mutate(method = "rk4")
#   plot.data <- bind_rows(temp1, temp2)
#   colnames(plot.data)[2] <- "biomass"
# 
#   p <- ggplot(plot.data, aes(x = time, y = biomass, color = method, shape = method), size = 3)
#   p <- p + geom_point()
#   # p <- p + scale_color_manual(values = col.lines)
#   p <- p + labs(title = taxon)
#   # p
#   list.plots[[taxon]] <- p
# 
# }
# 
# file.name <- paste0("EndTimeSeriesBiomass_","ComparMethods_rk4_lsoda")
# print.pdf.plots(list.plots = list.plots, width = 8, height = 6, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name,
#                 png = F)



# source("utilities.r")

print(vect.analysis.env.fact)
vect.analysis <- c("ReachID", "X", "Y", "MonitoringProgram", catch.variable, vect.analysis.env.fact)

# make dataframes (wide and long) with steady state, abundance and prob of occ
# df.generic <- data.env.res[,vect.analysis]
df.generic <- data.env.res[,vect.analysis]

df.generic[,vect.occ.inv] <- NA

vect.names.output <- c("SteadyState", "Abundance", "ThreshPresAbs", "ProbObs", "SamplePresAbs")
list.wide.df.all.res <- list()
for (output in vect.names.output) {
    list.wide.df.all.res[[output]] <- df.generic
}

# print warnings steady state
sink(file = paste0(dir.outputs, name.run, 
                   "warnings_calcSS_", 
                   ratio.tail, "tail_",
                   threshold.select, "thresh",
                   ".txt"))
cat("Ratio of the tail of the output time vector:", ratio.tail,
    "\nDistribution of invertebrates mass/size:", summary(mass.all.inv),
    "\nThreshold selected:", threshold.select)
    
start <- proc.time()

for (n in 1:length(list.results)) {
    # n <- 1
    temp.res <- as.matrix(list.results[[n]][["res"]])
    temp.parfix <- list.results[[n]][["par.fixc"]]
    temp.y.names <- list.results[[n]][["y.names"]]
    temp.inv <- list.results[[n]][["Invertebrates"]]
    temp.sites <- temp.y.names$reaches
    temp.res.ss <- calc.res.steadystate.update(res=temp.res, times=tout, par=temp.parfix, y.names=temp.y.names,
                                               ratio.tail = ratio.tail, threshold = threshold.select, warn.pom = F)
    D.mean.mod <- temp.res.ss$D.mean.mod
    cat("\nCatchment:", names(list.results)[n],
      "\nNumber of warnings:", length(temp.res.ss$Warnings),
        "\nWarnings steady state:\n")
    print(temp.res.ss$Warnings)
    
    for (site in temp.sites) {
        # site <- temp.sites[2]
        for (taxon in temp.inv) {
            # taxon <- temp.inv[2]
          occ.taxon <- paste0("Occurrence.", taxon)
          ind <- which(grepl(taxon, names(D.mean.mod)) 
                       & grepl(paste0("^", site, "_"), # "^" Asserts that we are at the start because many sites could have same pattern
                               names(D.mean.mod)))
          temp.ss <- D.mean.mod[ind]
          if(any(is.na(temp.ss))){ # print warning if NA in results
              cat("\nWARNING NA in results:\nTaxon:", taxon,
                  "\nSite:", site, "\n")
          }
          temp.m.taxon <- mass.all.inv[which(grepl(taxon, names(mass.all.inv)))]
          temp.abund <- temp.ss/temp.m.taxon
          if(par.adjust[["thresh.indiv"]][["Flag"]]){
            temp.thresh.occ <- ifelse(temp.abund > par.adjust[["thresh.indiv"]][["Value"]], 1, 0) # if there is more than 1 individuals it's present, otherwise absent
          } else {
            temp.thresh.occ <- ifelse(temp.abund > 0.5, 1, 0) # if there is more than 0.5 individuals it's present, otherwise absent
          }
          temp.prob.obs <- get.prob.obs(ss = temp.ss, M.taxon = temp.m.taxon)
          # print(temp.prob.obs)
          # rbinom(100,1,temp.prob.obs) # test rbinom
          temp.sample.occ <- rbinom(1,1,temp.prob.obs)
          
          for (output in vect.names.output) {
              # output <- vect.names.output[1]
              out <- ifelse(output == "SteadyState", temp.ss, 
                            ifelse(output == "Abundance", temp.abund, 
                                   ifelse(output == "ThreshPresAbs", temp.thresh.occ, 
                                          ifelse(output == "ProbObs", temp.prob.obs, 
                                                 ifelse(output == "SamplePresAbs", temp.sample.occ, 
                                                        NA)))))
              list.wide.df.all.res[[output]][which(list.wide.df.all.res[[output]][,"ReachID"] == site), occ.taxon] <- out
              
          }
        }
    }
}
duration <- proc.time()-start
cat("Time calculation steady state:", duration,"\n")
sink()

# recover observation data for selected sites
data.obs <- df.generic
data.inv.bdm$ReachID <- gsub("_", ".", data.inv.bdm$SiteId)
for (site in data.obs$ReachID) {
    # site <- df.generic$ReachID[1]
    temp.monit.prog <- data.obs[which(data.obs$ReachID == site), "MonitoringProgram"]
    # if(site %in% data.inv.bdm$ReachID){
        for(taxon in vect.inv){
          occ.taxon <- paste0("Occurrence.", taxon)
            if(temp.monit.prog == "BDM"){
                obs <- data.inv.bdm[which(data.inv.bdm$ReachID == site), which(grepl(taxon, colnames(data.inv.bdm)))]
            } else {
                cind.taxon <- which(grepl(taxon, colnames(data.env.res)))
                if(length(cind.taxon) > 0){ # if tayon exists in all monit. prog. dataset
                    obs <- data.env.res[which(data.env.res$ReachID == site), cind.taxon]
                } else {
                    obs <- NA
                }
            }
            if(length(obs) == 0){ obs <- NA}
            if(!is.na(obs) & obs == 0){
                data.obs[which(data.obs$ReachID == site), occ.taxon] <- "Absent"
            } else if(!is.na(obs) & obs == 1){
                data.obs[which(data.obs$ReachID == site), occ.taxon] <- "Present"
            }
        }
    # }
}

# add it to list with other outputs
vect.names.output <- c(vect.names.output, "Observation")
list.wide.df.all.res[["Observation"]] <- data.obs

# make list long dataframes
list.long.df.res.taxa <- list()
for (output in vect.names.output) {
    # output <- vect.names.output[1]
    temp.df.output <- list.wide.df.all.res[[output]]
    temp.long.df.output <- gather(temp.df.output, key = Taxon, value = Output, -all_of(vect.analysis))
    colnames(temp.long.df.output)[which(colnames(temp.long.df.output) == "Output")] <- output
    list.long.df.res.taxa[[output]] <- temp.long.df.output
}

# make long dataframe with all taxa (env. fact. and outputs stay wide)
long.df.res.taxa <- list.long.df.res.taxa %>% purrr::reduce(merge, by = c(vect.analysis, "Taxon"))
long.df.res.taxa$Taxon <- gsub("Occurrence.", "", long.df.res.taxa$Taxon)

# make long dataframe with all taxa and env. factors (outputs stay wide)
long.df.res.taxa.env.fact <- gather(long.df.res.taxa, 
                                key = Factor, value = Value, 
                                -c(which(!colnames(long.df.res.taxa) %in% vect.env.inputs)))

# make long dataframe with all taxa and outputs (env. fact. stay wide)
long.df.res.taxa.outputs <- gather(long.df.res.taxa, 
                                    key = Output, value = Result, 
                                    -c(which(!colnames(long.df.res.taxa) %in% vect.names.output)))

# make long dataframe with everything
long.df.all <- gather(long.df.res.taxa.env.fact,
                      key = Output, value = Result,
                      -c(which(!colnames(long.df.res.taxa.env.fact) %in% vect.names.output)))


# Write final datasets ####
for (output in vect.names.output) {
  # output <- vect.names.output[1]
  wide.df.output <- list.wide.df.all.res[[output]]
  file.name <- paste0(dir.outputs, name.run, "WideData_Results", output,".csv")
  write.table(wide.df.output, file = file.name, sep =";", row.names = F)
}

# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # analysis with another run
# 
# name.run.other <- "Ticino_20Sites_ICE_update.traits.hybrid_CorrectPar_100Yea_36501Steps_"
# list.res.wide.other <- list()
# long.df.all.other <- data.frame()
# for (output in vect.names.output) {
#   # output <- vect.names.output[1]
#   file.name <- paste0(dir.outputs, name.run.other, "WideData_Results", output,".csv")
#   wide.df.output <- read.csv(file.name, sep = ";")
#   list.res.wide.other[[paste0("wide.other_", output)]] <- wide.df.output
#   # long.df.output
#   temp.long.df.output1 <- gather(wide.df.output, key = Taxon, value = Result, -all_of(vect.analysis))
#   temp.long.df.output2 <- gather(temp.long.df.output1, key = Factor, value = Value,
#                                  -c(which(!colnames(temp.long.df.output1) %in% vect.env.inputs)))
#   temp.long.df.output2$Output <- output
#   # colnames(temp.long.df.output2)[which(colnames(temp.long.df.output2) == "Output")] <- output
#   long.df.all.other <- bind_rows(long.df.all.other, temp.long.df.output2)
# }
# 
# long.df.all.other$Taxon <- gsub("Occurrence.", "", long.df.all.other$Taxon)

# make dataframe with prevalence and taxonomy of results
print(vect.names.output[c(3,5)])
# select two different of the Presence/Absence outputs (threshold or sample) 
# to calculate prevalence in results
for (prev.output in vect.names.output[c(3,5)]) {
  print(prev.output)
  df.prev.results <- data.frame("Taxon" = vect.inv)
  df.prev.results$Occurrence.taxa <- vect.occ.inv
  df.prev.results$Prevalence <- NA
  df.prev.results$Prevalence.NoDisp <- NA
  df.prev.results$Taxonomic.level <- NA
  df.prev.results$Missing.values <- 0
  
  for (occ.taxon in df.prev.results$Occurrence.taxa) {
    # occ.taxon <- df.prev.results$Occurrence.taxa[3]
    rind <- which(df.prev.results$Occurrence.taxa == occ.taxon)
    taxon <- df.prev.results$Taxon[rind]
    nb.sites <- dim(list.wide.df.all.res[[prev.output]])[1]
    vect.pres.abs.tax <- list.wide.df.all.res[[prev.output]][,occ.taxon]
    nb.pres <- sum(na.omit(vect.pres.abs.tax))
    nb.abs <- sum(na.omit(vect.pres.abs.tax) == 0)
    nb.na <- sum(is.na(vect.pres.abs.tax))
    prev <- nb.pres/nb.sites
    if(nb.abs != 0){
      prev.nodisp <- nb.pres/(nb.abs + nb.pres) # sum abs and pres without NA
    } else { # if there are no abs
      prev.nodisp <- 1
    }
    tax.level <- data.prevalence[which(data.prevalence$Taxon == taxon), "Taxonomic.level"]
    df.prev.results[rind, "Prevalence"] <- prev
    df.prev.results[rind, "Prevalence.NoDisp"] <- prev.nodisp
    df.prev.results[rind, "Taxonomic.level"] <- tax.level
  }
  
  # write csv df prevalence
  file.name <- paste0(dir.outputs, name.run, "PrevalenceTaxonomy_", prev.output,".csv")
  write.table(df.prev.results, file = file.name, sep =";", row.names = F)
}

### ICE streambugs ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(ICE){
  
  # data.base.ice <- read.csv(paste0(dir.outputs, "WideDataInputs_ICE_50Sites.csv"), sep = ";")
  # data.results.ice <- read.csv(paste0(dir.outputs,
  # "10Catch_50Sites_ICE_30Steps_PolyInterp_20runC_update.traits.hybrid_10Yea_3651Steps_WideData_ResultsProbObs.csv"), sep = ";")
  # no.sites.ice <- 2
  # no.steps.ice <- 30
  
  data.results.ice <- list.wide.df.all.res$ProbObs
  select.env.fact <- vect.dir.env.fact[2]
  print(select.env.fact)
  vect.reaches.ice <- data.base.ice$ReachID
  
  # catch <- "Aare"
  # 
  # 
  # 
  # data.test <- data.results.ice %>%
  #   filter(Watershed == catch)
  # # filter(Watershed == "Ticino")
  # # filter(Watershed == "Rhone")
  # 
  # vect.reaches.test <- data.base.ice %>%
  #   filter(Watershed == catch)
  # vect.reaches.test <- vect.reaches.test$ReachID
  # 
  # list.problems.na <- list()
  # list.problems.optima <- list()
  # # test if some outputs are not monotously increasing/decreasing (multiple optima)
  # for(reach in vect.reaches.test){
  #   # reach <- vect.reaches.ice[1]
  #   # reach <- "SynthPoint158Ti"
  #   # print(reach)
  #   # watershed <- data.base.ice[which(data.base.ice$ReachID == reach), "Watershed"]
  #   # for(taxon in selected.taxa.analysis){
  #     taxon <- selected.taxa.analysis[8]
  #     temp.results <- data.test[which(grepl(reach, data.test$ReachID)), c("ReachID", select.env.fact, taxon)]
  #     vect.res <- temp.results[,taxon]
  #     ind.max <- which(ggpmisc::find_peaks(vect.res,  span = 3))
  #     ind.na <- which(is.na(vect.res))
  #     name.problem <- paste0(reach,"_",taxon)
  #     if(length(ind.na) > 0){
  #       list.problems.na[[name.problem]] <- c(problem = "NA",
  #                                          indices = ind.na,
  #                                          reach = reach,
  #                                          catch = watershed,
  #                                          taxon = taxon)
  # 
  #     } else if(length(ind.max) > 1){
  #       cat("Problem for", name.problem, "\n")
  #       plot(1:length(vect.res), vect.res, main = name.problem)
  #       points(ind.max, vect.res[ind.max], col = "red")
  #       list.problems.optima[[name.problem]] <- c(problem = "Optima",
  #                                          indices = ind.max,
  #                                          reach = reach,
  #                                          catch = watershed,
  #                                          taxon = taxon)
  #       }
  #   # }
  # }
  # 
  # rind.na <- which(is.na(data.results.ice$Occurrence.Gammaridae))
  # sites.na <- data.results.ice[rind.na, c("ReachID", env.factor[1:8])]
  
  # temp.vect.taxa <- names(selected.taxa.analysis)[names(selected.taxa.analysis) %in% vect.taxa]
  temp.vect.taxa <- vect.inv
  temp.occ.taxa <- paste0("Occurrence.", temp.vect.taxa)
  
  list.plots <- list()
  
  for(taxon in temp.occ.taxa){
    
    # taxon <- temp.occ.taxa[5] # beatisalpinus
    # taxon <- temp.occ.taxa[14] # gammaridae
    
    temp.reach.id <- stringr::str_sub(data.results.ice$ReachID[seq(1, length(data.results.ice$ReachID), no.steps.ice)], end =-2)

    pred.ice.streamb <- data.results.ice[,c("ReachID","Watershed", select.env.fact, taxon)] %>%
      mutate(observation_number = rep(1:no.sites.ice, each = no.steps.ice),
             reach_id = rep(temp.reach.id, each = no.steps.ice),
             column_label = 1,
             model = "Streambugs") %>%
      rename(pred = taxon)
    # pred.ice.streamb <- filter(pred.ice.streamb, Watershed == "Rhein", observation_number == 28)
    
    env.factor.sampled <- data.frame(variable = data.base.ice[,select.env.fact])
    
    pred.ice.mean <- pred.ice.streamb %>%
      group_by_at(select.env.fact[1]) %>%
      summarize(avg = median(na.omit(pred)))
    # pred.ice.mean.bounds  <- pred.ice.mean %>% 
    #   # group_by(model) %>%
    #   summarize(x.mean=max(select.env.fact),
    #             y.mean.min=min(avg),
    #             y.mean.max=max(avg))
    

    plot.data <- pred.ice.streamb
    # plot.data$model <- factor(plot.data$model, levels=c("GLM", "GAM", "RF"))
    
    plot.data.trait <- long.df.trait %>%
      filter(Taxon %in% gsub("Occurrence.", "", taxon),
             Factor %in% select.env.fact) %>%
      rename(PrefTrans = "ExpTrans_-10")

    p <- ggplot(data=plot.data) 
    p <- p + geom_line(aes(x=.data[[select.env.fact]],
                    y=pred,
                    group=observation_number,
                    color=as.character(reach_id)))#,
                # show.legend = FALSE) 
    # p <- p + geom_line(data=pred.ice.mean,
    #             aes(x=.data[[select.env.fact]], 
    #                 y=avg),size=1.5)
    p <- p + geom_rug(data = env.factor.sampled,
               aes(x=variable), 
               color="grey20",
               alpha=0.7,
               inherit.aes=F) 
    p <- p + geom_point(data = plot.data.trait, aes(x = Value, y = PrefTrans), shape = 4, color = "#7CAE00", size = 3)
    # p <- p + geom_segment(data=pred.ice.mean.bounds,
    #                inherit.aes = FALSE,
    #                lineend="round",
    #                linejoin="round",
    #                aes(x=x.mean,
    #                    y=y.mean.min,
    #                    xend=x.mean,
    #                    yend=y.mean.max),
    #                arrow=arrow(length = unit(0.3, "cm"),
    #                            ends = "both"))
    # p <- p + facet_wrap(~Watershed)
    p <- p + ylim(0,1)
    p <- p + theme_bw()
    # p <- p + theme(strip.background = element_rect(fill = "white"),
    #               legend.title = element_text(size=24),
    #               legend.text = element_text(size=20))
    p <- p + labs(title = taxon,
                   x = "Temperature",
                   y = "Predicted probability of occurrence",
                   color = "Index site")
    # p
    # ggplotly(p)
    list.plots[[taxon]] <- p
  }
  
  file.name <- paste0("ice_", ifelse(no.sites.ice == 1, select.sites, no.sites.ice), "_Sites_", no.steps.ice, "Steps_", select.env.fact)
  # file.name <- paste0(file.name, "_perCatch")
  # paste0(name.run, file.name)
  print.pdf.plots(list.plots = list.plots, width = 10, height = 10, 
                  dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
                  png = F)
  
  
  # # analysis baetis rhodani 
  # m.rhodani <- mass.all.inv[4]
  # vect.ss.rhodani <- list.wide.df.all.res$SteadyState$Occurrence.Baetisrhodani
  # plot(data.ICE$tempmaxC, vect.ss.rhodani)
  # vect.abund.rhodani <- vect.ss.rhodani / m.rhodani
  # plot(data.ICE$tempmaxC, vect.abund.rhodani)
  # vect.prob.rhodani <- sapply(vect.ss.rhodani, FUN = get.prob.obs, m.rhodani)
  # plot(data.ICE$tempmaxC, vect.prob.rhodani)

}

### hist pres/abs all taxa ####

# take back prevalence information from data
temp.df.prev <- data.prevalence[which(data.prevalence$Taxon %in% vect.inv), c("Taxon", "Prevalence")]
temp.df.prev$Taxon <- paste(sub("Occurrence.","", temp.df.prev$Taxon))
temp.df.prev <- arrange(temp.df.prev, desc(Prevalence))
temp.df.prev$TaxonPrev <- paste(temp.df.prev$Taxon, temp.df.prev$Prevalence, sep = "_")
temp.df.prev$ScaledPrevalence <- n.sites * temp.df.prev$Prevalence

plot.data <- long.df.res.taxa.outputs
plot.data <- filter(plot.data, grepl("PresAbs", plot.data$Output))
plot.data$Result <- ifelse(plot.data$Result == 1, "Present", ifelse(plot.data$Result == 0, "Absent", NA))
plot.data$Result <- factor(plot.data$Result, levels=c("Present", NA, "Absent"))

p <- ggplot(plot.data, aes(x = Taxon, fill = Result))
p <- p + geom_bar(stat="count", position='stack')
p <- p + facet_wrap(~ Output, ncol = 1)
p <- p + scale_fill_manual(values=vect.col.pres.abs)
p <- p + theme_bw()
p <- p + scale_x_discrete(limits = temp.df.prev$Taxon)
p <- p + theme(axis.text.x = element_text(angle=90))
# p

file.name <- "GeomBar_PresAbs_allTaxa"
print.pdf.plots(list.plots = list(p), width = 22, height = 12, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
                png = F)


### response to env fact ####

# # merge two different runs to compare them
# temp.long.df <- long.df.all
# temp.long.df$Run <- "50years"
# temp.long.df2 <- long.df.all.other
# temp.long.df2$Run <- "100years"
# temp.plot.data <- bind_rows(temp.long.df, temp.long.df2)

list.plots <- list()
for(taxon in vect.inv){
    # taxon <- vect.inv[4]
    print(taxon)
    plot.data <- long.df.all %>% # to comment if comparing different runs
      # plot.data <- temp.plot.data %>% # to compare different runs
      
      filter(Taxon %in% taxon) %>%
      filter(Factor %in% vect.dir.env.fact) %>%
      filter(Output %in% vect.names.output[c(2,4,5)])
    plot.data$Result <- as.numeric(plot.data$Result)
    plot.data$Dispersal <- ifelse(is.na(plot.data$Result), "Yes", "No")
    plot.data$Result[is.na(plot.data$Result)] <- 0
    plot.data$Factor <- factor(plot.data$Factor, levels=vect.dir.env.fact)
    # summary(plot.data$Result)
    
    # plot.data.trait <- filter(long.df.trait, Taxon %in% taxon)
    plot.data.trait <- long.df.trait %>%
      filter(Taxon %in% taxon) %>%
      rename(PrefTrans = "ExpTrans_-10",
             PrefOrig  = "Pref") %>%
      gather( key = "Transformed", value = "Pref", c("PrefTrans", "PrefOrig"))
    plot.data.trait$Factor <- factor(plot.data.trait$Factor, levels=vect.dir.env.fact)
    plot.data.trait$Output <- "ProbObs"
    
    
    p <- ggplot(plot.data, aes(x = Value, y = Result, shape = Dispersal), color = "orange")
    # p <- ggplot(plot.data, aes(x = Value, y = Result, color = Run)) # to compare different runs
    
    p <- p + geom_point(alpha = 0.5)
    p <- p + geom_point(data = plot.data.trait, aes(x = Value, y = Pref, color = Transformed), shape = 4, size = 1, alpha = 0.6)
    p <- p + facet_grid(Output ~ Factor, scales = "free")
    # p <- p + scale_y_continuous(limits = c(0,1))
    p <- p + labs(title = taxon)
    p <- p + theme_bw()
    p <- p + theme(strip.background = element_rect(fill = "white"))
    # p <- p + scale_color_manual(values=c("Yes" = "grey60", "No" = "orange")) # to comment if comparing different runs
    p <- p + scale_color_manual(values=c("PrefOrig" = "#619CFF", "PrefTrans" = "#F8766D")) # to comment if comparing different runs
    # ggplotly(p)
    # p
    list.plots[[taxon]] <- p
}
# list.plots$Ceratopogonidae

# file.name <- paste0("ResponseDirEnvFact_allTaxa_Compar50_100years") # to compare different runs
file.name <- paste0("ResponseDirEnvFact_allTaxa") # to compare different runs
print.pdf.plots(list.plots = list.plots, width = 8, height = 6, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
                png = F)

# extract selected taxa for analysis and plot only these ones
names.selected.taxa <- names(selected.taxa.analysis)
list.plots.selected.taxa <- Filter(Negate(is.null), list.plots[names.selected.taxa])

file.name <- paste0("ResponseDirEnvFact_SelectedTaxa")
print.pdf.plots(list.plots = list.plots.selected.taxa, width = 8, height = 6, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
                png = F)

# OTHER PLOTS ####

## histograms ###

# ## hist biomass per taxon ####
# print(vect.names.output)
# for (output in vect.names.output[1:2]) {
#     file.name <- paste0(dir.outputs, name.run, output, "_Hist", ".pdf")
#     pdf(file.name,width=8.25,height=11.75)
#     par(mfrow=c(5,4))
#     for (taxon in vect.inv) {
#         # taxon <- Invertebrates[1]
#         taxon.ss <- list.wide.df.all.res[[output]][,taxon]
#         hist(taxon.ss, main = taxon, 
#              xlab = output)
#     }
#     dev.off()
# }

# ### hist pres/abs per taxon ####
# list.plots <- list()
# 
# for (taxon in vect.inv) {
#     taxon <- vect.inv[1]
#     print(taxon)
#     plot.data <- long.df.res.taxa[,c(vect.analysis, "Taxon", "ThreshPresAbs", "Observation")]
#     plot.data$ThreshPresAbs <- ifelse(plot.data$ThreshPresAbs == 1, "Present", ifelse(plot.data$ThreshPresAbs == 0, "Absent", NA))
#     plot.data <- filter(plot.data, Taxon %in% taxon)
#     plot.data <- gather(plot.data, key = Output, value = Value, -all_of(c(vect.analysis, "Taxon")))
#     plot.data$Value <- factor(plot.data$Value, levels=c("Present", NA, "Absent"))
#     prev.taxon <- data.prevalence[which(data.prevalence$Taxon == taxon), "Prevalence"]
#     
#     p <- ggplot(plot.data, aes(x = Output, fill = Value))
#     p <- p + geom_bar(stat="count", position='stack')
#     p <- p + theme_bw()
#     p <- p + scale_fill_manual(values=c(vect.col.pres.abs, "NA" = "grey"))
#     p <- p + labs(title = taxon,
#                   subtitle = paste("Prevalence:", prev.taxon))
#     list.plots[[taxon]] <- p
# }
# list.plots$Baetisrhodani
# 
# file.name <- paste0("GeomBar_ThreshPresAbsvsObs_perTaxon")
# print.pdf.plots(list.plots = list.plots, width = 6, height = 5, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name,
#                 png = F)


# plot.data <- long.df.res.taxa[,c(vect.analysis, "Taxon", "ThreshPresAbs", "Observation")]
# plot.data$ThreshPresAbs <- ifelse(plot.data$ThreshPresAbs == 1, "Present", ifelse(plot.data$ThreshPresAbs == 0, "Absent", NA))
# # plot.data <- filter(plot.data, Taxon %in% taxon)
# plot.data <- gather(plot.data, key = Output, value = Value, -all_of(c(vect.analysis, "Taxon")))
# p <- ggplot(plot.data, aes(x = Taxon, fill = Value))
# p <- p + geom_bar(stat="count",position='dodge')
# p <- p + theme(axis.text.x = element_text(angle=90))
# p

## abund vs env fact ####


# ### resp abund to one env fact ####
# 
# list.plots <- list()
# for(taxon in vect.inv){
#     print(taxon)
#     plot.data <- long.df.res.taxa.env.fact %>%
#         filter(Taxon %in% taxon) %>%
#         filter(Factor %in% vect.analysis.env.fact)
#     
#     p <- ggplot(plot.data, aes(x = Value, y = Abundance))
#                                # , color = Observation))
#     p <- p + geom_point()
#     p <- p + facet_wrap(. ~ Factor, scales = "free_x")
#     # p <- p + scale_y_continuous(limits = c(0,1))
#     p <- p + labs(title = taxon)
#     p <- p + theme_bw()
#     # p <- p + scale_color_manual(values=vect.col.pres.abs)
#     # p
#     list.plots[[taxon]] <- p
# }
# # list.plots$Ceratopogonidae
# 
# file.name <- paste0("AbundanceVsEnvFact")
# print.pdf.plots(list.plots = list.plots, width = 12, height = 8, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
#                 png = F)


# ### resp prob occ to indirect env fact ####
# 
# list.plots <- list()
# for(taxon in vect.inv){
#     # taxon <- vect.inv[1]
#     print(taxon)
#     plot.data <- long.df.all %>%
#       filter(Taxon %in% taxon) %>%
#       filter(Factor %in% vect.indir.env.fact) %>%
#       filter(Output %in% vect.names.output[3:5])
#     plot.data$Result <- as.numeric(plot.data$Result)
#     plot.data.trait <- filter(long.df.trait, Taxon %in% taxon)
#     
#     p <- ggplot(plot.data, aes(x = Value, y = Result))
#     p <- p + geom_point(alpha = 0.7, color = "gray30")
#     p <- p + facet_grid(Output ~ Factor, scales = "free_x")
#     p <- p + scale_y_continuous(limits = c(0,1))
#     p <- p + labs(title = taxon)
#     p <- p + theme_bw()
#     p <- p + scale_color_manual(values=vect.col.pres.abs)
#     # p
#     list.plots[[taxon]] <- p
# }
# list.plots$Ceratopogonidae
# 
# file.name <- paste0("ThreshPresAbsVsIndirEnvFact")
# print.pdf.plots(list.plots = list.plots, width = 8, height = 6, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
#                 png = F)


# ### resp to two env fact
# print(vect.dir.env.fact)
# list.plots <- list()
# for (taxon in vect.inv) {
#     # taxon <- vect.inv[1]
#     # plot.data <- long.df.all %>%
#     #     filter(Taxon %in% taxon) %>%
#     #     filter(Factor %in% vect.dir.env.fact)
#     # plot.data$LimFact <- ifelse(plot.data$Factor %in% vect.dir.env.fact[1:2], "Growth", "Death")
#     # plot.data$Fact1 <- plot.data$Factor
#     # plot.data$Val1 <- plot.data$Value
#     # plot.data$Fact2 <- plot.data$Factor
#     # plot.data$Val2 <- plot.data$Value
#     # 
#     # p <- ggplot(plot.data, aes(x = Val1, y = Val2, color = Observation))
#     # p <- p + geom_point()
#     # p <- p + facet_grid(Fact1 ~ Fact2, scales = "free")
#     # # p <- p + scale_y_continuous(limits = c(0,1))
#     # # p <- p + labs(title = taxon)
#     # p <- p + theme_bw()
#     # p <- p + scale_color_manual(values=vect.col.pres.abs)
#     # p
#     
#     temp.list.plots <- list()
#     for (n in 1:nrow(df.comb.dir.env.fact)) {
#         plot.data <- long.df.res.taxa %>%
#             filter(Taxon %in% taxon)
#         # n <- 1
#         comb.env.fact <- df.comb.dir.env.fact[n,]
#         p <- ggplot(plot.data, aes_string(x = comb.env.fact[1], y = comb.env.fact[2], 
#                                    size = "Abundance", color = "Observation"))
#         p <- p + geom_point()
#         p <- p + scale_color_manual(values = vect.col.pres.abs)
#         p <- p + theme_bw()
#         p <- p + theme(legend.position="none")
#         
#         # p
#         temp.list.plots[[n]] <- p
#     }
#     # temp.list.plots[[2]]
#     q <- gridExtra::grid.arrange(grobs = temp.list.plots, ncol = 2, top = taxon)
#     
#     list.plots[[taxon]] <- q
# }
# 
# file.name <- paste0("AbundVsEnvFact_TwoFact")
# print.pdf.plots(list.plots = list.plots, width = 12, height = 8, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
#                 png = F)


# ### resp env fact multidim ####
# # taxon <- vect.inv[6]
# 
# list.plots <- list()
# for(taxon in vect.inv){
#     plot.data <- filter(long.df.res.taxa, Taxon %in% taxon)
#     plot.data$saprowqclass <- as.factor(plot.data$saprowqclass)
#     colnames(plot.data)
#     p <- ggplot(plot.data, aes(x = tempmaxC, y = currentms, 
#                                size = Abundance, color = orgmicropollTU, shape = saprowqclass))
#     p <- p + geom_point(alpha = 0.7)
#     p <- p + labs(title = taxon)
#     p <- p + theme_bw()
#     list.plots[[taxon]] <- p
# }
# file.name <- paste0("AbundanceVsEnvFact_MultiDimensional")
# print.pdf.plots(list.plots = list.plots, width = 10, height = 8, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
#                 png = F)


# ## prob occ vs obs Swiss map ####
# 
# map.inputs    <- map.inputs(directory = paste0(dir.inputs,"swiss.map.gdb"))
# 
# list.plots <- list()
# for(taxon in vect.inv){
#     # plot predicted probabilities on map
#     # taxon <- vect.survivors[1]
#     print(taxon)
#     plot.data <- long.df.all %>%
#         filter(Taxon %in% taxon )
#     prev.taxon <- data.prevalence[which(data.prevalence$Taxon == taxon), "Prevalence"]
#     
#     # Map geometries
#     g <- ggplot()
#     g <- g + geom_sf(data = map.inputs$ch, fill=NA, color="black")
#     g <- g + geom_sf(data = map.inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
#     g <- g + geom_sf(data = map.inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)
#     g <- g + geom_point(data = plot.data, aes(X, Y, size = ProbObs, 
#                                               color = Observation), 
#                         alpha = 0.7)
#     # g <- g + facet_wrap(~ Model, strip.position="top", ncol = 2)
#     
#     # Configure themes and labels
#     g <- g + theme_void()
#     g <- g + theme(plot.title = element_text(size = 14, vjust = 10),
#                    panel.grid.major = element_line(colour="transparent"),
#                    plot.margin = unit(c(1,1,1,1), "lines") # ,
#                    # legend.title = element_text(size=24),
#                    # legend.text = element_text(size=20),
#                    # strip.text = element_text(size=20)
#     )
#     
#     g <- g + labs(title = taxon,
#                   subtitle = paste0("Prevalence: ", prev.taxon),
#                   x = "",
#                   y = "",
#                   size = "Probability of\noccurrence",
#                   color = "Observation")
#     
#     # Configure legends and scales
#     g <- g + guides(size = guide_legend(override.aes = list(color="black", stroke=0), order=1))
#     # g <- g + scale_size(range = c(0.3,1.5))
#     g <- g + scale_radius(limits = c(0,1), breaks = seq(0, 1, 0.25), range = c(0.3, 2.5))
#     
#     # g <- g + scale_size_area(max_size = 2)
#     g <- g + scale_color_manual(values=vect.col.pres.abs)
#     g <- g + scale_shape_identity() # Plot the shape according to the data
#     g
#     
#     list.plots[[taxon]] <- g
# }
# Tsarabe7!

# file.name <- paste0("Maps_ProbObs_vs_Obs")
# print.pdf.plots(list.plots = list.plots, width = 10, height = 8, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
#                 png = F)
# 
# 
