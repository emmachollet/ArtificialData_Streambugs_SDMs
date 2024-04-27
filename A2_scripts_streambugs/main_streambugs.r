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
if ( !require("tidyverse") ) { install.packages("tidyverse"); library("tidyverse") }               # to sort, join, merge data
if ( !require("ggplot2") ) { install.packages("ggplot2"); library("ggplot2") }         # to do nice plots
if ( !require("plotly") ) { install.packages("plotly"); library("plotly") }            # to do nice "interactive" plots
if ( !require("sf") ) { install.packages("sf"); library("sf") }                        # to read layers for plotting maps

## Directory and file definitions ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dir.inputs      <- "../A1_inputs_streambugs/"
dir.outputs     <- "../A3_outputs_streambugs/"
dir.utilities   <- "../00_utilities/"

file.midat.all    <- "All_2339samples_134taxa_13envfact_ProcessedDataset.csv" # output of: PreprocessSwissInvertebrateEnvironmentalModellingData
file.midat.bdm    <- "BDM_950samples_294taxa_13envfact_ProcessedDataset.csv" # output of: PreprocessSwissInvertebrateEnvironmentalModellingData
file.taxonomy     <- "invertebrates_taxonomy_2023-08-03.dat"  # output of: PreprocessSwissInvertebrateEnvironmentalModellingData
file.par.update   <- "T1d_BDM_CH_Tvsappestsubst_maxpost_trait_pars_2023-12-15.dat"
file.synth.points <- "Synthetic_environmental_data_2024-03-19.dat" # output of: PreprocessSwissInvertebrateEnvironmentalModellingData
file.selected.taxa      <- "selected_taxa_analysis.csv"

source(paste0(dir.utilities, "utilities_global.r"))
source("utilities_streambugs.r")
source("functions_run_streambugs.r")
source("library/application specific NIS/load_library.r")

# prepare inputs for plotting results on swiss maps
map.inputs        <- map.inputs(directory = paste0(dir.utilities,"swiss.map.gdb"))

## Analysis options ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# run simulations to produce ICE
ICE <- T # if False, simulations run for whole dataset
no.sites.ice <- 1
no.steps.ice <- 60
env.fact.ice <- "tempmaxC"

# global options
run.C        <- T # run C version of the rhs of Streambugs
res.add      <- F # calculate additional output (rates, limitation factors) when running Streambugs
plot.biomass <- T # plot time evolution of biomass (heavier computation for dataframes and plots if many sites and many taxa)

# has to be defined globally but is used in fct to construct parameters
inp <- NA # set NA for time-dependent input and parameters

# set-up output times
n.years       <- 20 # number of years
n.days.appart <- 1  # number of days from one time step to another (can be smaller than one day, e.g., 0.5)
tout          <- c(seq(0, n.years, by = n.days.appart/365)) # set up time vector of [n.years*365/n.days.appart] steps
n.time.steps  <- length(tout) # number of time steps simulated
cat("Simulating", n.years, "years every", n.days.appart, "day, total:", n.time.steps, "time steps")

# catchment scale variable by which we select the taxa pool (ask Rosi Siber for different catchments delimitation)
catch.variable <- "Watershed" 
# catch.variable <- "Watershed_Reg"

# tune parameters to get stronger presence/absence response to environmental factors, necessary for our application
# if Flag True, then Value is applied to the tuned parameter
par.adjust <- list("round.inv.traits"    = list("Flag" = F, "Value" = 0.1),  # round current and temperature preference if trait below Value (e.g., 0.1), to get stronger taxa (pres/abs) response
                   "curve.curr.temp"     = list("Flag" = T, "Value" = -15),  # assign Value (e.g., -10) to curve parameter of limitation factor of current and temperature, to get stronger (pres/abs) taxa
                   "curve.orgmic.sapro"  = list("Flag" = F, "Value" = 10),   # assign Value (e.g., -10) to curve parameter of limitation factor of orgmicropoll and saproby, to get stronger (pres/abs) taxa
                   "interc.orgmic.sapro" = list("Flag" = F, "Value" = 4),    # assign Value (e.g., -10) to curve parameter of limitation factor of orgmicropoll and saproby, to get stronger (pres/abs) taxa
                   "hdens"               = list("Flag" = F, "Value" = 0.05), # multiply hdens by small number (e.g., 0.05) to get strong self-inhibition, main impact on abundance and not pres/abs 
                   "thresh.indiv"        = list("Flag" = F, "Value" = 1))    # define presence data points above specific threshold Value (e.g., 1) of number of individuals

# catchment selection for simulation (can be specific or all to analyze taxa pool and results)
# select.catch <- "all"
# select.catch <- "Aare"
# select.catch <- "Doubs"
# select.catch <- "Limmat"
# select.catch <- "Reuss"
# select.catch <- "RheinabBS"
# select.catch <- "Rhone"
# select.catch <- "RhoneLeman"
select.catch <- "Ticino"

# specify how we select the taxa pool
# metric: based on prevalence or number of presence points in the catchment
# threshold: threshold above which taxa is selected
taxa.selection <- c("metric" = "Prevalence", 
                    "threshold" = 0.15)
# taxa.selection <- c("metric" = "NbPresPoints", 
#                     "threshold" = 7)

# select taxonomic levels to filter taxa pool
select.taxonomy = c("species", "genus", "family") 

# select colors for present/absent plots
vect.col.pres.abs <- c( "Present" = "#00BFC4", "Absent" = "#F8766D")
scales::show_col(vect.col.pres.abs)

# select specific site(s) for debugging and analysis
# select.sites   <- "random"
# select.sites <- "SynthPoint2612Aa"
select.sites <- "SynthPoint2533Ti"
# select.sites   <- "CH076ZGRe" # e.g., site with problem in Reuss (when running all reaches at once, 10 sites, 3 years)

# select number of sites per catchment for simulation
if(ICE){
  n.sites.per.catch    <- no.sites.ice * no.steps.ice
} else {
  n.sites.per.catch <- "all"
}
sites.selection      <- c("n.sites" = n.sites.per.catch, 
                          "select.sites" = select.sites)
correct.shade        <- 0.6 # factor to correct shading (fshade) for now, otherwise algae is too limited

# define run name for all output files
name.par.adjust <- ""
for (i in 1:length(par.adjust)) {
  name.par.adjust <- paste0(name.par.adjust, ifelse(par.adjust[[i]][["Flag"]], paste0(names(par.adjust)[i], par.adjust[[i]][["Value"]], "_"), ""))
}

name.run <- paste0(
  ifelse(ICE, "ICE_", ""),
  name.par.adjust,
  # ifelse(run.C, "runC_", "runR_"),
  # ifelse(res.add, "ResAdd_", ""),
  "Prev", taxa.selection[2], "_",
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
data.selected.taxa  <- read.csv(paste0(dir.utilities, file.selected.taxa), header = T, sep = ";", stringsAsFactors = F)


## Environmental data ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### preprocess dataset ####

# remove duplicated sites in data from midat
colnames(data.env.midat)[which(colnames(data.env.midat) == "SiteId")] <- "ReachID" # change column name of site id to match streambugs
data.env.midat <- data.env.midat[!duplicated(data.env.midat$ReachID), ] # remove duplicated sites

# bind env data from midat and synthetic points
colnames(data.env.synth)[which(colnames(data.env.synth) == "SiteId")] <- "ReachID" # change column name of site id to match streambugs
rind <- sample(1:dim(data.env.synth)[1], 1800) # subselect randomly some synthetic points otherwise too many of them
data.env.synth <- data.env.synth[rind,]
data.env.synth$MonitoringProgram <- "SyntheticPoint" # specify "fake" monitoring program
data.env.all <- bind_rows(data.env.midat, data.env.synth)
remove(data.env.synth) # remove entire synthetic sample points data from environement because slows down the system

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

# ### plot hist and maps ####
# 
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
  if(select.sites != "random"){
    rind.ice.site <- which(grepl(select.sites, data.env.inputs$ReachID))
  } else {
    rind.ice.site <- sample(1:nrow(data.env.inputs), no.sites.ice)
  }
  data.base.ice <- data.env.inputs[rind.ice.site,]
  
  # get range of env.factor
  min.env.fact <- min(data.env.inputs[,env.fact.ice])
  max.env.fact <- max(data.env.inputs[,env.fact.ice])
  
  # create range of all env.factor
  range.size <- max.env.fact-min.env.fact
  step.size  <- range.size/(no.steps.ice-1)
  vect.env.fact.ice <- seq(from=min.env.fact,
                           to=max.env.fact,
                           by=step.size)
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

n.sites <- ifelse(n.sites.per.catch == "all", 
                  dim(data.env.inputs)[1], 
                  length(vect.catch.select) * n.sites.per.catch) # compute nb total of sites for catchments selected
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
                                          catch.variable = catch.variable, vect.catch.select = vect.catch.select, dir.plots = dir.outputs)

data.prevalence <- list.df.prev.pres$Prevalence
data.nb.pres <- list.df.prev.pres$NbPresPoints
all.taxa <- data.prevalence$Taxon

if(taxa.selection["metric"] == "Prevalence"){
    data.taxa.selection <- data.prevalence
} else if(taxa.selection["metric"] == "NbPresPoints"){
    data.taxa.selection <- data.nb.pres
}

### get parameters for all taxa ####

# get parameters from different sources and make dataframes

# clean data of updated parameters from Vermeiren et al., 2021
for (i in 1:ncol(data.par.update)) {
    # i <- 1
    old.name <- colnames(data.par.update)[i]
    new.name <- if(old.name == "X"){ "Taxon" 
    } else if(old.name == "vco"){ "tempmaxtolval_class1"
    } else if(old.name == "cod"){ "tempmaxtolval_class2"
    } else if(old.name == "mod"){ "tempmaxtolval_class3"
    } else if(old.name == "war"){ "tempmaxtolval_class4"
    } else if(old.name == "st"){ "currenttolval_class1"
    } else if(old.name == "slow"){ "currenttolval_class2"
    } else if(old.name == "mod.1"){ "currenttolval_class3"
    } else if(old.name == "high"){ "currenttolval_class4"
    } else { old.name }
    colnames(data.par.update)[i] <- new.name
}

# extract trait scores and mass of all taxa
file.name <- paste0(
    ifelse(par.adjust[["round.inv.traits"]][["Flag"]], "round.inv.traits", "orig"), "_",
    "list.df.preferences_266Taxa.rds")
if(file.exists(paste0(dir.outputs, file.name))){
    list.df.preferences <- readRDS(paste0(dir.outputs, file.name))
    mass.all.inv <- list.df.preferences[["mass"]]
    vect.all.inv <- names(mass.all.inv)
} else {
    list.par.all.taxa <- construct.variables.par.catch(catch = taxa.selection["metric"], data.env.inputs = data.env.inputs, 
                                                       data.par.update = data.par.update, par.adjust = par.adjust,
                                                       data.taxa.selection = data.taxa.selection, taxa.selection = taxa.selection,
                                                       sites.selection = sites.selection, select.taxonomy = select.taxonomy,
                                                       catch.variable = catch.variable, name.run = name.run,
                                                       dir.inputs = dir.inputs, dir.plots = dir.outputs)
    # list.par.all.taxa$par.fixc[which(grepl("Limoniidae", names(list.par.all.taxa$par.fixc)) & grepl("class", names(list.par.all.taxa$par.fixc)))]
    system.def <- streambugs.get.sys.def(y.names=list.par.all.taxa$y.names$y.names, par=list.par.all.taxa$par.fixc)
    vect.all.inv <- list.par.all.taxa$Invertebrates
    list.df.preferences <- write.df.preferences(system.def = system.def,
                                                Invertebrates = vect.all.inv,
                                                dir.outputs = dir.outputs)
    mass.all.taxa <- system.def$par.taxaprop.direct$parvals[,"M"] # mean individual biomass of taxon i
    for (taxon in list.par.all.taxa$y.names$taxa) {
        names(mass.all.taxa)[which(grepl(taxon, names(mass.all.taxa)))] <- taxon
    }
    mass.all.taxa <- mass.all.taxa[which(!duplicated(names(mass.all.taxa)))]
    mass.all.inv <- mass.all.taxa[-c(1:3)]
    list.df.preferences[["mass"]] <- mass.all.inv
    saveRDS(list.df.preferences, file = paste0(dir.outputs, 
                                               ifelse(par.adjust[["round.inv.traits"]][["Flag"]], "round.inv.traits", "orig"), "_",
                                               "list.df.preferences_", length(vect.all.inv), "Taxa", ".rds"))
}

# make long dataframe trait preference
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
        long.df.trait <- rbind(long.df.trait, temp.df)
    }
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SIMULATIONS ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Run streambugs ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(name.run)

# construct variables and parameters per catchment
list.variables.par.catch <- lapply(vect.catch.select, FUN = construct.variables.par.catch,
                                   data.env.inputs = data.env.inputs, 
                                   data.par.update = data.par.update, par.adjust = par.adjust,
                                   data.taxa.selection = data.taxa.selection, 
                                   taxa.selection = taxa.selection, 
                                   sites.selection = sites.selection, 
                                   select.taxonomy = select.taxonomy, 
                                   catch.variable = catch.variable, 
                                   plot.foodweb = T, name.run = name.run,
                                   dir.inputs = dir.inputs, dir.plots = dir.outputs)

# streambugs.debug = 3
# # print metadata
# sink(file = paste0(dir.outputs, name.run, "debug.txt"))

# run streambugs
list.results <- lapply(list.variables.par.catch, FUN = run.streambugs.catch,
                       tout = tout, res.add = res.add, name.run = name.run, 
                       dir.output = dir.outputs, write.plot.results = F#,
                       # atol = 1e-06,
                       # rtol = 1e-06
                       )# ,
                       # hmax = 0.1)
# sink()

# print metadata
sink(file = paste0(dir.outputs, name.run, "metadata.txt"))
for (i in 1:length(list.results)) {
    # list.results$Aare$list.metadata.catch
    print(list.results[[i]][["list.metadata.catch"]])
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
}
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

if(plot.biomass){
    
    # make list of long dataframe biomass for each taxon
    list.df.res.taxa <- list()
    
    for (taxon in temp.vect.taxa) {
        # taxon <- temp.vect.taxa[1]
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
            q <- q + theme(legend.position="none")
        # }
        return(q)
    })
    # list.plots.all.sites[[2]]
    
    # print pdf with plots
    file.name <- paste0("_IndivBiomass_ggplot_allTaxa_gridSites")
    print.pdf.plots(list.plots = list.plots.sep.sites, width = 15, height = 12, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name,
                    png = F)
    file.name <- paste0("IndivBiomass_ggplot_allTaxa_allSites")
    # file.name <- paste0("IndivBiomass_ggplot_allTaxa_", select.sites)
    # file.name <- paste0("IndivBiomass_ggplot_allTaxa_", "CH070TGRh")
    
    print.pdf.plots(list.plots = list.plots.all.sites, width = 8, height = 6, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name,
                    png = F)
}


# Env data of modeled sites ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# extract env.data selected for running streambugs
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
# if(dim(data.env.res)[1] < 10){
#     p <- p + geom_text(size = 3, vjust = 0, nudge_y = 0.05)
# }
p <- p + theme_bw()
# p
file.name <- "EnvData_MultiDim"
print.pdf.plots(list.plots = list(p), 
                width = 10, height = 8, 
                dir.output = dir.outputs, info.file.name = name.run, file.name = file.name,
                png = F)

q <- ggplot(plot.data, aes(x = tempmaxC, y = currentms, 
                           color = orgmicropollTU, shape = saprowqclass)) #, label = ReachID))
q <- q + geom_point(alpha = 0.6)
q <- q + scale_color_continuous(trans='reverse')
q <- q + facet_wrap(. ~ Watershed)
q <- q + theme_bw()
q

file.name <- "EnvData_MultiDim_perCatch"
print.pdf.plots(list.plots = list(q), 
                width = 10, height = 8, 
                dir.output = dir.outputs, info.file.name = name.run, file.name = file.name,
                png = F)

# Additional results ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(res.add){
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
    mycolors <- sample(scales::hue_pal()(8), length(unlist(list.factors.type)), replace = T)
    
    processed.add.res <- get.plot.data.add.res(catch.results = list.results[[1]],
                                             list.factors.type = list.factors.type)
    list.df.res.add.taxa <- processed.add.res$list.df.add.res.taxa
    plot.data.all.taxa.sites <- processed.add.res$plot.data
    
    # plot all rates and limitation factors
    temp.vect.taxa <- unique(plot.data.all.taxa.sites$Taxon)
    # temp.vect.sites <- unique(plot.data$Site)
    list.plots <- list()
    for (taxon in temp.vect.taxa) {
        plot.data <- plot.data.all.taxa.sites %>%
            filter(Taxon == taxon) # %>%
        # filter(Site != 107)
        p <- ggplot(plot.data, aes(x = time, y = Value, color = Factor))
        p <- p + geom_line(size=1)
        p <- p + facet_grid(FactorType ~ Site, scales = "free")
        if(!grepl("Algae", taxon)){ p <- p + scale_color_manual(values = mycolors)}
        p <- p + theme_bw()
        p <- p + labs(title = taxon)
        # ggplotly(p)
        
        list.plots[[taxon]] <- p
    }
    file.name <- paste0("AddRes_", length(temp.vect.taxa), "Taxa")
    print.pdf.plots(list.plots = list.plots, width = 23, height = 8, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
                    png = F)
}

# Merge results in dataframes ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# source("utilities.r")

print(vect.analysis.env.fact)
vect.analysis <- c("ReachID", "X", "Y", "MonitoringProgram", catch.variable, vect.analysis.env.fact)

# make dataframes (wide and long) with steady state, abundance and prob of occ
# df.generic <- data.env.res[,vect.analysis]
df.generic <- data.env.inputs[,vect.analysis]

df.generic[,vect.occ.inv] <- NA

vect.names.output <- c("SteadyState", "Abundance", "ThreshPresAbs", "ProbObs", "SamplePresAbs")
list.wide.df.all.res <- list()
for (output in vect.names.output) {
    list.wide.df.all.res[[output]] <- df.generic
}

for (n in 1:length(list.results)) {
    # n <- 5
    temp.res <- as.matrix(list.results[[n]][["res"]])
    temp.parfix <- list.results[[n]][["par.fixc"]]
    temp.y.names <- list.results[[n]][["y.names"]]
    temp.inv <- list.results[[n]][["Invertebrates"]]
    temp.sites <- temp.y.names$reaches
    temp.res.ss <- calc.res.steadystate(res=temp.res, times=tout, par=temp.parfix, y.names=temp.y.names)
    for (site in temp.sites) {
        # site <- temp.sites[2]
        for (taxon in temp.inv) {
            # taxon <- temp.inv[2]
          occ.taxon <- names(vect.inv)[which(vect.inv == taxon)]
          ind <- which(grepl(taxon, names(temp.res.ss)) 
                       & grepl(paste0("^", site, "_"), # "^" Asserts that we are at the start because many sites could have same pattern
                               names(temp.res.ss)))
          temp.ss <- temp.res.ss[ind]
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

# recover observation data for selected sites
data.obs <- df.generic
data.inv.bdm$ReachID <- gsub("_", ".", data.inv.bdm$SiteId)
for (site in data.obs$ReachID) {
    # site <- df.generic$ReachID[1]
    temp.monit.prog <- data.obs[which(data.obs$ReachID == site), "MonitoringProgram"]
    # if(site %in% data.inv.bdm$ReachID){
        for(taxon in vect.inv){
          occ.taxon <- names(vect.inv)[which(vect.inv == taxon)]
            if(temp.monit.prog == "BDM"){
                obs <- data.inv.bdm[which(data.inv.bdm$ReachID == site), which(grepl(taxon, colnames(data.inv.bdm)))]
            } else {
                cind.taxon <- which(grepl(taxon, colnames(data.env.inputs)))
                if(length(cind.taxon) > 0){ # if tayon exists in all monit. prog. dataset
                    obs <- data.env.inputs[which(data.env.inputs$ReachID == site), cind.taxon]
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
long.df.res.taxa <- list.long.df.res.taxa %>% reduce(merge, by = c(vect.analysis, "Taxon"))
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

# ICE streambugs ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# data.base.ice <- read.csv(paste0(dir.outputs, "WideDataInputs_ICE_100Sites.csv"), sep = ";")
# data.results.ice <- read.csv(paste0(dir.outputs,
# "10Catch_2000Sites_ICE_Prev0.15_10Yea_3651Steps_WideData_ResultsProbObs.csv"), sep = ";")

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

temp.vect.taxa <- names(selected.taxa.analysis)[names(selected.taxa.analysis) %in% vect.taxa]
temp.occ.taxa <- paste0("Occurrence.", temp.vect.taxa)

list.plots <- list()

for(taxon in temp.occ.taxa){
  
  # taxon.under.obs <- names(taxa.colnames)[4]
  
  pred.ice.streamb <- data.results.ice[,c("ReachID","Watershed", select.env.fact, taxon)] %>%
    mutate(observation_number = rep(1:no.sites.ice, each = no.steps.ice),
           column_label = 1,
           model = "Streambugs") %>%
    rename(pred = taxon)
  
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
  
  fig3 <- ggplot(data=plot.data) +
    geom_line(aes(x=.data[[select.env.fact]],
                  y=pred,
                  group=observation_number, 
                  color=as.character(observation_number)),
              show.legend = FALSE) +
    # geom_line(data=pred.ice.mean,
    #           aes(x=.data[[select.env.fact]], y=avg),
    #           size=1.5) +
    geom_rug(data = env.factor.sampled,
             aes(x=variable), 
             color="grey20",
             alpha=0.7,
             inherit.aes=F) + 
    # geom_segment(data=pred.ice.mean.bounds,
    #              inherit.aes = FALSE,
    #              lineend="round",
    #              linejoin="round",
    #              aes(x=x.mean,
    #                  y=y.mean.min,
    #                  xend=x.mean,
    #                  yend=y.mean.max),
    #              arrow=arrow(length = unit(0.3, "cm"),
    #                          ends = "both")) +
    # facet_wrap(~factor(model, levels = c("GLM", "GAM", "RF"))) +
  # facet_wrap(~Watershed) +
  ylim(0,1) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          legend.title = element_text(size=24),
          legend.text = element_text(size=20)) +
    labs(title = taxon,
         x = "Temperature",
         y = "Predicted probability of occurrence")
  # fig3
  list.plots[[taxon]] <- fig3
}

file.name <- paste0("ice_", ifelse(no.sites.ice == 1, select.sites, no.sites.ice), "_Sites_", no.steps.ice, "Steps_", select.env.fact)
print.pdf.plots(list.plots = list.plots, width = 6, height = 6, 
                dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
                png = F)


# # !!!! analysis baetis rhodani ####
# m.rhodani <- mass.all.inv[4]
# vect.ss.rhodani <- list.wide.df.all.res$SteadyState$Occurrence.Baetisrhodani
# plot(data.ICE$tempmaxC, vect.ss.rhodani)
# vect.abund.rhodani <- vect.ss.rhodani / m.rhodani
# plot(data.ICE$tempmaxC, vect.abund.rhodani)
# vect.prob.rhodani <- sapply(vect.ss.rhodani, FUN = get.prob.obs, m.rhodani)
# plot(data.ICE$tempmaxC, vect.prob.rhodani)


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
# p <- ggplot(plot.data, aes(x = Taxon, y = Abundance))
# p <- p + geom_boxplot()
# p <- p + geom_point(data = temp.df.prev, aes(x = Taxon, y = ScaledPrevalence), color = "green", size = 2)
# p <- p + theme_bw()
# p <- p + scale_x_discrete(limits = temp.df.prev$Taxon)
# p <- p + theme(axis.text.x = element_text(angle=90))
# p


## abund vs env fact ####


### resp abund to one env fact ####

list.plots <- list()
for(taxon in vect.inv){
    print(taxon)
    plot.data <- long.df.res.taxa.env.fact %>%
        filter(Taxon %in% taxon) %>%
        filter(Factor %in% vect.analysis.env.fact)
    
    p <- ggplot(plot.data, aes(x = Value, y = Abundance))
                               # , color = Observation))
    p <- p + geom_point()
    p <- p + facet_wrap(. ~ Factor, scales = "free_x")
    # p <- p + scale_y_continuous(limits = c(0,1))
    p <- p + labs(title = taxon)
    p <- p + theme_bw()
    # p <- p + scale_color_manual(values=vect.col.pres.abs)
    # p
    list.plots[[taxon]] <- p
}
# list.plots$Ceratopogonidae

file.name <- paste0("AbundanceVsEnvFact")
print.pdf.plots(list.plots = list.plots, width = 12, height = 8, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
                png = F)


### resp prob occ to direct env fact ####

list.plots <- list()
for(taxon in vect.inv){
    # taxon <- vect.inv[30]
    print(taxon)
    plot.data <- long.df.all %>%
      filter(Taxon %in% taxon) %>%
      filter(Factor %in% vect.dir.env.fact) %>%
      filter(Output %in% vect.names.output[3:5])
    plot.data$Result <- as.numeric(plot.data$Result)
    plot.data$Dispersal <- ifelse(is.na(plot.data$Result), "Yes", "No")
    plot.data$Result[is.na(plot.data$Result)] <- 0
    plot.data$Factor <- factor(plot.data$Factor, levels=vect.dir.env.fact)
    # summary(plot.data$Result)
    
    plot.data.trait <- filter(long.df.trait, Taxon %in% taxon)
    plot.data.trait$Factor <- factor(plot.data.trait$Factor, levels=vect.dir.env.fact)
    
    p <- ggplot(plot.data, aes(x = Value, y = Result, color = Dispersal))
    p <- p + geom_point(alpha = 0.5)
    p <- p + geom_point(data = plot.data.trait, aes(x = Value, y = Pref), shape = 4, color = "#7CAE00", size = 3)
    p <- p + facet_grid(Output ~ Factor, scales = "free_x")
    p <- p + scale_y_continuous(limits = c(0,1))
    p <- p + labs(title = taxon)
    p <- p + theme_bw()
    p <- p + scale_color_manual(values=c("Yes" = "grey60", "No" = "orange"))
    # p
    list.plots[[taxon]] <- p
}
# list.plots$Ceratopogonidae

file.name <- paste0("PresAbsVsDirEnvFact")
print.pdf.plots(list.plots = list.plots, width = 8, height = 6, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
                png = F)

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
# 
# file.name <- paste0("Maps_ProbObs_vs_Obs")
# print.pdf.plots(list.plots = list.plots, width = 10, height = 8, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
#                 png = F)
# 
# 
