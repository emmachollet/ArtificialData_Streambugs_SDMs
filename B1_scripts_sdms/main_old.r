###############################################################################|
#                                                                              #
#  --- Study of the influence of noise on overfitting of machine learning ---  #
#            --- and statistical species distribution models. ---              #
#                                                                              #
#                          --- February 2023 ---                               #
#                                                                              #
#                   --- Emma Chollet, Gaspard Fragnière ---                    #
#                --- Andreas Scheidegger and Nele Schuwirth ---                #
#                                                                              #
#                      --- emma.chollet@eawag.ch ---                           #
#                                  l                                           #
###############################################################################|

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# PRELIMINARIES ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Sys.setenv(LANG="EN")
set.seed(13)   # Always set seed to a lucky number
getwd()        # Show working directory. It needs to be the location of 'main.r'
rm(list=ls())  # Free work space
graphics.off() # Clean graphics display


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Libraries ####

if (!require("dplyr")){install.packages("dplyr"); library("dplyr")}                    # to sort, join, merge data
if ( !require("tidyr") ) { install.packages("tidyr"); library("tidyr") }               # to sort, join, merge data
if (!require("readr")){install.packages("readr"); library("readr")}
if (!require("jsonlite")){install.packages("jsonlite"); library("jsonlite")}
if (!require("pROC")){install.packages("pROC"); library("pROC")}  # to compute AUC                  
if ( !require("sf") ) { install.packages("sf"); library("sf") }                        # to read layers for plotting maps

# specific for Neural Network
if (!require("reticulate")){install.packages("reticulate"); library("reticulate")}
#install_miniconda()              # run this the very first time reticulate is installed
#install.packages("tensorflow")
library("tensorflow")
#install_tensorflow()             # run this line only when opening new R session
#install.packages("keras")
library("keras")
#install_keras()                  # run this line only when opening new R session
#use_condaenv()
# reticulate::install_python(version = '<version>')

# caret has to be loaded at the end to not cache function 'train'
if (!require("caret")){install.packages("caret"); library("caret")}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Directories, files and functions ####

# define directories
dir.input.data          <- "../A3_outputs_streambugs/"
dir.output              <- "../B2_outputs_sdms/"
dir.utilities           <- "../00_utilities/"
# dir.create(dir.output)

# define files
name.streambugs.run     <- "10Catch_3081Sites_Prev0.15_10Yea_3651Steps"
file.input.data         <- paste0(name.streambugs.run, "_WideData_ResultsThreshPresAbs.csv")
file.prev.taxa          <- paste0(name.streambugs.run, "_PrevalenceTaxonomy_ThreshPresAbs.csv")
file.selected.taxa      <- "selected_taxa_analysis.csv"



# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Load data and functions ####

# read data
data.env.taxa          <- read.csv(paste0(dir.input.data, file.input.data), header = T, sep = ";", stringsAsFactors = F)
data.prev.taxa      <- read.csv(paste0(dir.input.data, file.prev.taxa), header = T, sep = ";", stringsAsFactors = F)
data.selected.taxa  <- read.csv(paste0(dir.utilities, file.selected.taxa), header = T, sep = ";", stringsAsFactors = F)

# load functions
source(paste0(dir.utilities, "utilities_global.r"))
source("performances_assessment.r")
source("training_pipeline.r")
source("ml_models.r")
source("plotting.r")

# prepare inputs for plotting results on swiss maps
map.inputs      <- map.inputs(directory = paste0(dir.utilities,"swiss.map.gdb"))


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Inputs, training options and noise scenarios ####

# environmental factors used for streambugs
vect.dir.env.fact <- c("Temperature"      = "tempmaxC", 
                       "Flow velocity"    = "currentms", 
                       "Toxic units"      = "orgmicropollTU", 
                       "Saproby"          = "saprowqclass")
vect.indir.env.fact <- c("Shade"          = "shade", 
                         "Litter input"   = "Lit_Inp" ,        
                         "Phosphorus"     = "C_P", 
                         "Nitrate"        = "C_N")
vect.info <- colnames(data.env.taxa)[1:5]
print(vect.info)

# environmental factors selected for sdms
env.factor <- c(vect.dir.env.fact, vect.indir.env.fact)
env.factor.full <- env.factor
no.env.fact <- length(env.factor)

# taxa list full
cind.taxa <- which(grepl("Occurrence.", colnames(data.env.taxa)))
vect.taxa.full <- colnames(data.env.taxa)[cind.taxa]
names(vect.taxa.full) <- gsub("Occurrence.", "", vect.taxa.full)
no.taxa.full <- length(vect.taxa.full)
print(vect.taxa.full)

# taxa list selected for analysis and above prevalence threshold
prev.threshold <- 0.1
temp.occ.taxa <- data.prev.taxa[which(data.prev.taxa$Prevalence > prev.threshold), "Occurrence.taxa"]
selected.taxa.analysis <- data.selected.taxa$Occurrence.taxa
taxa.colnames <- selected.taxa.analysis[which(selected.taxa.analysis %in% temp.occ.taxa)]
names(taxa.colnames) <- gsub("Occurrence.", "", taxa.colnames)
no.taxa <- length(taxa.colnames)
print(taxa.colnames)

# select colors for present/absent plots
vect.col.pres.abs <- c( "Present" = "#00BFC4", "Absent" = "#F8766D")
scales::show_col(vect.col.pres.abs)

# sdms training options
number.split            <- 3
split.criterion         <- "ReachID"
number.sample           <- 1000
sdm.models              <- list("GLM" = "glm",
                                "GAM" = "gamLoess",
                                "RF"  = "rf")#,
                                # "ann")
no.models <- length(sdm.models)
models <- append(list("null"), sdm.models)

# noise scenarios
list.noise <- list(
  
  "noise.disp"     = list("type"       = "dispersal",
                         "target"     = NULL,
                         "amount"     = NULL,
                         "parameters" = NULL) # ,

  # "noise.temp"     = list("type"       = "gaussian",
  #                         "target"     = "temperature",
  #                         "amount"     = 5,
  #                         "parameters" = list("min"=0, "max"=35)),
  # 
  # "noise.gamm"     = list("type"       = "missdetection",
  #                         "target"     = "Occurrence.Gammaridae",
  #                         "amount"     = 0.1,
  #                         "parameters" = NULL),
  # 
  # "noise.add.fact" = list("type"       = "add_factor",
  #                         "target"     = "random1", # new factor name
  #                         "amount"     = NULL,
  #                         "parameters" = NULL),
  # 
  # "noise.rem.fact" = list("type"       = "remove_factor",
  #                         "target"     = "temperature",
  #                         "amount"     = NULL,
  #                         "parameters" = NULL),
  
)

# define experiment name
experiment.name <- paste0(
  number.sample, "Sites_",
  no.taxa, "Taxa_",
  no.models, "SDMs_",
  no.env.fact, "EnvFact_",
  paste(names(list.noise), collapse = "_")
)
print(experiment.name)

# create output directory for experiment results
dir.experiment          <- paste0(dir.output, experiment.name, "/")
dir.metadata            <- paste0(dir.experiment, "metadata.json")
dir.create(dir.experiment)

# saving metadata to JSON file
metadata.list           <- list("input_data"       = file.input.data,
                                "prev_taxa"        = file.prev.taxa,
                                "experiment.name"  = experiment.name,
                                "number_split"     = number.split,
                                "split.criterion"  = split.criterion,
                                "number.sample"    = number.sample,
                                "models"           = models,
                                "env_factor"       = env.factor,
                                "env_factor_full"  = env.factor.full,
                                "noise"            = list.noise) 

metadata.json           <- toJSON(metadata.list)
write(metadata.json, dir.metadata)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# PREPROCESS DATA ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# preprocess input dataset
data.input <- data.env.taxa[,c(vect.info,env.factor,taxa.colnames)] # subselect columns of interest
data.input[which(is.na(data.input$MonitoringProgram)), "MonitoringProgram"] <- "SyntheticPoint" # replace empty monitoring programs
catch.variable <- "Watershed"

# sanity check of NAs
sum(is.na(data.input[,c(vect.info,env.factor)])) # should be zero 
sum(is.na(data.input[,c(taxa.colnames)])) # could have a lot, decide to replace by zero or not

if("noise.disp" %in% names(list.noise)){
  data.input[is.na(data.input)] <- 0 # replace NAs by 0
}

# replace 0/1 by absent/present
for (taxon in taxa.colnames) {
  data.input[which(data.input[,taxon] == 0),taxon] <- "Absent"
  data.input[which(data.input[,taxon] == 1),taxon] <- "Present"
  data.input[,taxon] = as.factor(data.input[,taxon])
}

# make long dataframe for analysis
long.df.input.taxa <- gather(data.input, key = Taxon, value = Observation, all_of(taxa.colnames))
long.df.input.taxa$Taxon <- gsub("Occurrence.", "", long.df.input.taxa$Taxon)
long.df.input.env.taxa <- gather(long.df.input.taxa, key = Factor, value = Value, all_of(env.factor))


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Prevalence analysis ####

# plot bars pres/abs/NA
plot.data <- long.df.input.taxa
plot.data$Observation <- factor(plot.data$Observation, levels=c("Present", NA, "Absent"))
temp.taxa <- unique(plot.data$Taxon)
temp.catch <- unique(plot.data$Watershed)

# analyze distribution of prevalence and NAs across catchments
p <- ggplot(plot.data, aes(x = Taxon, fill = Observation))
p <- p + geom_bar(stat="count", position='stack')
p <- p + facet_wrap(~ Watershed, scales = "free_y")
p <- p + scale_fill_manual(values=vect.col.pres.abs)
p <- p + theme_bw()
p <- p + scale_x_discrete(limits = temp.taxa)
p <- p + theme(axis.text.x = element_text(angle=90))
p

# analyze distribution of prevalence and NAs across taxa
q <- ggplot(plot.data, aes(x = Watershed, fill = Observation))
q <- q + geom_bar(stat = "count", position = "stack")
q <- q + facet_wrap(~ Taxon, scales = "free_y")
q <- q + scale_fill_manual(values=vect.col.pres.abs)
q <- q + theme_bw()
q <- q + scale_x_discrete(limits = temp.catch)
q <- q + theme(axis.text.x = element_text(angle=90))
q

# print plots in pdf
file.name <- "GeomBar_PresAbs_TaxaCatch"
print.pdf.plots(list.plots = list(q,p), width = 12, height = 9, dir.output = dir.experiment, info.file.name = "", file.name = file.name, 
                png = F)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Split and standardize data ####

preprocessed.data.cv  <- preprocess.data(data=data.input,
                                         env.fact=env.factor,
                                         dir=dir.experiment,
                                         split.type="CV",
                                         nb.split=number.split,
                                         splitting.criterion=split.criterion)

preprocessed.data.fit <- preprocess.data(data=data.input,
                                         env.fact=env.factor,
                                         dir=dir.experiment,
                                         split.type="FIT")


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# APPLY MODELS ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Cross-validation ####

models.cv <- lapply(sdm.models, FUN = apply.caret.model, 
                    preprocessed.data.cv, split.type = "CV", 
                    taxa.colnames, env.factor, list.noise)

save.models(models=models.cv,
            path=dir.experiment,
            split.type="CV")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Fit entire dataset ####

models.fit <- lapply(sdm.models, FUN = apply.caret.model, 
                     preprocessed.data.fit, split.type = "FIT", 
                     taxa.colnames, env.factor, list.noise)

save.models(models=models.fit,
            path=dir.experiment,
            split.type="FIT")

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# PLOTS EXPERIMENT ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# models.fit      <- load.models(path=dir.experiment,
#                                split.type="FIT")
# models.cv       <- load.models(path=dir.experiment, 
#                                split.type="CV")

std.const.cv    <- readRDS(file=paste0(dir.experiment, "standardization_constant_CV.rds"))
std.const.fit   <- readRDS(file=paste0(dir.experiment, "standardization_constant_FIT.rds"))
prev.taxa <- data.prev.taxa

# Variable for plots ----
all.results <- summarize.all.results(models.cv, prev.taxa)

ggplot.all.results <- restructure.all.results(all.results)

color.map <- c('GLM'           = 'deepskyblue',   # Generalized Linear Model
               'GAM'      = 'green',         # Generalized Additive Model
               # 'ann'           = 'orange',        # Artificial Neural Network
               'RF'            = 'red')#,           # Random Forest
               # 'null'          = 'black')         # Null Model


# Fig. 1 ----
gg.plot.dev.infos <-  ggplot.all.results %>%
  group_by(fit_pred, model) %>%
  summarise(mean=mean(dev),
            max=max(dev),
            min=min(dev),
            median=median(dev))

median.null.pred <- min((gg.plot.dev.infos %>%
                           filter(fit_pred=="pred" & model=="null"))["median"])
median.best.pred <- min((gg.plot.dev.infos %>%
                           filter(fit_pred=="pred"))["median"])

fig1 <- ggplot(data=ggplot.all.results,
               aes(x=model, y=dev)) +
  geom_boxplot(aes(fill=fit_pred)) +
  # geom_hline(aes(yintercept=median.null.pred,
  #                color="median null model (pred)"),
  #            linetype="dashed") + 
  geom_hline(aes(yintercept=median.best.pred,
                 color="median best model (pred)"),
             linetype="dashed") +
  scale_colour_manual(values = c('green','black')) +
  theme_minimal() +
  theme(legend.title=element_blank())

pdf(paste0(dir.experiment, "boxplot_performance.pdf"))
print(fig1)
dev.off()


# Fig. 2 ----
fig2 <- ggplot(data=subset(ggplot.all.results, model %in% "null"),
               aes(x=prevalence, y=dev, color=model)) +
  geom_smooth(se = FALSE) + 
  geom_point(data=ggplot.all.results) +
  facet_wrap(~fit_pred) +
  scale_color_manual(values=color.map) +
  theme_minimal() + 
  theme(legend.title=element_blank())

pdf(paste0(dir.experiment, "bellplot_performance.pdf"))
print(fig2)
dev.off() 

# Fig. 3 ----

list.plots <- list()

for(taxon.under.obs in names(taxa.colnames)){
  
  cat("Plotting ICE for:", taxon.under.obs)
  
  # env.fact.under.obs <- "temperature"
  env.fact.under.obs <- "tempmaxC"
  
  # taxon.under.obs <- names(taxa.colnames)[1]
  
  # input.env.factors <- get.input.env.factors(experiment.name)
  input.env.factors <- env.factor
  
  ice.dfs <- plot.ice(models.performance=models.fit,
                      env.factor=env.fact.under.obs,
                      taxa=taxon.under.obs,
                      standardization.constant=std.const.fit[[1]],
                      observations = preprocessed.data.fit$`Entire dataset`,
                      nb.sample=100,
                      resolution=200,
                      input.env.factors=input.env.factors)
  
  observations              <- ice.dfs[["observations"]]
  env.factor.sampled        <- ice.dfs[["env.factor.sampled"]]
  observations.mean         <- observations %>%
    group_by(across(all_of(env.fact.under.obs)), model) %>%
    summarise(avg = mean(pred))
  observations.mean.bounds  <- observations.mean %>% group_by(model) %>%
    summarise(x.mean=max(across(all_of(env.fact.under.obs))),
              y.mean.min=min(avg),
              y.mean.max=max(avg))
  
  fig3 <- ggplot(data=observations) +
    geom_line(aes(x=.data[[env.fact.under.obs]],
                  y=pred,
                  group=observation_number, 
                  color=as.character(observation_number)),
              show.legend = FALSE) +
    geom_line(data=observations.mean,
              aes(x=.data[[env.fact.under.obs]], y=avg),
              size=1.5) +
    geom_rug(data = env.factor.sampled,
             aes(x=variable), 
             color="grey20",
             alpha=0.7,
             inherit.aes=F) + 
    geom_segment(data=observations.mean.bounds,
                 inherit.aes = FALSE,
                 lineend="round",
                 linejoin="round",
                 aes(x=x.mean,
                     y=y.mean.min,
                     xend=x.mean,
                     yend=y.mean.max),
                 arrow=arrow(length = unit(0.3, "cm"),
                             ends = "both")) +
    facet_wrap(~model)
  
  list.plots[[taxon.under.obs]] <- fig3
  
  pdf(paste0(dir.experiment, "ice_", taxon.under.obs, ".pdf"))
  print(fig3)
  dev.off()
}  

file.name <- paste0("ice_alltaxa")
print.pdf.plots(list.plots = list.plots, width = 15, height = 8, 
                dir.output = dir.experiment, info.file.name = "", file.name = file.name, 
                png = F)



# # Fig 4
# noised.dataset <- models.fit[[1]][[1]][[1]][[1]][["observation"]]
# taxa.col <- colnames(noised.dataset)[which(grepl("Occurrence.", colnames(noised.dataset)))]
# 
# noised.prevalence <- lapply(taxa.col, FUN=function(taxon){
#   occurence  <- noised.dataset[[taxon]]
#   prevalence <- length(occurence[occurence=="present"])/length(occurence) 
#   return(prevalence)
# })  
# 
# noised.prev.taxa <- prev.taxa[which(prev.taxa$Occurrence.taxa %in% taxa.colnames),]  
# noised.prev.taxa[["Prevalence"]] <- noised.prevalence
# 
# prev.taxa[["Status"]] <- "unnoised"
# noised.prev.taxa[["Status"]] <- "noised"
# 
# full.prev.taxa <- rbind(prev.taxa, noised.prev.taxa)
# full.prev.taxa[["Occurrence.taxa"]] <- as.factor(full.prev.taxa[["Occurrence.taxa"]])
# full.prev.taxa$Prevalence <- unlist(full.prev.taxa$Prevalence)
# full.prev.taxa$Occurrence.taxa <- gsub("Occurrence.", "", full.prev.taxa$Occurrence.taxa)
# 
# fig4 <- ggplot(data=full.prev.taxa,
#                aes(x=Prevalence,
#                    y=Occurrence.taxa,
#                    color=Status)) + 
#   geom_point() +
#   scale_y_discrete(limits = rev(gsub("Occurrence.", "", prev.taxa$Occurrence.taxa))) +
#   labs(x="Prevalence", y="Taxa")
# theme_minimal() + 
#   theme(legend.title=element_blank())
# 
# 
# pdf(paste0(dir.experiment, "change_prevalence.pdf"))
# print(fig4)
# dev.off()


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# PLOTS COMPARISON ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Creating the plots for the different experiments ----


# compare experiment dispersal noise
list.exp.disp     <- list("4envfact"          = "1000Sites_15Taxa_3SDMs_4EnvFact_",
                          "8envfact"          = "1000Sites_15Taxa_3SDMs_8EnvFact_",
                          "4envfact disp noise"         = "1000Sites_15Taxa_3SDMs_4EnvFact_noise.disp",
                          "8envfact disp noise"        = "1000Sites_15Taxa_3SDMs_8EnvFact_noise.disp")

create.comparison.plots("gauss_temp", list.exp.gauss, 
                        file.prev.taxa="8Catch_1416Sites_curve.curr.temp-10_interc.orgmic.sapro4_10Yea_3651Steps_PrevalenceTaxonomy_ThreshPresAbs.csv",
                        taxon="Gammaridae",
                        env.fact="tempmaxC")















# Noise: baseline (no noise) ----
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


experiment.name         <- paste0("baseline_",
                                  format(Sys.time(), "%d_%m_%Y_%Hh%M"))
number.split            <- 3
split.criterion         <- "ReachID"
number.sample           <- 1500
models                  <- list("null", 
                                "glm",
                                "gamloess",
                                "rf")#,
                                # "ann")

# four direct environmental factors streambugs
env.factor <- c("Temperature"      = "tempmaxC", 
                "Flow velocity"    = "currentms", 
                "Toxic units"      = "orgmicropollTU", 
                "Saproby"          = "saprowqclass")
env.factor.full <- env.factor

# original environmental factors
# env.factor              <- c("Temperature"                      = "temperature",       # Temp
#                              "Flow velocity"                    = "velocity",          # FV
#                              "Riparian agriculture"             = "A10m",              # A10m
#                              "Livestock unit density"           = "cow.density",       # LUD
#                              "Insecticide application rate"     = "IAR",               # IAR
#                              "Urban area"                       = "urban.area",        # Urban
#                              "Forest-river intersection"        = "FRI",               # FRI
#                              "Forest-river intersection buffer" = "bFRI",              # bFRI
#                              "Width variability"                = "width.variability") # WV
# 
# env.factor.full         <- c(env.factor,
#                              "Temperature2"                     = "temperature2",
#                              "Velocity2"                        = "velocity2")

# all environmental factors streambugs
# env.factor <- c("Temp_average", 
#                 "tempmaxC", 
#                 "L", "w", "I0", "shade",
#                 "C_P", "C_N", "C_SusPOM", "currentms", 
#                 "orgmicropollTU", 
#                 "saprowqclass", "Lit_Inp" , 
#                 "microhabaf_type1", "microhabaf_type2", "microhabaf_type3", "microhabaf_type4"
# )



noise1 <- list("type"       = "gaussian",
               "target"     = "temperature",
               "amount"     = 5,
               "parameters" = list("min"=0, "max"=35))

noise2 <- list("type"       = "missdetection",
               "target"     = "Occurrence.Gammaridae",
               "amount"     = 0.1, 
               "parameters" = NULL)

noise3 <- list("type"       = "add_factor",
               "target"     = "random1", # new factor name
               "amount"     = NULL, 
               "parameters" = NULL)

noise4 <- list("type"       = "remove_factor",
               "target"     = "temperature",
               "amount"     = NULL, 
               "parameters" = NULL)

noise  <- list("noise1"     = noise1,
               "noise2"     = noise2,
               "noise3"     = noise3,
               "noise4"     = noise4)

noise  <- list()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# training model
training.pipeline(file.input.data=file.input.data,
                  file.prev.taxa=file.prev.taxa,
                  experiment.name=experiment.name,
                  number.split=number.split,
                  split.criterion=split.criterion,
                  number.sample=number.sample,                   
                  models=models,
                  noise=noise,
                  env.factor=env.factor,
                  env.factor.full=env.factor.full)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# plotting models and saving model
create.plots(experiment.name=experiment.name,
             file.prev.taxa=file.prev.taxa)



# Noise: adding gaussian noise to environmental factor "temperature" ----
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
for (sd in c(1,3,5)){
  
  
  experiment.name         <- paste0("gaussian_temp_", sd, "_",
                                    format(Sys.time(), "%d_%m_%Y_%Hh%M"))
  number.split            <- 3
  split.criterion         <- "ReachID"
  number.sample           <- 1500
  models                  <- list("null", 
                                  "glm",
                                  "gamloess",
                                  "rf")#,
  # "ann")
  
  # four direct environmental factors streambugs
  env.factor <- c("Temperature"      = "tempmaxC", 
                  "Flow velocity"    = "currentms", 
                  "Toxic units"      = "orgmicropollTU", 
                  "Saproby"          = "saprowqclass")
  env.factor.full <- env.factor
  
  noise1 <- list("type"       = "gaussian",
                 "target"     = "tempmaxC",
                 "amount"     = sd,
                 "parameters" = list("min"=0, "max"=35))
  
  noise  <- list("noise1"     = noise1)

  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # training model
  training.pipeline(file.input.data=file.input.data,
                    file.prev.taxa=file.prev.taxa,
                    experiment.name=experiment.name,
                    number.split=number.split,
                    split.criterion=split.criterion,
                    number.sample=number.sample,                   
                    models=models,
                    noise=noise,
                    env.factor=env.factor,
                    env.factor.full=env.factor.full)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # plotting models and saving model
  create.plots(experiment.name=experiment.name,
               file.prev.taxa=file.prev.taxa)
  
}


# Noise: missdetecting taxon Gammaridae ----
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
for (p in c(10,25,50)){
  
  
  experiment.name         <- paste0("gammaridae_missdetection_", p, "_",
                                    format(Sys.time(), "%d_%m_%Y_%Hh%M"))
  number.split            <- 3
  split.criterion         <- "SiteId"
  number.sample           <- 10000
  models                  <- list("null", "glm", "gamloess", "rf", "ann")
  env.factor              <- c("Temperature"                      = "temperature",       # Temp
                               "Flow velocity"                    = "velocity",          # FV
                               "Riparian agriculture"             = "A10m",              # A10m
                               "Livestock unit density"           = "cow.density",       # LUD
                               "Insecticide application rate"     = "IAR",               # IAR
                               "Urban area"                       = "urban.area",        # Urban
                               "Forest-river intersection"        = "FRI",               # FRI
                               "Forest-river intersection buffer" = "bFRI",              # bFRI
                               "Width variability"                = "width.variability") # WV
  
  env.factor.full         <- c(env.factor,
                               "Temperature2"                     = "temperature2",
                               "Velocity2"                        = "velocity2")
  
  
  noise1 <- list("type"       = "missdetection",
                 "target"     = "Occurrence.Gammaridae",
                 "amount"     = p/100.0, 
                 "parameters" = NULL)

  
  noise  <- list("noise1"     = noise1)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # training model
  training.pipeline(file.input.data=file.input.data,
                    file.prev.taxa=file.prev.taxa,
                    experiment.name=experiment.name,
                    number.split=number.split,
                    split.criterion=split.criterion,
                    number.sample=number.sample,                   
                    models=models,
                    noise=noise,
                    env.factor=env.factor,
                    env.factor.full=env.factor.full)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # plotting models and saving model
  create.plots(experiment.name=experiment.name,
               file.prev.taxa=file.prev.taxa)
  
}


# Noise: missdetecting all taxa ----
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
for (p in c(10,25,50)){
  
  
  experiment.name         <- paste0("all_taxa_missdetection_", p, "_",
                                    format(Sys.time(), "%d_%m_%Y_%Hh%M"))
  number.split            <- 3
  split.criterion         <- "SiteId"
  number.sample           <- 10000
  models                  <- list("null", "glm", "gamloess", "rf", "ann")
  env.factor              <- c("Temperature"                      = "temperature",       # Temp
                               "Flow velocity"                    = "velocity",          # FV
                               "Riparian agriculture"             = "A10m",              # A10m
                               "Livestock unit density"           = "cow.density",       # LUD
                               "Insecticide application rate"     = "IAR",               # IAR
                               "Urban area"                       = "urban.area",        # Urban
                               "Forest-river intersection"        = "FRI",               # FRI
                               "Forest-river intersection buffer" = "bFRI",              # bFRI
                               "Width variability"                = "width.variability") # WV
  
  env.factor.full         <- c(env.factor,
                               "Temperature2"                     = "temperature2",
                               "Velocity2"                        = "velocity2")
  
  
  noise <- lapply(TAXA.COLNAMES, FUN=function(taxon){
    noise_taxon <- list("type"       = "missdetection",
                        "target"     = taxon,
                        "amount"     = p/100.0, 
                        "parameters" = NULL)
    
    return(noise_taxon)
  })
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # training model
  training.pipeline(file.input.data=file.input.data,
                    file.prev.taxa=file.prev.taxa,
                    experiment.name=experiment.name,
                    number.split=number.split,
                    split.criterion=split.criterion,
                    number.sample=number.sample,                   
                    models=models,
                    noise=noise,
                    env.factor=env.factor,
                    env.factor.full=env.factor.full)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # plotting models and saving model
  create.plots(experiment.name=experiment.name,
               file.prev.taxa=file.prev.taxa)
  
}


# Noise: Reducing the number of datapoints ----
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
for (s in c(2000,1000,500)){
  
  
  experiment.name         <- paste0("subset_", s, "_",
                                    format(Sys.time(), "%d_%m_%Y_%Hh%M"))
  number.split            <- 3
  split.criterion         <- "SiteId"
  number.sample           <- s
  models                  <- list("null", "glm", "gamloess", "rf", "ann")
  env.factor              <- c("Temperature"                      = "temperature",       # Temp
                               "Flow velocity"                    = "velocity",          # FV
                               "Riparian agriculture"             = "A10m",              # A10m
                               "Livestock unit density"           = "cow.density",       # LUD
                               "Insecticide application rate"     = "IAR",               # IAR
                               "Urban area"                       = "urban.area",        # Urban
                               "Forest-river intersection"        = "FRI",               # FRI
                               "Forest-river intersection buffer" = "bFRI",              # bFRI
                               "Width variability"                = "width.variability") # WV
  
  env.factor.full         <- c(env.factor,
                               "Temperature2"                     = "temperature2",
                               "Velocity2"                        = "velocity2")
  
  noise <-list()
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # training model
  training.pipeline(file.input.data=file.input.data,
                    file.prev.taxa=file.prev.taxa,
                    experiment.name=experiment.name,
                    number.split=number.split,
                    split.criterion=split.criterion,
                    number.sample=number.sample,                   
                    models=models,
                    noise=noise,
                    env.factor=env.factor,
                    env.factor.full=env.factor.full)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # plotting models and saving model
  create.plots(experiment.name=experiment.name,
               file.prev.taxa=file.prev.taxa)
  
}


# Noise: removing environmental factor(s) ----
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
for (index in c(1,2,3)){
  
  l <- list(list("velocity"),
            list("velocity",     
                 "A10m",
                 "cow.density", 
                 "IAR"),
            list("velocity",     
                 "A10m",
                 "cow.density", 
                 "IAR",      
                 "urban.area",
                 "FRI",
                 "bFRI",            
                 "width.variability"))
  
  to_remove <- l[[index]]
  
  
  experiment.name         <- paste0("env_fact_removed_", length(to_remove), "_",
                                    format(Sys.time(), "%d_%m_%Y_%Hh%M"))
  number.split            <- 3
  split.criterion         <- "SiteId"
  number.sample           <- 10000
  models                  <- list("null", "glm", "gamloess", "rf", "ann")
  env.factor              <- c("Temperature"                      = "temperature",       # Temp
                               "Flow velocity"                    = "velocity",          # FV
                               "Riparian agriculture"             = "A10m",              # A10m
                               "Livestock unit density"           = "cow.density",       # LUD
                               "Insecticide application rate"     = "IAR",               # IAR
                               "Urban area"                       = "urban.area",        # Urban
                               "Forest-river intersection"        = "FRI",               # FRI
                               "Forest-river intersection buffer" = "bFRI",              # bFRI
                               "Width variability"                = "width.variability") # WV
  
  env.factor.full         <- c(env.factor,
                               "Temperature2"                     = "temperature2",
                               "Velocity2"                        = "velocity2")
  
  
  noise <- lapply(to_remove, FUN=function(target){
    noise.env.fact <- list("type"       = "remove_factor",
                           "target"     = target,
                           "amount"     = NULL, 
                           "parameters" = NULL)
    
    return(noise.env.fact)
  })
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # training model
  training.pipeline(file.input.data=file.input.data,
                    file.prev.taxa=file.prev.taxa,
                    experiment.name=experiment.name,
                    number.split=number.split,
                    split.criterion=split.criterion,
                    number.sample=number.sample,                   
                    models=models,
                    noise=noise,
                    env.factor=env.factor,
                    env.factor.full=env.factor.full)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # plotting models and saving model
  create.plots(experiment.name=experiment.name,
               file.prev.taxa=file.prev.taxa)
  
}


# Noise: adding random environmental factor(s) ----
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
for (nb in c(1,3,9)){
  
  experiment.name         <- paste0("add_env_fact_", nb, "_",
                                    format(Sys.time(), "%d_%m_%Y_%Hh%M"))
  number.split            <- 3
  split.criterion         <- "SiteId"
  number.sample           <- 10000
  models                  <- list("null", "glm", "gamloess", "rf", "ann")
  env.factor              <- c("Temperature"                      = "temperature",       # Temp
                               "Flow velocity"                    = "velocity",          # FV
                               "Riparian agriculture"             = "A10m",              # A10m
                               "Livestock unit density"           = "cow.density",       # LUD
                               "Insecticide application rate"     = "IAR",               # IAR
                               "Urban area"                       = "urban.area",        # Urban
                               "Forest-river intersection"        = "FRI",               # FRI
                               "Forest-river intersection buffer" = "bFRI",              # bFRI
                               "Width variability"                = "width.variability") # WV
  
  env.factor.full         <- c(env.factor,
                               "Temperature2"                     = "temperature2",
                               "Velocity2"                        = "velocity2")
  
  
  noise <- lapply(as.list(seq(1, nb)), FUN=function(n){
    noise.fact <- list("type"       = "add_factor",
                       "target"     = paste0("random",n), # new factor name
                       "amount"     = NULL, 
                       "parameters" = NULL)
    
    return(noise.fact)
  })
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # training model
  training.pipeline(file.input.data=file.input.data,
                    file.prev.taxa=file.prev.taxa,
                    experiment.name=experiment.name,
                    number.split=number.split,
                    split.criterion=split.criterion,
                    number.sample=number.sample,                   
                    models=models,
                    noise=noise,
                    env.factor=env.factor,
                    env.factor.full=env.factor.full)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # plotting models and saving model
  create.plots(experiment.name=experiment.name,
               file.prev.taxa=file.prev.taxa)
  
}



# Noise: adding gaussian noise on an environmental factor (flow velocity)  ----
experiment.name         <- paste0("gaussian_velocity_3",
                                  format(Sys.time(), "%d_%m_%Y_%Hh%M"))
number.split            <- 3
split.criterion         <- "SiteId"
number.sample           <- 10000
models                  <- list("null", "glm", "gamloess", "rf", "ann")
env.factor              <- c("Temperature"                      = "temperature",       # Temp
                             "Flow velocity"                    = "velocity",          # FV
                             "Riparian agriculture"             = "A10m",              # A10m
                             "Livestock unit density"           = "cow.density",       # LUD
                             "Insecticide application rate"     = "IAR",               # IAR
                             "Urban area"                       = "urban.area",        # Urban
                             "Forest-river intersection"        = "FRI",               # FRI
                             "Forest-river intersection buffer" = "bFRI",              # bFRI
                             "Width variability"                = "width.variability") # WV

env.factor.full         <- c(env.factor,
                             "Temperature2"                     = "temperature2",
                             "Velocity2"                        = "velocity2")

noise1 <- list("type"       = "gaussian",
               "target"     = "velocity",
               "amount"     = 3,
               "parameters" = list("min"=-Inf, "max"=Inf))



noise  <- list("noise"     = noise1)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# training model
training.pipeline(file.input.data=file.input.data,
                  file.prev.taxa=file.prev.taxa,
                  experiment.name=experiment.name,
                  number.split=number.split,
                  split.criterion=split.criterion,
                  number.sample=number.sample,                   
                  models=models,
                  noise=noise,
                  env.factor=env.factor,
                  env.factor.full=env.factor.full)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# plotting models and saving model
create.plots(experiment.name=experiment.name,
             file.prev.taxa=file.prev.taxa)





# Noise: adding gaussian noise on two environmental factors (temperature and flow velocity)  ----
experiment.name         <- paste0("gaussian_temp_velocity_3",
                                  format(Sys.time(), "%d_%m_%Y_%Hh%M"))
number.split            <- 3
split.criterion         <- "SiteId"
number.sample           <- 10000
models                  <- list("null", "glm", "gamloess", "rf", "ann")
env.factor              <- c("Temperature"                      = "temperature",       # Temp
                             "Flow velocity"                    = "velocity",          # FV
                             "Riparian agriculture"             = "A10m",              # A10m
                             "Livestock unit density"           = "cow.density",       # LUD
                             "Insecticide application rate"     = "IAR",               # IAR
                             "Urban area"                       = "urban.area",        # Urban
                             "Forest-river intersection"        = "FRI",               # FRI
                             "Forest-river intersection buffer" = "bFRI",              # bFRI
                             "Width variability"                = "width.variability") # WV

env.factor.full         <- c(env.factor,
                             "Temperature2"                     = "temperature2",
                             "Velocity2"                        = "velocity2")

noise1 <- list("type"       = "gaussian",
               "target"     = "temperature",
               "amount"     = 3,
               "parameters" = list("min"=0, "max"=35))

noise2 <- list("type"       = "gaussian",
               "target"     = "velocity",
               "amount"     = 3,
               "parameters" = list("min"=-Inf, "max"=Inf))



noise  <- list("noise1"     = noise1,
               "noise2"     = noise2)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# training model
training.pipeline(file.input.data=file.input.data,
                  file.prev.taxa=file.prev.taxa,
                  experiment.name=experiment.name,
                  number.split=number.split,
                  split.criterion=split.criterion,
                  number.sample=number.sample,                   
                  models=models,
                  noise=noise,
                  env.factor=env.factor,
                  env.factor.full=env.factor.full)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# plotting models and saving model
create.plots(experiment.name=experiment.name,
             file.prev.taxa=file.prev.taxa)
  


