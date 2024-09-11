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
#                                                                              #
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
if (!require("tidyr") ) { install.packages("tidyr"); library("tidyr") }               # to sort, join, merge data
if (!require("ggplot2")){install.packages("ggplot2"); library("ggplot2")}                    # to sort, join, merge data
if (!require("readr")){install.packages("readr"); library("readr")}
if (!require("jsonlite")){install.packages("jsonlite"); library("jsonlite")}
if (!require("pROC")){install.packages("pROC"); library("pROC")}  # to compute AUC    
if (!require("vtable") ) { install.packages("vtable"); library("vtable") }               # to do tables of summary statistics
if (!require("sf") ) { install.packages("sf"); library("sf") }                        # to read layers for plotting maps
if (!require("ggpattern")){install.packages("ggpattern"); library("ggpattern")}        # to add patterns in boxplots

# specific for Neural Network
if (!require("reticulate")){install.packages("reticulate"); library("reticulate")}    # links R to python
# install_miniconda()              # run this the very first time reticulate is installed # makes a separeta python installation
# reticulate::install_python(version = '3.12')

# reticulate::install_python(version = '3.8')
# install.packages("tensorflow")
# remove.packages("tensorflow")
# virtualenv_list()
# virtualenv_exists(envname = NULL)
# virtualenv_create("r-tensorflow", python = "python3.8")  # or python3.9
# conda_create("r-tensorflow", python_version = "3.8")
library("tensorflow")
# use_condaenv("r-tensorflow", required = TRUE)
# tensorflow::tf_config()
# conda_list()
# virtualenv_list()
# virtualenv_remove("r-tensorflow") 
# install_tensorflow(method = "conda", envname = "r-tensorflow")
# path <-"C:/Users/cholleem/AppData/Local/r-miniconda/envs/r-reticulate/python.exe"
# tensorflow::tf_config() # Check if TensorFlow is installed and accessible
# virtualenv_install("C:/Users/ClientAdmin/Documents/.virtualenvs/r-tensorflow", "tensorflow==2.16")
# install_tensorflow(envname = "C:/Users/ClientAdmin/Documents/.virtualenvs/r-tensorflow")
# virtualenv_remove("r-tensorflow")
# install_tensorflow(method = "conda")             # run this line only when opening new R session
# install_tensorflow()
# install_tensorflow(envname = "C:/Users/ClientAdmin/Documents/.virtualenvs/r-tensorflow")
# install.packages("keras")
# remove.packages("keras")
library("keras")
install_keras()
# install_keras(method = "conda")                  # run this line only when opening new R session
# install_tensorflow(method = "conda", envname = "r-tensorflow")
# use_condaenv()
# reticulate::install_python(version = '<version>')
# Sys.setenv(RETICULATE_PYTHON = path)
# virtualenv_create("r-tensorflow")

# caret has to be loaded at the end to not cache function 'train'
if (!require("mgcv")){install.packages("mgcv"); library("mgcv")}
if (!require("caret")){install.packages("caret"); library("caret")}


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Directories and files ####

# define directories
dir.input.data          <- "../A3_outputs_streambugs/"
dir.output              <- "../B2_outputs_sdms/"
dir.utilities           <- "../00_utilities/"
# dir.create(dir.output)

# define files
name.streambugs.run     <- "8Catch_3009Sites_PolyInterp30_runC_100Yea_36501Steps"
file.input.data         <- paste0(name.streambugs.run, "_WideData_ResultsSamplePresAbs.csv")
file.prev.taxa          <- paste0(name.streambugs.run, "_PrevalenceTaxonomy_SamplePresAbs.csv")
file.selected.taxa      <- "selected_taxa_analysis.csv"


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Load data and functions ####

# read data
data.env.taxa       <- read.csv(paste0(dir.input.data, file.input.data), header = T, sep = ";", stringsAsFactors = F)
data.prev.taxa      <- read.csv(paste0(dir.input.data, file.prev.taxa), header = T, sep = ";", stringsAsFactors = F)
data.selected.taxa  <- read.csv(paste0(dir.utilities, file.selected.taxa), header = T, sep = ";", stringsAsFactors = F)

# load functions
source(paste0(dir.utilities, "utilities_global.r"))
source("data_noising.r")
source("performances_assessment.r")
source("training_pipeline.r")
source("ml_models.r")
source("global_variables.r")
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
# env.factor <- vect.dir.env.fact
env.factor <- c(vect.dir.env.fact, vect.indir.env.fact)
env.factor.full <- c(env.factor, "Temp2" = "tempmaxC2", "FV2" = "currentms2")
no.env.fact <- length(env.factor)

# taxa list full
cind.taxa <- which(grepl("Occurrence.", colnames(data.env.taxa)))
vect.taxa.full <- colnames(data.env.taxa)[cind.taxa]
names(vect.taxa.full) <- gsub("Occurrence.", "", vect.taxa.full)
no.taxa.full <- length(vect.taxa.full)
print(vect.taxa.full)

# update number of NAs in input data
for (taxon in vect.taxa.full) {
    # taxon <- vect.taxa.full[1]
    data.prev.taxa[which(data.prev.taxa$Occurrence.taxa == taxon), "Missing.values"] <- sum(is.na(data.env.taxa[,taxon]))
}
summary(data.prev.taxa$Missing.values)

# taxa list selected for analysis and above prevalence threshold
prev.threshold <- 0.1
na.threshold <- 1000
temp.occ.taxa <- data.prev.taxa[which(data.prev.taxa$Prevalence.NoDisp > prev.threshold
                                      & data.prev.taxa$Prevalence.NoDisp < 1 -prev.threshold
                                      & data.prev.taxa$Missing.values < na.threshold), "Occurrence.taxa"]

selected.taxa.analysis <- data.selected.taxa$Occurrence.taxa
# taxa.colnames <- selected.taxa.analysis[which(selected.taxa.analysis %in% temp.occ.taxa)]
taxa.colnames <- temp.occ.taxa

names(taxa.colnames) <- gsub("Occurrence.", "", taxa.colnames)
no.taxa <- length(taxa.colnames)
print(no.taxa)
print(taxa.colnames)

# select colors for present/absent plots
vect.col.pres.abs <- c( "Present" = "#00BFC4", "Absent" = "#F8766D")
scales::show_col(vect.col.pres.abs)

# sdms training options
number.split            <- 3
split.criterion         <- "ReachID"
number.sample           <- 3000
number.sample           <- ifelse(dim(data.env.taxa)[1] < number.sample, dim(data.env.taxa)[1], number.sample)
sdm.models              <- c(
                                # "GLM" = "glm",
                                # # "GAM" = "gam",
                                # "GAM" = "gamSpline",
                                # "RF"  = "rf",
                                "ANN" = "ann")#,
# "ann")
no.models <- length(sdm.models)
models <- append(c("Null" = "null"), sdm.models)

# noise scenarios
list.noise <- list(
  
  # "noise.disp"     = list("type"       = "dispersal",
  #                        "target"     = NULL,
  #                        "amount"     = NULL,
  #                        "parameters" = NULL) # ,

  # "noise.temp"     = list("type"       = "gaussian",
  #                         "target"     = "tempmaxC",
  #                         "amount"     = 3,
  #                         "parameters" = list("min"=0, "max"=30)) #,
  # # 
  # "noise.gamm"     = list("type"       = "missdetection",
  #                         "target"     = "Occurrence.Gammaridae",
  #                         "amount"     = 0.1,
  #                         "parameters" = NULL)#,
  
  # # 
  # "noise.add.fact" = list("type"       = "add_factor",
  #                         "target"     = "random1", # new factor name
  #                         "amount"     = NULL,
  #                         "parameters" = NULL),
  # 
  # "noise.rem.fact" = list("type"       = "remove_factor",
  #                         "target"     = "currentms",
  #                         "amount"     = NULL,
  #                         "parameters" = NULL)#,
  
)
name.list.noise <- paste(names(list.noise), collapse = "_")

# p <- 0.1 # probability at which a presence is turned into an absence
# 
# list.noise <- lapply(taxa.colnames, FUN=function(taxon){
#     noise_taxon <- list("type"       = "missdetection",
#                         "target"     = taxon,
#                         "amount"     = p,
#                         "parameters" = NULL)
# 
#     return(noise_taxon)
# })
# name.list.noise <- paste0("misdetection.all.taxa", p)
  

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# PREPROCESS DATA ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Add noise ####

# preprocess input dataset
data.input <- data.env.taxa
data.input$tempmaxC2 <- data.input$tempmaxC^2
data.input$currentms2 <- data.input$currentms^2
data.input <- data.input[,c(vect.info,env.factor.full,taxa.colnames)] # subselect columns of interest
catch.variable <- "Watershed"

# add noise
data.input.noised        <- add.noise(data = data.input,
                                      number.sample,
                                      noise = list.noise,
                                      env.fact = env.factor,
                                      env.fact.full = env.factor.full)
data.input      <- data.input.noised$`noised data`
env.factor      <- data.input.noised$`env fact`
env.factor.full <- data.input.noised$`env fact full`
no.env.fact <- length(env.factor)

# reduce dataset to selected number of samples
rind.select <- runif(n=number.sample, min=1, max=dim(data.input)[1])
data.input <- data.input[rind.select,]

# sanity check of NAs
sum(is.na(data.input[,c(vect.info,env.factor)])) # should be zero 
sum(is.na(data.input[,c(taxa.colnames)])) # could have a lot, decide to replace by zero or not

# replace NAs by 0 to simulate noise due to dispersal limitation
if("noise.disp" %in% names(list.noise)){
  data.input[is.na(data.input)] <- 0 # replace NAs by 0
}

# replace 0/1 by absent/present
for (taxon in taxa.colnames) {
  data.input[which(data.input[,taxon] == 0),taxon] <- "Absent"
  data.input[which(data.input[,taxon] == 1),taxon] <- "Present"
  data.input[,taxon] = as.factor(data.input[,taxon])
}

# # make long dataframe for analysis
# long.df.input.taxa <- gather(data.input, key = Taxon, value = Observation, all_of(taxa.colnames))
# long.df.input.taxa$Taxon <- gsub("Occurrence.", "", long.df.input.taxa$Taxon)
# long.df.input.env.taxa <- gather(long.df.input.taxa, key = Factor, value = Value, all_of(env.factor))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Write metadata ####

# define experiment name
experiment.name <- paste0(
    number.sample, "Sites_",
    no.taxa, "Taxa_",
    no.models, "SDMs_",
    no.env.fact, "EnvFact_",
    name.list.noise
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
## Prevalence analysis ####

# # plot bars pres/abs/NA
# plot.data <- long.df.input.taxa
# plot.data$Observation <- factor(plot.data$Observation, levels=c("Present", NA, "Absent"))
# temp.taxa <- unique(plot.data$Taxon)
# temp.catch <- unique(plot.data$Watershed)
# 
# # analyze distribution of prevalence and NAs across catchments
# p <- ggplot(plot.data, aes(x = Taxon, fill = Observation))
# p <- p + geom_bar(stat="count", position='stack')
# p <- p + facet_wrap(~ Watershed, scales = "free_y")
# p <- p + scale_fill_manual(values=vect.col.pres.abs)
# p <- p + theme_bw()
# p <- p + scale_x_discrete(limits = temp.taxa)
# p <- p + theme(axis.text.x = element_text(angle=90))
# p
# 
# # analyze distribution of prevalence and NAs across taxa
# q <- ggplot(plot.data, aes(x = Watershed, fill = Observation))
# q <- q + geom_bar(stat = "count", position = "stack")
# q <- q + facet_wrap(~ Taxon, scales = "free_y")
# q <- q + scale_fill_manual(values=vect.col.pres.abs)
# q <- q + theme_bw()
# q <- q + scale_x_discrete(limits = temp.catch)
# q <- q + theme(axis.text.x = element_text(angle=90))
# q
# 
# # print plots in pdf
# file.name <- "GeomBar_PresAbs_TaxaCatch"
# print.pdf.plots(list.plots = list(q,p), width = 12, height = 9, dir.output = dir.experiment, info.file.name = "", file.name = file.name, 
#                 png = F)
# 
# plot.data <- filter(plot.data, Taxon %in% "Gammaridae")
# q <- ggplot(plot.data, aes(x = Watershed, fill = Observation))
# q <- q + geom_bar(stat = "count", position = "stack")
# q <- q + scale_fill_manual(values=vect.col.pres.abs)
# q <- q + theme_bw()
# q <- q + scale_x_discrete(limits = temp.catch)
# q <- q + theme(axis.text.x = element_text(angle=90))
# q
# file.name <- "GeomBar_PresAbs_GammaridaeCatch"
# print.pdf.plots(list.plots = list(q), width = 7, height = 5, dir.output = dir.experiment, info.file.name = "", file.name = file.name, 
#                 png = F)

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

# recover Streambugs data for ICE
data.base.ice       <- read.csv(paste0(dir.input.data, "WideDataInputs_ICE_50Sites.csv"), sep = ";")
data.base.ice$tempmaxC2 <- data.base.ice$tempmaxC^2
data.base.ice$currentms2 <- data.base.ice$currentms^2
name.ice.streambugs <- "8Catch_50Sites_ICE_50Steps_PolyInterp30_runC_100Yea_36501Steps_"
ice.df.streambugs   <- read.csv(paste0(dir.input.data, name.ice.streambugs,
                                       "WideData_ResultsProbObs.csv"), sep = ";")
no.sites <- as.numeric(stringr::str_match(name.ice.streambugs, "Catch_(.+)Sites")[2])
no.steps <- as.numeric(stringr::str_match(name.ice.streambugs, "ICE_(.+)Steps_Poly")[2])

# check balance presence/absence
for(i.split in 1:number.split){
    cat("Split number:", i.split, "\n")
    for (taxon in taxa.colnames) {
        cat(taxon)
        print(summary(preprocessed.data.cv[[i.split]][["Training data"]][,taxon]))
    }
}
print(summary(preprocessed.data.cv[[1]][["Training data"]][,"Occurrence.Baetisalpinus"]))
# Absent Present 
# 1149     836

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# APPLY MODELS ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Cross-validation ####

models.cv  <- apply.ml.models(data=preprocessed.data.cv,
                              models=models,
                              split.type="CV", 
                              taxa.colnames = taxa.colnames, 
                              env.factor = env.factor, 
                              env.factor.full = env.factor.full)

save.models(models=models.cv,
            path=dir.experiment,
            split.type="CV")

# models.cv <- lapply(sdm.models, FUN = apply.caret.model,
#                     preprocessed.data.cv, split.type = "CV",
#                     taxa.colnames, env.factor, list.noise)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Fit entire dataset ####

models.fit <- apply.ml.models(data=preprocessed.data.fit,
                              models=models,
                              split.type="FIT",
                              taxa.colnames, env.factor, env.factor.full)

save.models(models=models.fit,
            path=dir.experiment,
            split.type="FIT")


# models.fit <- lapply(sdm.models, FUN = apply.caret.model, 
#                      preprocessed.data.fit, split.type = "FIT", 
#                      taxa.colnames, env.factor, list.noise)


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
prev.taxa$Prevalence.DispLim <- prev.taxa$Prevalence
prev.taxa$Prevalence <- prev.taxa$Prevalence.NoDisp

## Variable for plots ----
all.results <- summarize.all.results(models.cv, prev.taxa)

# write summarized results in csv
file.name <- paste0(experiment.name, "allresults.csv")
write.table(all.results, file = paste0(dir.experiment, file.name), sep = ";", row.names = F)

ggplot.all.results <- restructure.all.results(all.results, models)
ggplot.all.results$model <- factor(ggplot.all.results$model, levels= names(models))

color.map <- c('GLM'           = 'deepskyblue',   # Generalized Linear Model
               'GAM'           = 'green',         # Generalized Additive Model
               'ANN'           = 'orange',        # Artificial Neural Network
               'RF'            = 'red',           # Random Forest
               'Null'          = 'black')         # Null Model
scales::show_col(color.map)

## Fig. 1 : boxplot perf ----
gg.plot.dev.infos <-  ggplot.all.results %>%
  group_by(fit_pred, model) %>%
  summarise(mean=mean(dev),
            max=max(dev),
            min=min(dev),
            median=median(dev))

median.null.pred <- min((gg.plot.dev.infos %>%
                           filter(fit_pred=="pred" & model=="Null"))["median"])
median.best.pred <- min((gg.plot.dev.infos %>%
                           filter(fit_pred=="pred"))["median"])

fig1 <- ggplot(data=ggplot.all.results,
               aes(x=model, y=dev)) +
    geom_boxplot(aes(fill=fit_pred)) +
    geom_hline(aes(yintercept=median.null.pred,
                 color="median null model (pred)"),
             linetype="dashed") +
    geom_hline(aes(yintercept=median.best.pred,
                 color="median best model (pred)"),
             linetype="dashed") +
    scale_fill_manual(values=c(fit = "#998ec3", pred = "#f1a340")) +
    scale_colour_manual(values = c('green','black')) +
    theme_minimal() +
    theme(legend.title=element_blank()) +
    labs(x = "Models",
         y = "Standardized deviance\n ",
         title = "")

pdf(paste0(dir.experiment, "boxplot_performance.pdf"), width = 8, height = 5)
print(fig1)
dev.off()



## Fig. 2 : bellplot ----
fig2 <- ggplot(data=subset(ggplot.all.results, model %in% "Null"),
               aes(x=prevalence, y=dev, color=model)) +
  geom_smooth(se = FALSE) + 
  geom_point(data=ggplot.all.results) +
  facet_wrap(~fit_pred) +
  scale_color_manual(values=color.map) +
  theme_minimal() + 
  theme(legend.title=element_blank()) +
    labs(x = "Prevalence",
         y = "Standardized deviance\n ",
         title = "")


pdf(paste0(dir.experiment, "bellplot_performance.pdf"), width = 8, height = 5)
print(fig2)
dev.off() 

## Fig. 3 : ICE ----


dir.exp.ice <- paste0(dir.experiment, "ICE/")
dir.create(dir.exp.ice)
preprocessed.data.ice <- preprocess.data(data=data.base.ice,
                                         env.fact=env.factor,
                                         dir=dir.exp.ice,
                                         split.type="FIT")
std.const.ice   <- readRDS(file=paste0(dir.exp.ice, "standardization_constant_FIT.rds"))


# select.env.fact <- "temperature"
select.env.fact <- "tempmaxC"

list.plots <- list()

for(taxon.under.obs in names(taxa.colnames)){
  
  cat("Plotting ICE for:", taxon.under.obs)

  # taxon.under.obs <- names(taxa.colnames)[1]
  
  pred.ice.streamb <- ice.df.streambugs[,c("ReachID", "Watershed", select.env.fact, paste0("Occurrence.", taxon.under.obs))] %>%
    mutate(observation_number = rep(1:no.sites, each = no.steps),
           column_label = 1,
           model = "Streambugs") %>%
    rename(pred = paste0("Occurrence.", taxon.under.obs))
  
  env.factor.sampled <- data.frame(variable = data.base.ice[,select.env.fact])
  
  pred.ice.mean <- pred.ice.streamb %>%
    group_by_at(select.env.fact[1]) %>%
    summarize(avg = mean(na.omit(pred))) %>%
    mutate(model = "Streambugs")
  pred.ice.mean.bounds  <- pred.ice.mean %>%
    # group_by(model) %>%
    summarize(x.mean=max(select.env.fact),
              y.mean.min=min(avg),
              y.mean.max=max(avg)) %>%
    mutate(model = "Streambugs")

  # input.env.factors <- get.input.env.factors(experiment.name)
  input.env.factors <- env.factor
  
  ice.dfs <- plot.ice(models.performance = models.fit,
                      select.env.fact=select.env.fact,
                      taxa=taxon.under.obs,
                      standardization.constant=std.const.ice[[1]],
                      observations = preprocessed.data.ice[["Entire dataset"]][,input.env.factors],
                      nb.sample=no.sites,
                      resolution=no.steps,
                      input.env.factors=input.env.factors)
  
  pred.models              <- ice.dfs[["observations"]]
  env.factor.sampled       <- ice.dfs[["env.factor.sampled"]]
  pred.models.mean         <- pred.models %>%
    group_by(across(all_of(select.env.fact)), model) %>%
    summarise(avg = mean(pred))
  pred.models.mean.bounds  <- pred.models.mean %>% group_by(model) %>%
    summarise(x.mean=max(across(all_of(select.env.fact))),
              y.mean.min=min(avg),
              y.mean.max=max(avg))
  # colnames(pred.ice.streamb)
  # colnames(pred.models)
  col.names <- c(select.env.fact, "pred", "observation_number", "column_label", "model")
  plot.data <- rbind(pred.ice.streamb[,col.names], pred.models[,col.names] )
  plot.data.mean  <- rbind(pred.ice.mean, pred.models.mean)       
  plot.data.mean.bounds  <- rbind(pred.ice.mean.bounds, pred.models.mean.bounds)
  # plot.data$model <- factor(plot.data$model, levels=c("GLM", "GAM", "RF"))
  
  fig3 <- ggplot(data=plot.data) +
    geom_line(aes(x=.data[[select.env.fact]],
                  y=pred,
                  group=observation_number, 
                  color=as.character(observation_number)),
              show.legend = FALSE) +
    geom_line(data=plot.data.mean,
              aes(x=.data[[select.env.fact]], y=avg),
              size=1.5) +
    geom_rug(data = env.factor.sampled,
             aes(x=variable), 
             color="grey20",
             alpha=0.7,
             inherit.aes=F) + 
    scale_x_continuous(limits = c(min(env.factor.sampled), max(env.factor.sampled))) +
    # geom_segment(data=plot.data.mean.bounds,
    #              inherit.aes = FALSE,
    #              lineend="round",
    #              linejoin="round",
    #              aes(x=x.mean,
    #                  y=y.mean.min,
    #                  xend=x.mean,
    #                  yend=y.mean.max),
    #              arrow=arrow(length = unit(0.3, "cm"),
    #                          ends = "both")) +
    facet_wrap(~factor(model, levels = c("Streambugs", names(sdm.models))), ncol = 4) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          legend.title = element_text(size=24),
          legend.text = element_text(size=20)) +
    labs(title = taxon.under.obs,
         x = "Temperature",
         y = "Predicted probability of occurrence")
  # fig3
  
  list.plots[[taxon.under.obs]] <- fig3
  
  # pdf(paste0(dir.experiment, "ice_", experiment.name, taxon.under.obs, ".pdf"))
  # print(fig3)
  # dev.off()
}  

file.name <- paste0("ice_alltaxa_4x1")
print.pdf.plots(list.plots = list.plots, width = 10, height = 5, 
                dir.output = dir.experiment, info.file.name = experiment.name, file.name = file.name, 
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
## Options for plots comparison ----

# create directory for comparison plots
dir.compar.plots <- paste0(dir.output, "comparison_plots/")
dir.create(dir.compar.plots)

# set color map for models
model.color.map <- c('GLM'     = "#619CFF",  # 'deepskyblue',   # Generalized Linear Model
                     'GAM'     = "#00BA38", # 'green',         # Generalized Additive Model
                     'ANN'     = 'orange',        # Artificial Neural Network
                     'RF'      = "#F8766D") # 'red')#,           # Random Forest
                     # 'Null'          = 'black')         # Null Model
scales::show_col(model.color.map)

# select taxon and env. factor
taxon.under.obs <- names(taxa.colnames)[5]
select.env.fact <- env.factor[1]
name.select.env.fact <- names(select.env.fact)
print(taxon.under.obs)
print(select.env.fact)

# select experiments to compare
list.exp     <- list("Best case scenario"                  = "3000Sites_45Taxa_4SDMs_8EnvFact_",
                     "Reduce dataset size"                 = "300Sites_45Taxa_4SDMs_8EnvFact_",
                     "Remove environmental predictor"      = "3000Sites_45Taxa_4SDMs_7EnvFact_noise.rem.fact",
                     "Noise on temperature"                = "3000Sites_45Taxa_4SDMs_8EnvFact_noise.temp",
                     "Misdetection"                        = "3000Sites_45Taxa_4SDMs_8EnvFact_misdetection.all.taxa0.1",
                     "Noise dispersal limitation"          = "3000Sites_45Taxa_4SDMs_8EnvFact_noise.disp")

# test.colors <- RColorBrewer::brewer.pal(8, "Set1")
# print(test.colors)
# scales::show_col(test.colors)

# color.map <- RColorBrewer::brewer.pal(length(list.exp), "Dark2")
# scales::show_col(color.map)

color.map <- c('1'            = '#4daf4a',
               '2'            = '#377eb8',
               '3'            = '#984ea3',
               '4'            = '#e41a1c',
               '5'            = '#ff7f00',
               '6'            = '#ffd92f')

color.map <- color.map[1:length(list.exp)]
names(color.map) <- names(list.exp)
scales::show_col(color.map)

# set options for plots
# width and height of an A4 page in inches
width.a4 = 8.3
height.a4 = 11.7

mytheme <- theme_bw() +
    theme(text = element_text(size = 14),
          title = element_text(size = 14),
          # plot.title = element_text(face = "bold"),
          # strip.text = element_text(size = 14), #, face = "bold"),
          strip.background = element_rect(fill = "white"),
          # axis.title = element_text(size = 14), #, face = "bold"),
          # axis.text.x = element_text(vjust = 10),
          # axis.text.y = element_text(hjust = 10),
          # axis.title.x.bottom = element_text(vjust = -15),
          # axis.title.x = element_text(margin = margin(t = 0, r = 0, b = -30, l = 0)), # put some space between the axis title and the numbers
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
          legend.position="bottom"
          # legend.text = element_text(size = 14),
          # legend.title = element_text(size = 14) #, face= "bold") #,
          # panel.grid.major = element_line(colour = NA),
          # panel.grid.minor = element_line(colour = NA)
    )

theme_set(mytheme)

# create file names for saving plots
file.name.exp <- paste0("comparison_300points_rem.FV_noise.temp_misdet_disp.lim_", 
                        length(list.exp), "exp_")
print(file.name.exp)
file.name.tax <- paste0(file.name.exp,
                        taxon.under.obs, "_",
                        select.env.fact)
print(file.name.tax)

# create.comparison.plots("gauss_temp", list.exp.gauss, 
#                         file.prev.taxa="8Catch_1416Sites_curve.curr.temp-10_interc.orgmic.sapro4_10Yea_3651Steps_PrevalenceTaxonomy_ThreshPresAbs.csv",
#                         taxon="Gammaridae",
#                         env.fact="tempmaxC")


## ICE and PDP data ----


multi.ice <- lapply(list.exp, FUN=function(name){
  
  # name <- list.exp[[1]]
  dir.experiment          <- paste0(dir.output, name, "/")
  cat("\nLoading models and computing ICE for taxon", taxon.under.obs, "and experiment:", name, "\n")
  
  models.fit      <- load.models(path=dir.experiment, split.type="FIT")
  std.const.fit   <- readRDS(file=paste0(dir.experiment, "standardization_constant_FIT.rds"))
  correct.names <- names(models)[which(names(models.fit) %in% models | names(models.fit) %in% names(models))]
  names(models.fit) <- correct.names
  file.metadata <- "metadata.json"
  metadata.exp <- fromJSON(paste(readLines(paste0(dir.experiment, file.metadata)), collapse=""))
  
  input.env.factors <- metadata.exp$env_factor
  
  ice.dfs <- plot.ice(models.performance=models.fit,
                      select.env.fact=select.env.fact,
                      taxa=taxon.under.obs,
                      standardization.constant=std.const.fit[[1]],
                      observations = preprocessed.data.fit$`Entire dataset`,
                      nb.sample=no.sites,
                      resolution=no.steps,
                      input.env.factors=input.env.factors)
  
  observations              <- ice.dfs[["observations"]]
  
  observations.mean         <- observations %>%
    group_by(across(all_of(select.env.fact)), model) %>%
    summarise(avg = mean(pred))
  
  observations.mean["noise"] <- name
  observations["noise"]      <- name
  
  return(list(observations, observations.mean))
})

long.multi.ice <- bind_rows(lapply(multi.ice, "[[", 1), .id = "column_label_noise")
long.multi.ice$model <- factor(long.multi.ice$model, levels = names(sdm.models))
long.multi.ice$column_label_noise <- factor(long.multi.ice$column_label_noise, levels = unlist(names(list.exp)))

multi.mean <- lapply(multi.ice, "[[", 2)
final.multi.ice <- bind_rows(multi.mean, .id = "column_label")
final.multi.ice$model <- factor(final.multi.ice$model, levels = names(sdm.models))
final.multi.ice$column_label <- factor(final.multi.ice$column_label, levels = unlist(names(list.exp)))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fig 0: ICE comparison

plot.data <- long.multi.ice


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fig 1 & Fig 2: PDP comparison


# compute boundaries
min.boundaries <- lapply(multi.mean, FUN=function(ice){
  min(ice[[name.select.env.fact]])
})

max.boundaries <- lapply(multi.mean, FUN=function(ice){
  max(ice[[name.select.env.fact]])
})

lb <- max(unlist(min.boundaries))
hb <- min(unlist(max.boundaries))

plot.data <- final.multi.ice

pred.streambugs.mean <- ice.df.streambugs[,c("ReachID","Watershed", select.env.fact, paste0("Occurrence.", taxon.under.obs))] %>%
  mutate(observation_number = rep(1:no.sites, each = no.steps),
         column_label = 1) %>%
  rename(pred = paste0("Occurrence.", taxon.under.obs)) %>%
  group_by_at(select.env.fact[1]) %>%
  summarize(avg = mean(na.omit(pred)))

fig1 <- ggplot(data=plot.data) +
  geom_line(aes(x=.data[[name.select.env.fact]],
                y=avg,
                group=model,
                colour=model),
            size=1, alpha = 0.8) +
  geom_line(data=pred.streambugs.mean,
            aes(x=.data[[name.select.env.fact]], y=avg),
            size=1, color = "darkgrey", alpha = 0.8) +
  facet_wrap(~factor(column_label, levels=unlist(names(list.exp))), ncol = 3)+
  #facet_wrap(~column_label) +
  xlim(lb, hb) +
  scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values=model.color.map) +
    # theme(axis.text.x = element_text(hjust = -1)) +
  labs(x = name.select.env.fact,
       y = "Predicted probability of occurrence",
       colour="Models",
       title = taxon.under.obs)
# fig1

pdf(paste0(dir.compar.plots, file.name.tax, "_pdp_per_scenario_2x3.pdf"), width = width.a4, height = height.a4/2)
print(fig1)
dev.off()

# rm(fig1)

fig2 <- ggplot(data=plot.data) +
  geom_line(aes(x=.data[[name.select.env.fact]],
                y=avg,
                group=column_label,
                colour=column_label),
            size=1) +
  geom_line(data=pred.streambugs.mean,
            aes(x=.data[[name.select.env.fact]], y=avg),
            size=1, color = "darkgrey", alpha = 0.8) +
  scale_color_manual(values=color.map) +
  xlim(lb, hb) +
  # scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~model) +
  # facet_wrap(~factor(model, levels = c("GLM", "GAM", "RF"))) +
  # theme_bw() +
  # theme(strip.background = element_rect(fill = "white"),
  #       legend.position = "bottom") +#,
  # legend.title = element_text(size=24),
  # legend.text = element_text(size=20)) +
  labs(x = name.select.env.fact,
       y = "Predicted probability of occurrence",
       colour="Scenario", 
       title = taxon.under.obs)
# fig2

pdf(paste0(dir.compar.plots, file.name.tax, "_pdp_per_model.pdf"), width = width.a4, height = height.a4*2/3)
print(fig2)
dev.off()

# rm(fig2)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fig 3: plotting the score of given taxon

# load and process all results
multi.all.results <- lapply(list.exp, FUN=function(name){
  
    # name <- list.exp[[2]]
  dir.experiment <- paste0(dir.output, name, "/")
  # models.cv      <- load.models(path=dir.experiment, split.type="CV")
  file.name <- paste0(dir.experiment, name, "allresults.csv")
  
  all.results <- read.table(file.name, header = T, sep = ";")
  for(i in 1:ncol(all.results)){
      # i <- 1
      correct.name <- gsub(models[1], names(models)[1],
                           gsub(models[2], names(models)[2],
                                gsub(models[3], names(models)[3],
                                     gsub(models[4], names(models)[4],
                                     colnames(all.results)[i]))))
      colnames(all.results)[i] <- correct.name
  }
  # all.results    <- summarize.all.results(models.cv, data.prev.taxa)
  # rm(models.cv) 
  
  restr.all.results <- restructure.all.results(all.results, models)
  restr.all.results["noise"] <- name
  
  return(list(all.results, restr.all.results))
})

select.results <- c("dev_pred", "auc_pred", "likelihood_ratio")
summary.metric <- "median"
table.summary.metrics <- bind_rows(lapply(names(list.exp), FUN = function(name, multi.all.results, select.results, summary.metric){

    # name <- names(list.exp)[1]
    temp.res <- multi.all.results[[name]][[1]]
    temp.names <- apply(expand.grid(names(models), select.results), 1, paste, collapse="_")
    temp.summary <- temp.res[, temp.names] %>%
        summarise_if(is.numeric, summary.metric, na.rm = TRUE) %>%
        mutate(scenario = name) %>% 
        mutate(across(is.numeric, round, digits=2))
    return(temp.summary)
    
    }, multi.all.results, select.results, summary.metric), .id = "id")

# write summarized results in csv
file.name <- paste0(file.name.exp, summary.metric, "_results.csv")
write.table(table.summary.metrics, file = paste0(dir.compar.plots, file.name), sep = ",", row.names = F)

final.multi.all.results <- bind_rows(lapply(multi.all.results, "[[", 2), .id = "column_label")
final.multi.all.results$model <- factor(final.multi.all.results$model, levels= names(models))
final.multi.all.results$column_label <- factor(final.multi.all.results$column_label, levels = names(list.exp))

filtered.multi.all.results <- final.multi.all.results %>%
  filter(taxa == taxon.under.obs)

fig3 <- ggplot(data=filtered.multi.all.results) +
  geom_point(aes(x=column_label,
                 y=dev,
                 shape=fit_pred, 
                 color=model),
             size=3,
             alpha=0.7) +
  scale_x_discrete(limits=names(list.exp)) +
  scale_color_manual(values=model.color.map) +
  labs(x="Scenario",
       y="Standardized deviance",
       title = taxon.under.obs,
       color = "Model",
       shape = "")
  #facet_wrap(~fit_pred) + 
  # theme_minimal() +
  # theme(legend.title=element_blank(),
  #       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


pdf(paste0(dir.compar.plots, file.name.tax, "_taxa_score.pdf"), width = width.a4*1.5, height = height.a4/2)
print(fig3)
dev.off()

rm(fig3)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fig 6: multi-scenario bellplots ----

plot.data <- final.multi.all.results %>%
    filter(fit_pred == "pred")

fig6 <- ggplot(data=plot.data,
               aes(x=prevalence, y=dev, color=column_label)) + 
    geom_point() +
    facet_wrap(~model)  +
    scale_color_manual(values=color.map) + 
    scale_y_continuous(limits = c(0,2)) +
    labs(x = "Prevalence",
         y = "Standardized deviance") +
    # scale_x_discrete(names(taxa.colnames)) +
    theme_minimal() + 
    theme(legend.title=element_blank())
# fig6

pdf(paste0(dir.compar.plots, file.name.exp, "_bellplot_pred_per_scenario.pdf"), height = 6, width = 10)
print(fig6)
dev.off() 

fig7 <- ggplot(data=plot.data,
               aes(x=prevalence, y=dev, color=model)) + 
    geom_point() +
    facet_wrap(~column_label)  +
    scale_color_manual(values=model.color.map) + 
    scale_y_continuous(limits = c(0,2)) +
    labs(x = "Prevalence",
         y = "Standardized deviance",
         color = "Model")# +
    # scale_x_discrete(names(taxa.colnames)) +
    # theme_minimal() + 
    # theme(legend.title=element_blank())

# fig7

pdf(paste0(dir.compar.plots, file.name.exp, "_bellplot_pred_per_model.pdf"), height = height.a4/2, width = width.a4)
print(fig7)
dev.off() 

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fig 5: multi box plot ----

plot.data <- final.multi.all.results
# plot.data$model <- factor(plot.data$model, levels= names(models))
plot.data$pattern <- ifelse(plot.data$fit_pred == "fit", "Training", "Testing")
plot.data$pattern <- factor(plot.data$pattern, levels= c("Training", "Testing"))


fig5 <- ggplot(data=plot.data, aes(x=model, y=dev)) +
    # geom_boxplot(aes(x=column_label,
    #                  y=dev,
    #                  fill=model,
    #                  alpha = pattern)) +
    geom_boxplot_pattern(aes(fill = model, pattern = pattern),
                         pattern_colour = 'black',
                         pattern_fill = 'black',
                         pattern_density = 0.1,
                         # pattern_spacing = 0.01,
                         pattern_key_scale_factor = 1.5
                         # pattern = 'stripe', 
                         # pattern_size = 0.2
                         ) +
    # geom_boxplot_pattern(aes(fill = model, pattern = pattern, pattern_fill = pattern)) +
                         # , pattern_density = 0.35, outlier.shape = NA) +
    # scale_x_discrete(limits=rev(names(models))) +
    scale_x_discrete(limits=names(models)) +
    scale_y_continuous(limits = c(0, 2)) +
    facet_wrap(~column_label, ncol = 3) +
               # , strip.position="left") +
    # theme_minimal() +
    # scale_alpha_manual(values = c("fit" = 0.3, "pred" = 1)) +
    scale_pattern_manual(values= c("Training" = "stripe", "Testing" = "none")) + #,
                         # guide = guide_legend(reverse = TRUE)) + # manually assign pattern
    # scale_pattern_continuous(limits = c("fit" = "Training", "pred" = "Testing")) +
    scale_fill_manual(values = c("Null" = "grey", model.color.map)) +
    # theme(legend.title=element_blank()) +
          # axis.text.y = element_blank()) + #,
          # axis.text.x = element_text(angle = 15, # vjust = 0.5, 
          #                            hjust=1)) +
    labs(x = "Model",
         y = "Standardized deviance",
         fill = "",
         pattern = "")
fig5


pdf(paste0(dir.compar.plots, file.name.exp, "boxplot_colormodels_stripesfit_2x3scenarios.pdf"), width = width.a4*1.2, height = height.a4*2/3)
print(fig5)
dev.off() 

# rm(fig5)

fig6 <- ggplot(data=plot.data) +
    geom_boxplot(aes(x=model,
                     y=dev,
                     fill=column_label)) +
    scale_x_discrete(limits=names(models))+
    scale_y_continuous(limits = c(0, 2)) +
    facet_wrap(~fit_pred, nrow = 2) + 
    theme_minimal() +
    scale_fill_manual(values = color.map) +
    theme(legend.title=element_blank()) +
    labs(x="Model",
         y="Standardized deviance")
# fig6


pdf(paste0(dir.compar.plots, file.name.general, "_boxplot_per_model.pdf"), width = 12, height = 8)
print(fig6)
dev.off() 


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fig 4: dual ICE ----
# exp.baseline <- list.exp[[1]]
# # exp.extreme  <- list.exp[[length(list.exp)]]
# no.select.exp <- 3
# exp.extreme  <- list.exp[[no.select.exp]]

list.dual.exp <- list.exp[c(1,3)]

# taxon.under.obs <- paste0("Occurrence.", taxon.under.obs)

# list.dual.exp <- list("Baseline" = exp.baseline,
#                       names(list.exp)[no.select.exp] = exp.extreme)

dual.ice <- lapply(list.dual.exp, FUN=function(name){
  
    # name <- list.dual.exp[[1]]
  dir.experiment          <- paste0(dir.output, name, "/")
  print(dir.experiment)
  models.fit      <- load.models(path=dir.experiment, split.type="FIT")
  std.const.fit   <- readRDS(file=paste0(dir.experiment, "standardization_constant_FIT.rds"))
  correct.names <- names(models)[which(names(models.fit) %in% models | names(models.fit) %in% names(models))]
  names(models.fit) <- correct.names
  
  # input.env.factors <- get.input.env.factors(exp.name=name)
  input.env.factors <- input.env.factors
  
  ice.dfs <- plot.ice(models.performance=models.fit,
                      select.env.fact=select.env.fact,
                      taxa=taxon.under.obs,
                      standardization.constant=std.const.fit[[1]],
                      observations = preprocessed.data.fit$`Entire dataset`,
                      nb.sample=no.sites,
                      resolution=no.steps,
                      input.env.factors=input.env.factors)
  
  return(ice.dfs)
})

obs1      <- dual.ice[[1]][["observations"]]
env.fact1 <- dual.ice[[1]][["env.factor.sampled"]]
obs2      <- dual.ice[[2]][["observations"]]
env.fact2 <- dual.ice[[2]][["env.factor.sampled"]]

merged.obs       <- bind_rows(list(obs1, obs2),
                              .id = "column_label")
merged.env.fact  <- bind_rows(list(env.fact1, env.fact2),
                              .id = "column_label")

observations.mean         <- merged.obs %>%
  group_by(across(all_of(select.env.fact)), model, column_label) %>%
  summarise(avg = mean(pred))
observations.mean.bounds  <- observations.mean %>% group_by(model, column_label) %>%
  summarise(x.mean=max(across(all_of(name.select.env.fact))),
            y.mean.min=min(avg),
            y.mean.max=max(avg))

list.obs.mean <- list(observations.mean %>% filter(column_label==1),
                      observations.mean %>% filter(column_label==2))

min.boundaries <- lapply(list.obs.mean, FUN=function(ice){
  min(ice[[name.select.env.fact]])
})

max.boundaries <- lapply(list.obs.mean, FUN=function(ice){
  max(ice[[name.select.env.fact]])
})

lb <- max(unlist(min.boundaries))
hb <- min(unlist(max.boundaries))

names.scenarios <- names(list.dual.exp)

merged.obs["column_label"] <- ifelse(merged.obs[["column_label"]]==1, names.scenarios[1], names.scenarios[2]) 
observations.mean["column_label"] <- ifelse(observations.mean[["column_label"]]==1, names.scenarios[1], names.scenarios[2])
merged.env.fact["column_label"] <- ifelse(merged.env.fact[["column_label"]]==1, names.scenarios[1], names.scenarios[2])
observations.mean.bounds["column_label"] <- ifelse(observations.mean.bounds[["column_label"]]==1, names.scenarios[1], names.scenarios[2])

plot.data <- merged.obs
plot.data$model <- factor(plot.data$model, levels= names(models))

fig4 <- ggplot(data=plot.data) +
  geom_line(aes(x=.data[[select.env.fact]],
                y=pred,
                group=observation_number, 
                color=as.character(observation_number)),
            alpha=0.4,
            show.legend = FALSE) +
  geom_line(data=observations.mean,
            aes(x=.data[[name.select.env.fact]], y=avg),
            size=1.5) +
  geom_rug(data = merged.env.fact,
           aes(x=variable),
           alpha=0.7,
           inherit.aes=F) + 
  geom_segment(data=observations.mean.bounds,
               inherit.aes = FALSE,
               lineend="round",
               linejoin="round",
               aes(x=hb,
                   y=y.mean.min,
                   xend=hb,
                   yend=y.mean.max),
               arrow=arrow(length = unit(0.3, "cm"),
                           ends = "both"),
               alpha=0.9) +
  xlim(lb, hb) +
  facet_wrap(~column_label+model) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +#,
  # legend.title = element_text(size=24),
  # legend.text = element_text(size=20)) +
  labs(x =name.select.env.fact,
       y = "Predicted probability of occurrence", 
       title = taxon.under.obs)
# fig4

pdf(paste0(dir.compar.plots, filename, "_ice.pdf"))
print(fig4)
dev.off()

rm(fig4)

# repeat fig4 for each model separately
for (model.name in c("GLM", "GAM", "RF")){ # c("glm", "gamloess", "rf", "ann")){
  
  filtered.merged.obs               <- merged.obs               %>% filter(model == model.name)
  filtered.observations.mean        <- observations.mean        %>% filter(model == model.name)
  filtered.observations.mean.bounds <- observations.mean.bounds %>% filter(model == model.name)
  
  fig <- ggplot(data=filtered.merged.obs) +
    geom_line(aes(x=.data[[select.env.fact]],
                  y=pred,
                  group=observation_number, 
                  color=as.character(observation_number)),
              show.legend = FALSE) +
    geom_line(data=filtered.observations.mean,
              aes(x=.data[[name.select.env.fact]], y=avg),
              size=1.5) +
    geom_rug(data = merged.env.fact,
             aes(x=variable), 
             color="grey20",
             alpha=0.7,
             inherit.aes=F) + 
    geom_segment(data=filtered.observations.mean.bounds,
                 inherit.aes = FALSE,
                 lineend="round",
                 linejoin="round",
                 aes(x=hb,
                     y=y.mean.min,
                     xend=hb,
                     yend=y.mean.max),
                 arrow=arrow(length = unit(0.3, "cm"),
                             ends = "both")) +
    xlim(lb, hb) +
    facet_wrap(~column_label) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +#,
    # legend.title = element_text(size=24),
    # legend.text = element_text(size=20)) +
    labs(x =name.select.env.fact,
         y = "Predicted probability of occurrence")
  
  # dir <- "../output_data/comparison_plots/"
  
  pdf(paste0(dir.compar.plots, filename, "_ice_", model.name, ".pdf"), width = 8, height = 4)
  print(fig)
  dev.off()
  
  rm(fig)
  
}

# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Fig 5: multi box plot ---
# multi.all.results <- lapply(list.exp, FUN=function(name){
#   
#   dir.experiment <- paste0(dir.output, name, "/")
#   models.cv      <- load.models(path=dir.experiment, split.type="CV")
#   
#   all.results    <- summarize.all.results(models.cv, data.prev.taxa)
#   
#   rm(models.cv) 
#   
#   all.results <- restructure.all.results(all.results, models)
#   
#   all.results["noise"] <- name
#   
#   return(all.results)
# })
# 
# final.multi.all.results <- bind_rows(multi.all.results, .id = "column_label")
# 
# fig5 <- ggplot(data=final.multi.all.results) +
#   geom_boxplot(aes(x=column_label,
#                    y=dev,
#                    fill=fit_pred)) +
#   scale_x_discrete(limits=names(list.exp)) +
#   facet_wrap(~model) + 
#   theme_minimal() +
#   theme(legend.title=element_blank())
# 
# 
# pdf(paste0(dir.compar.plots, filename, "_boxplot.pdf"))
# print(fig5)
# dev.off() 
# 
# rm(fig5)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fig 6: multi-scenario bellplots ----
multi.all.results <- lapply(list.exp, FUN=function(name){
  
  dir.experiment <- paste0(dir.output, name, "/")
  models.cv      <- load.models(path=dir.experiment, split.type="CV")
  
  all.results    <- summarize.all.results(models.cv, data.prev.taxa)
  
  rm(models.cv) 
  
  all.results <- restructure.all.results(all.results, models)
  
  all.results["noise"] <- name
  
  return(all.results)
})

final.multi.all.results <- bind_rows(multi.all.results, .id = "column_label")


fig6 <- ggplot(data=final.multi.all.results,
               aes(x=prevalence, y=dev, color=column_label)) + 
  geom_point(data=final.multi.all.results) +
  facet_wrap(~model + fit_pred)  +
  scale_color_manual(values=color.map) + 
  theme_minimal() + 
  theme(legend.title=element_blank())

pdf(paste0(dir.compar.plots, filename, "_bellplot.pdf"))
print(fig6)
dev.off() 

rm(fig6)

