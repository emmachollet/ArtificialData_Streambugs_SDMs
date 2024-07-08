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
name.streambugs.run     <- "10Catch_3081Sites_update.traits.hybrid_rk4Met_50Yea_18251Steps"
file.input.data         <- paste0(name.streambugs.run, "_WideData_ResultsThreshPresAbs.csv")
file.prev.taxa          <- paste0(name.streambugs.run, "_PrevalenceTaxonomy_ThreshPresAbs.csv")
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
env.factor <- c(vect.dir.env.fact, vect.indir.env.fact, "Temp2" = "tempmaxC2", "FV2" = "currentms2")
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
  
  # "noise.disp"     = list("type"       = "dispersal",
  #                        "target"     = NULL,
  #                        "amount"     = NULL,
  #                        "parameters" = NULL) # ,

  # "noise.temp"     = list("type"       = "gaussian",
  #                         "target"     = "tempmaxC",
  #                         "amount"     = 1,
  #                         "parameters" = list("min"=0, "max"=30)) #,
  # # 
  # "noise.gamm"     = list("type"       = "missdetection",
  #                         "target"     = "Occurrence.Gammaridae",
  #                         "amount"     = 0.1,
  #                         "parameters" = NULL)#,
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

name.list.noise <- paste(names(list.noise), collapse = "_")
  
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
# PREPROCESS DATA ####
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# preprocess input dataset
data.input <- data.env.taxa
data.input$tempmaxC2 <- data.input$tempmaxC^2
data.input$currentms2 <- data.input$currentms^2
data.input <- data.input[,c(vect.info,env.factor,taxa.colnames)] # subselect columns of interest
data.input$Watershed <- gsub("Doux", "Doubs", data.input$Watershed) # correct cute mistake
catch.variable <- "Watershed"

# add noise
data.input.noised        <- add.noise(data.input,
                                      number.sample,
                                      list.noise,
                                      env.factor,
                                      env.factor.full)
data.input <- data.input.noised$`noised data`

# summary(data.input$tempmaxC)
# reduce dataset to selected number of samples
rind.select <- runif(n=number.sample, min=1, max=dim(data.input)[1])
data.input <- data.input[rind.select,]

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

plot.data <- filter(plot.data, Taxon %in% "Gammaridae")
q <- ggplot(plot.data, aes(x = Watershed, fill = Observation))
q <- q + geom_bar(stat = "count", position = "stack")
q <- q + scale_fill_manual(values=vect.col.pres.abs)
q <- q + theme_bw()
q <- q + scale_x_discrete(limits = temp.catch)
q <- q + theme(axis.text.x = element_text(angle=90))
q
file.name <- "GeomBar_PresAbs_GammaridaeCatch"
print.pdf.plots(list.plots = list(q), width = 7, height = 5, dir.output = dir.experiment, info.file.name = "", file.name = file.name, 
                png = F)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Split and standardize data ####

# Example for Johanna
# ~~~~~~~~~~~~~~~~~~~

# # define number of splits from the beginning and use it in all functions and loops (to avoid "magic numbers")
# no.splits <- 5
# 
# # make list of 5 equal folds based on column "ReachID", with no interescting indices
# list.folds <- createFolds(data.input$ReachID, no.splits) 
# 
# # make empty list of splits and name them, that we fill later in a loop
# list.splits <- vector(mode = "list", length = no.splits)
# names(list.splits) <- paste0("Split", 1:no.splits)
# 
# # make a loop over splits to make and insert training and testing datasets
# for (n in 1:no.splits) {
#   # n <- 1 # useful to have a commented indices to go through the loop line by line once to be sure it does what you want
#   rind.test <- list.folds[[n]] # extract the row indices from the folds list
#   
#   temp.test <- data.input[rind.test,] # take the row indices of fold for a temporary test set
#   temp.train <- data.input[-rind.test,] # take all other row indices for a temporary train set
#   
#   list.splits[[n]][["Training data"]] <- temp.train # insert your training set as Training data in the nth split
#   list.splits[[n]][["Testing data"]] <- temp.test # insert your test set as Testing data in the nth split
# }

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

## Variable for plots ----
all.results <- summarize.all.results(models.cv, prev.taxa)

ggplot.all.results <- restructure.all.results(all.results)

color.map <- c('GLM'           = 'deepskyblue',   # Generalized Linear Model
               'GAM'      = 'green',         # Generalized Additive Model
               # 'ann'           = 'orange',        # Artificial Neural Network
               'RF'            = 'red')#,           # Random Forest
               # 'null'          = 'black')         # Null Model


## Fig. 1 : boxplot perf ----
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


## Fig. 2 : bellplot ----
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

## Fig. 3 : ICE ----


model.fit <- models.fit$GLM$entire_dataset$training$Baetisalpinus$model
pdp::partial(model.fit, pred.var = "tempmaxC", plot = TRUE, rug = TRUE)



list.plots <- list()

for(taxon.under.obs in names(taxa.colnames)){
  
  cat("Plotting ICE for:", taxon.under.obs)
  
  # select.env.fact <- "temperature"
  select.env.fact <- "tempmaxC"
  
  # taxon.under.obs <- names(taxa.colnames)[1]
  
  # input.env.factors <- get.input.env.factors(experiment.name)
  input.env.factors <- env.factor
  
  ice.dfs <- plot.ice(models.performance=models.fit,
                      select.env.fact=select.env.fact,
                      taxa=taxon.under.obs,
                      standardization.constant=std.const.fit[[1]],
                      observations = preprocessed.data.fit[["Entire dataset"]][,input.env.factors],
                      nb.sample=100,
                      resolution=200,
                      input.env.factors=input.env.factors)
  
  observations              <- ice.dfs[["observations"]]
  env.factor.sampled        <- ice.dfs[["env.factor.sampled"]]
  observations.mean         <- observations %>%
    group_by(across(all_of(select.env.fact)), model) %>%
    summarise(avg = mean(pred))
  observations.mean.bounds  <- observations.mean %>% group_by(model) %>%
    summarise(x.mean=max(across(all_of(select.env.fact))),
              y.mean.min=min(avg),
              y.mean.max=max(avg))
  
  plot.data <- observations
  # plot.data$model <- factor(plot.data$model, levels=c("GLM", "GAM", "RF"))
  
  fig3 <- ggplot(data=plot.data) +
    geom_line(aes(x=.data[[select.env.fact]],
                  y=pred,
                  group=observation_number, 
                  color=as.character(observation_number)),
              show.legend = FALSE) +
    geom_line(data=observations.mean,
              aes(x=.data[[select.env.fact]], y=avg),
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
    facet_wrap(~factor(model, levels = c("GLM", "GAM", "RF"))) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          legend.title = element_text(size=24),
          legend.text = element_text(size=20)) +
    labs(title = taxon.under.obs,
         x = "Temperature",
         y = "Predicted probability of occurrence")
  
  list.plots[[taxon.under.obs]] <- fig3
  
  pdf(paste0(dir.experiment, "ice_", experiment.name, taxon.under.obs, ".pdf"))
  print(fig3)
  dev.off()
}  

file.name <- paste0("ice_alltaxa")
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

# options for plots comparison

# create directory for comparison plots
dir.compar.plots <- paste0(dir.output, "comparison_plots/")
dir.create(dir.compar.plots)

model.color.map <- c('GLM'     = "#619CFF",  # 'deepskyblue',   # Generalized Linear Model
                     'GAM'     = "#00BA38", # 'green',         # Generalized Additive Model
                     # 'ANN'           = 'orange',        # Artificial Neural Network
                     'RF'      = "#F8766D") # 'red')#,           # Random Forest
                     # 'Null'          = 'black')         # Null Model

# select taxon and env. factor
taxon.under.obs <- "Baetisalpinus" # names(taxa.colnames)[7]
select.env.fact <- env.factor[1]
name.select.env.fact <- names(select.env.fact)
# name.select.env.fact <- names(env.factor)
# set.seed(97)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Recover Streambugs response shape ----

data.base.ICE <- read.csv(paste0(dir.input.data, "WideDataInputs_ICE_100Sites.csv"), sep = ";")
ice.df.streambugs <- read.csv(paste0(dir.input.data, 
                                     "10Catch_2000Sites_ICE_Prev0.15_10Yea_3651Steps_WideData_ResultsProbObs.csv"), sep = ";")

rind.na <- which(is.na(ice.df.streambugs$Occurrence.Gammaridae))
sites.na <- ice.df.streambugs[rind.na, c("ReachID", env.factor[1:8])]

list.plots <- list()

for(taxon.under.obs in names(taxa.colnames)){
  
  # taxon.under.obs <- names(taxa.colnames)[4]
  
  pred.ice.streamb <- ice.df.streambugs[,c("ReachID","Watershed", select.env.fact, paste0("Occurrence.", taxon.under.obs))] %>%
    mutate(observation_number = rep(1:100, each = 20),
           column_label = 1,
           model = "Streambugs") %>%
    rename(pred = paste0("Occurrence.", taxon.under.obs))
  
  env.factor.sampled <- data.frame(variable = data.base.ICE[,select.env.fact])
  
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
    geom_line(data=pred.ice.mean,
              aes(x=.data[[name.select.env.fact]], y=avg),
              size=1.5) +
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
  facet_wrap(~Watershed) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white"),
          legend.title = element_text(size=24),
          legend.text = element_text(size=20)) +
    labs(title = taxon.under.obs,
         x = "Temperature",
         y = "Predicted probability of occurrence")
  # fig3
  list.plots[[taxon.under.obs]] <- fig3
}

file.name <- paste0("ice_streambugs_alltaxa")
print.pdf.plots(list.plots = list.plots, width = 6, height = 6, 
                dir.output = dir.output, info.file.name = "", file.name = file.name, 
                png = F)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Comparison different noise scenarios ----

# compare experiment dispersal noise
list.exp     <- list("Baseline"                            = "3000Sites_15Taxa_3SDMs_10EnvFact_",
                     "Reduce dataset size"                 = "1000Sites_15Taxa_3SDMs_10EnvFact_",
                     "Missdetection"                       = "3000Sites_15Taxa_3SDMs_10EnvFact_noise.gamm",
                     "Absence from dispersal limitation"   = "3000Sites_15Taxa_3SDMs_10EnvFact_noise.disp")

# test.colors <- RColorBrewer::brewer.pal(8, "Set1")
# print(test.colors)
# scales::show_col(test.colors)

# color.map <- RColorBrewer::brewer.pal(length(list.exp), "Dark2")
# scales::show_col(color.map)

color.map <- c('1'            = 'deepskyblue',
               '2'            = 'green',
               '3'            = 'orange',
               '4'            = 'red')

color.map <- color.map[1:length(list.exp)]
names(color.map) <- names(list.exp)


# create.comparison.plots("gauss_temp", list.exp.gauss, 
#                         file.prev.taxa="8Catch_1416Sites_curve.curr.temp-10_interc.orgmic.sapro4_10Yea_3651Steps_PrevalenceTaxonomy_ThreshPresAbs.csv",
#                         taxon="Gammaridae",
#                         env.fact="tempmaxC")


## Fig 1 & Fig 2: PDP comparison ----


multi.ice <- lapply(list.exp, FUN=function(name){
  
  # name <- list.exp[[1]]
  dir.experiment          <- paste0(dir.output, name, "/")
  
  models.fit      <- load.models(path=dir.experiment, split.type="FIT")
  std.const.fit   <- readRDS(file=paste0(dir.experiment, "standardization_constant_FIT.rds"))
  
  input.env.factors <- env.factor
  
  ice.dfs <- plot.ice(models.performance=models.fit,
                      select.env.fact=select.env.fact,
                      taxa=taxon.under.obs,
                      standardization.constant=std.const.fit[[1]],
                      observations = preprocessed.data.fit$`Entire dataset`,
                      nb.sample=100,
                      resolution=200,
                      input.env.factors=input.env.factors)
  
  observations              <- ice.dfs[["observations"]]
  observations.mean         <- observations %>%
    group_by(across(all_of(select.env.fact)), model) %>%
    summarise(avg = mean(pred))
  
  observations.mean["noise"] <- name
  
  return(observations.mean)
})

final.multi.ice <- bind_rows(multi.ice, .id = "column_label")

# compute boundaries
min.boundaries <- lapply(multi.ice, FUN=function(ice){
  min(ice[[name.select.env.fact]])
})

max.boundaries <- lapply(multi.ice, FUN=function(ice){
  max(ice[[name.select.env.fact]])
})

lb <- max(unlist(min.boundaries))
hb <- min(unlist(max.boundaries))

plot.data <- final.multi.ice
plot.data$model <- factor(plot.data$model, levels = c("GLM", "GAM", "RF"))
plot.data$column_label <- factor(plot.data$column_label, levels = unlist(names(list.exp)))

fig1 <- ggplot(data=plot.data) +
  geom_line(aes(x=.data[[name.select.env.fact]],
                y=avg,
                group=model,
                colour=model),
            size=1) +
  facet_wrap(~factor(column_label, levels=unlist(names(list.exp))), ncol = 4)+
  #facet_wrap(~column_label) +
  xlim(lb, hb) +
  # scale_y_continuous(limits = c(0,1)) +
  scale_color_manual(values=model.color.map) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +#,
        # legend.title = element_text(size=24),
        # legend.text = element_text(size=20)) +
  labs(x = name.select.env.fact,
       y = "Predicted probability of occurrence",
       colour="Models")

filename <- paste0("comparison_", 
                   length(list.exp), "exp_",
                   taxon.under.obs, "_",
                   select.env.fact)
pdf(paste0(dir.compar.plots, filename, "_pdp_per_scenario.pdf"), width = 11, height = 4)
print(fig1)
dev.off()

rm(fig1)



fig2 <- ggplot(data=plot.data) +
  geom_line(aes(x=.data[[name.select.env.fact]],
                y=avg,
                group=column_label,
                colour=column_label),
            size=1) +
  scale_color_manual(values=color.map) +
  xlim(lb, hb) +
  # scale_y_continuous(limits = c(0,1)) +
  facet_wrap(~factor(model, levels = c("GLM", "GAM", "RF"))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +#,
  # legend.title = element_text(size=24),
  # legend.text = element_text(size=20)) +
  labs(x =name.select.env.fact,
       y = "Predicted probability of occurrence",
       colour="Scenario")

pdf(paste0(dir.compar.plots, filename, "_pdp_per_model.pdf"), width = 9, height = 4)
print(fig2)
dev.off()

rm(fig2)


# Fig 3: plotting the score of given taxon
multi.all.results <- lapply(list.exp, FUN=function(name){
  
  dir.experiment <- paste0(dir.output, name, "/")
  models.cv      <- load.models(path=dir.experiment, split.type="CV")
  
  all.results    <- summarize.all.results(models.cv, data.prev.taxa)
  
  rm(models.cv) 
  
  all.results <- restructure.all.results(all.results)
  
  all.results["noise"] <- name
  
  return(all.results)
})

final.multi.all.results <- bind_rows(multi.all.results, .id = "column_label")

filtered.multi.all.results <- final.multi.all.results %>%
  filter(taxa == taxon)

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
       y="Standardized deviance") +
  #facet_wrap(~fit_pred) + 
  theme_minimal() +
  theme(legend.title=element_blank())


pdf(paste0(dir.compar.plots, filename, "_taxa_score.pdf"))
print(fig3)
dev.off()

rm(fig3)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fig 4: dual ICE ----
exp.baseline <- list.exp[[1]]
exp.extreme  <- list.exp[[length(list.exp)]]


taxon.under.obs <- paste0("Occurrence.", taxon.under.obs)

list.dual.exp <- list("Baseline" = exp.baseline,
                      "Dispersal limitation"  = exp.extreme)

dual.ice <- lapply(list.dual.exp, FUN=function(name){
  
  dir.experiment          <- paste0(dir.output, name, "/")
  
  models.fit      <- load.models(path=dir.experiment, split.type="FIT")
  std.const.fit   <- readRDS(file=paste0(dir.experiment, "standardization_constant_FIT.rds"))
  
  # input.env.factors <- get.input.env.factors(exp.name=name)
  input.env.factors <- input.env.factors
  
  ice.dfs <- plot.ice(models.performance=models.fit,
                      select.env.fact=select.env.fact,
                      taxa=taxon.under.obs,
                      standardization.constant=std.const.fit[[1]],
                      observations = preprocessed.data.fit$`Entire dataset`,
                      nb.sample=100,
                      resolution=200,
                      input.env.factors=input.env.factors)
  
  return(ice.dfs)
})

obs1      <- dual.ice[["Baseline"]][["observations"]]
env.fact1 <- dual.ice[["Baseline"]][["env.factor.sampled"]]
obs2      <- dual.ice[["Dispersal limitation"]][["observations"]]
env.fact2 <- dual.ice[["Dispersal limitation"]][["env.factor.sampled"]]

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

merged.obs["column_label"] <- ifelse(merged.obs[["column_label"]]==1, "Baseline", "High noise") 
observations.mean["column_label"] <- ifelse(observations.mean[["column_label"]]==1, "Baseline", "High noise")
merged.env.fact["column_label"] <- ifelse(merged.env.fact[["column_label"]]==1, "Baseline", "High noise")
observations.mean.bounds["column_label"] <- ifelse(observations.mean.bounds[["column_label"]]==1, "Baseline", "High noise")


fig4 <- ggplot(data=merged.obs) +
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
       y = "Predicted probability of occurrence")

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

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fig 5: multi box plot ----
multi.all.results <- lapply(list.exp, FUN=function(name){
  
  dir.experiment <- paste0(dir.output, name, "/")
  models.cv      <- load.models(path=dir.experiment, split.type="CV")
  
  all.results    <- summarize.all.results(models.cv, data.prev.taxa)
  
  rm(models.cv) 
  
  all.results <- restructure.all.results(all.results)
  
  all.results["noise"] <- name
  
  return(all.results)
})

final.multi.all.results <- bind_rows(multi.all.results, .id = "column_label")

fig5 <- ggplot(data=final.multi.all.results) +
  geom_boxplot(aes(x=column_label,
                   y=dev,
                   fill=fit_pred)) +
  scale_x_discrete(limits=names(list.exp)) +
  facet_wrap(~model) + 
  theme_minimal() +
  theme(legend.title=element_blank())


pdf(paste0(dir.compar.plots, filename, "_boxplot.pdf"))
print(fig5)
dev.off() 

rm(fig5)


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Fig 6: multi-scenario bellplots ----
multi.all.results <- lapply(list.exp, FUN=function(name){
  
  dir.experiment <- paste0(dir.output, name, "/")
  models.cv      <- load.models(path=dir.experiment, split.type="CV")
  
  all.results    <- summarize.all.results(models.cv, data.prev.taxa)
  
  rm(models.cv) 
  
  all.results <- restructure.all.results(all.results)
  
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

