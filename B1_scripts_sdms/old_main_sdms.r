# ###############################################################################|
# #                                                                              #
# #  --- Study of the influence of noise on overfitting of machine learning ---  #
# #            --- and statistical species distribution models. ---              #
# #                                                                              #
# #                          --- Octobre 2024 ---                               #
# #                                                                              #
# #                   --- Emma Chollet, Gaspard Fragnière ---                    #
# #                --- Andreas Scheidegger and Nele Schuwirth ---                #
# #                                                                              #
# #                      --- emma.chollet@eawag.ch ---                           #
# #                                                                              #
# ###############################################################################|
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # PRELIMINARIES ####
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# Sys.setenv(LANG="EN")
# set.seed(13)   # Always set seed to a lucky number
# getwd()        # Show working directory. It needs to be the location of 'main.r'
# rm(list=ls())  # Free work space
# graphics.off() # Clean graphics display
# 
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#  ## Libraries ####
# 
# if (!require("dplyr")){install.packages("dplyr"); library("dplyr")}                    # to sort, join, merge data
# if (!require("tidyr") ) { install.packages("tidyr"); library("tidyr") }               # to sort, join, merge data
# if (!require("ggplot2")){install.packages("ggplot2"); library("ggplot2")}                    # to do nice plots
# if (!require("ggpubr")){install.packages("ggpubr"); library("ggpubr")}                    # to do nice plots
# if (!require("gridExtra")){install.packages("gridExtra"); library("gridExtra")}                    # to do nice plots
# if (!require("readr")){install.packages("readr"); library("readr")}
# if (!require("jsonlite")){install.packages("jsonlite"); library("jsonlite")}
# if (!require("pROC")){install.packages("pROC"); library("pROC")}  # to compute AUC    
# if (!require("vtable") ) { install.packages("vtable"); library("vtable") }               # to do tables of summary statistics
# if (!require("sf") ) { install.packages("sf"); library("sf") }                        # to read layers for plotting maps
# if (!require("ggpattern")){install.packages("ggpattern"); library("ggpattern")}        # to add patterns in boxplots
# 
# # specific for Neural Network
# # if (!require("reticulate")){install.packages("reticulate"); library("reticulate")}    # links R to python
# library("reticulate")
# # reticulate::install_python(version = "3.10")  # or your desired version
# # install_miniconda()              # run this the very first time reticulate is installed # makes a separeta python installation
# # remotes::install_github("rstudio/keras")
# # keras3::install_keras(envname = "C:/Users/cholleem/.virtualenvs/r-keras", backend = "tensorflow")
# use_virtualenv("C:/Users/cholleem/.virtualenvs/r-keras", required = TRUE)
# # keras3::install_keras(backend = "tensorflow")
# 
# library("keras3")
# 
# # specific for Generalize Additive Model
# if (!require("mgcv")){install.packages("mgcv"); library("mgcv")}
# if (!require("gam")){install.packages("gam"); library("gam")}
# if (!require("splines")){install.packages("splines"); library("splines")}
# if (!require("foreach")){install.packages("foreach"); library("foreach")}
# 
# # caret has to be loaded at the end to not cache function 'train'
# if (!require("caret")){install.packages("caret"); library("caret")}
# 
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Directories and files ####
# 
# # define directories
# dir.input.data          <- "../A3_outputs_streambugs/"
# dir.output              <- "../B2_outputs_sdms/"
# dir.utilities           <- "../00_utilities/"
# # dir.create(dir.output)
# 
# # define files
# name.streambugs.run     <- "8Catch_3009Sites_PolyInterp30_runC_100Yea_36501Steps"
# file.input.data         <- paste0(name.streambugs.run, "_WideData_ResultsSamplePresAbs.csv")
# file.prev.taxa          <- paste0(name.streambugs.run, "_PrevalenceTaxonomy_SamplePresAbs.csv")
# file.selected.taxa      <- "selected_taxa_analysis.csv"
# 
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Load data and functions ####
# 
# # read data
# data.env.taxa       <- read.csv(paste0(dir.input.data, file.input.data), header = T, sep = ";", stringsAsFactors = F)
# data.prev.taxa      <- read.csv(paste0(dir.input.data, file.prev.taxa), header = T, sep = ";", stringsAsFactors = F)
# data.selected.taxa  <- read.csv(paste0(dir.utilities, file.selected.taxa), header = T, sep = ";", stringsAsFactors = F)
# 
# # load functions
# source(paste0(dir.utilities, "utilities_global.r"))
# source("data_noising.r")
# source("performances_assessment.r")
# source("training_pipeline.r")
# source("ml_models.r")
# source("global_variables.r")
# source("plotting.r")
# 
# # prepare inputs for plotting results on swiss maps
# map.inputs      <- map.inputs(directory = paste0(dir.utilities,"swiss.map.gdb"))
# 
# # set options for plots
# # width and height of an A4 page in inches
# width.a4 = 8.3
# height.a4 = 11.7
# plot.suffix <- c(".png", ".pdf")
# 
# mytheme <- theme_bw() +
#     theme(text = element_text(size = 14),
#           title = element_text(size = 14),
#           # plot.title = element_text(face = "bold"),
#           # strip.text = element_text(size = 14), #, face = "bold"),
#           strip.background = element_rect(fill = "white"),
#           strip.text.x = element_text(size = 14),
#           # axis.title = element_text(size = 14), #, face = "bold"),
#           # axis.text.x = element_text(vjust = 10),
#           # axis.text.y = element_text(hjust = 10),
#           # axis.title.x.bottom = element_text(vjust = -15),
#           # axis.title.x = element_text(margin = margin(t = 0, r = 0, b = -30, l = 0)), # put some space between the axis title and the numbers
#           axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
#           legend.position="bottom"
#           # legend.text = element_text(size = 14),
#           # legend.title = element_text(size = 14) #, face= "bold") #,
#           # panel.grid.major = element_line(colour = NA),
#           # panel.grid.minor = element_line(colour = NA)
#     )
# 
# theme_set(mytheme)
# 
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Noise scenarios ####
# list.scenarios <- list(
#     
#     "dataset.size"         = list("flag"  = T,
#                                   "range" = c("best" = 3000, 
#                                               "good" = 2100, 
#                                               "mid"  = 1200, 
#                                               "bad"  = 300)),
#     
#     "nb.predictors"        = list("flag"  = F,
#                                   "range" = c("best" = 8, 
#                                               "good" = 6, 
#                                               "mid"  = 4, 
#                                               "bad"  = 2)),
#     
#     "noise.temperature"    = list("flag"  = F,
#                                   "range" = c("best" = 0, 
#                                               "good" = 1.5, 
#                                               "mid"  = 3, 
#                                               "bad"  = 4.5)),
#     
#     "misdetection"         = list("flag"  = F,
#                                   "range" = c("best" = 0, 
#                                               "good" = 15, 
#                                               "mid"  = 30, 
#                                               "bad"  = 45))
#     
# )
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Inputs, training options ####
# 
# # environmental factors used for streambugs
# vect.dir.env.fact <- c("Temperature"              = "tempmaxC",
#                        "Flow velocity"            = "currentms",
#                        "Insecticide pollution"    = "orgmicropollTU",
#                        "Water quality"            = "saprowqclass")
# vect.indir.env.fact <- c("Shade"          = "shade",
#                          "Litter input"   = "Lit_Inp" ,
#                          "Phosphorus"     = "C_P",
#                          "Nitrogen"       = "C_N")
# vect.info <- colnames(data.env.taxa)[1:5]
# print(vect.info)
# 
# # environmental factors selected for sdms
# env.factor <- c(vect.dir.env.fact, vect.indir.env.fact)
# env.factor.full <- c(env.factor, "Temp2" = "tempmaxC2", "FV2" = "currentms2")
# no.env.fact <- length(env.factor)
# 
# # taxa list full
# cind.taxa <- which(grepl("Occurrence.", colnames(data.env.taxa)))
# vect.taxa.full <- colnames(data.env.taxa)[cind.taxa]
# names(vect.taxa.full) <- gsub("Occurrence.", "", vect.taxa.full)
# no.taxa.full <- length(vect.taxa.full)
# print(vect.taxa.full)
# 
# # update number of NAs in input data
# for (taxon in vect.taxa.full) {
#     # taxon <- vect.taxa.full[11]
#     data.prev.taxa[which(data.prev.taxa$Occurrence.taxa == taxon), "Missing.values"] <- sum(is.na(data.env.taxa[,taxon]))
# }
# summary(data.prev.taxa$Missing.values)
# 
# # taxa list selected for analysis and above prevalence threshold
# prev.threshold <- 0.1
# na.threshold <- 800
# temp.occ.taxa <- data.prev.taxa[which(data.prev.taxa$Prevalence.NoDisp > prev.threshold
#                                       & data.prev.taxa$Prevalence.NoDisp < 1 -prev.threshold
#                                       & data.prev.taxa$Missing.values < na.threshold), "Occurrence.taxa"]
# 
# # selected.taxa.analysis <- data.selected.taxa$Occurrence.taxa
# # taxa.colnames <- selected.taxa.analysis[which(selected.taxa.analysis %in% temp.occ.taxa)] # to select taxa listed manually and with a prev. and NA threshold
# taxa.colnames <- temp.occ.taxa
# names(taxa.colnames) <- gsub("Occurrence.", "", taxa.colnames)
# names.taxa <- names(taxa.colnames)
# no.taxa <- length(taxa.colnames)
# print(no.taxa)
# # print(taxa.colnames)
# 
# # select taxon and env. factor for later analysis
# taxon.under.obs <- names(taxa.colnames)[5]
# select.env.fact <- env.factor[1]
# name.select.env.fact <- names(select.env.fact)
# cat("Taxon selected for analysis:", taxon.under.obs)
# cat("Environmental predictor selected for analysis:", select.env.fact)
# 
# # select colors for present/absent plots
# vect.col.pres.abs <- c( "Present" = "#00BFC4", "Absent" = "#F8766D")
# scales::show_col(vect.col.pres.abs)
# 
# # sdms training options
# number.split            <- 3
# split.criterion         <- "ReachID"
# number.sample           <- 3000
# number.sample           <- ifelse(dim(data.env.taxa)[1] < number.sample, dim(data.env.taxa)[1], number.sample)
# sdm.models              <- c(
#                                 "GLM" = "glm",
#                                 # "GAM" = 'gamSpline',
#                                 "GAM" = "gamSpline",
#                                 "RF"  = "rf",
#                                 # "downRF" = "downrf")#,
#                                 "ANN" = "ann")#,
# # "ann")
# no.models <- length(sdm.models)
# # no.models  <- 4 # test
# models <- append(c("Null" = "null"), sdm.models)
# 
# # tune ann hyperparameters
# tune.ann <- F
# 
# # noise scenarios
# list.noise <- list(
#   
#   # "noise.disp"     = list("type"       = "dispersal",
#   #                        "target"     = NULL,
#   #                        "amount"     = NULL,
#   #                        "parameters" = NULL) # ,
# 
#   # "noise.temp"     = list("type"       = "gaussian",
#   #                         "target"     = "tempmaxC",
#   #                         "amount"     = 3,
#   #                         "parameters" = list("min"=0, "max"=30)) #,
#   # # 
#   # "noise.gamm"     = list("type"       = "missdetection",
#   #                         "target"     = "Occurrence.Gammaridae",
#   #                         "amount"     = 0.1,
#   #                         "parameters" = NULL)#,
#   
#   # # 
#   # "noise.add.fact" = list("type"       = "add_factor",
#   #                         "target"     = "random1", # new factor name
#   #                         "amount"     = NULL,
#   #                         "parameters" = NULL),
#   # 
#   # "noise.rem.fact" = list("type"       = "remove_factor",
#   #                         "target"     = "currentms",
#   #                         "amount"     = NULL,
#   #                         "parameters" = NULL)#,
#   
# )
# name.list.noise <- paste(names(list.noise), collapse = "_")
# 
# # p <- 0.1 # probability at which a presence is turned into an absence
# # 
# # list.noise <- lapply(taxa.colnames, FUN=function(taxon){
# #     noise_taxon <- list("type"       = "missdetection",
# #                         "target"     = taxon,
# #                         "amount"     = p,
# #                         "parameters" = NULL)
# # 
# #     return(noise_taxon)
# # })
# # name.list.noise <- paste0("misdetection.all.taxa", p)
#   
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # PREPROCESS DATA ####
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Add noise ####
# 
# # preprocess input dataset
# data.input <- data.env.taxa
# data.input$tempmaxC2 <- data.input$tempmaxC^2
# data.input$currentms2 <- data.input$currentms^2
# data.input <- data.input[,c(vect.info,env.factor.full,taxa.colnames)] # subselect columns of interest
# catch.variable <- "Watershed"
# 
# # add noise
# data.input.noised        <- add.noise(data = data.input,
#                                       number.sample,
#                                       noise = list.noise,
#                                       env.fact = env.factor,
#                                       env.fact.full = env.factor.full)
# data.input      <- data.input.noised$`noised data`
# env.factor      <- data.input.noised$`env fact`
# env.factor.full <- data.input.noised$`env fact full`
# no.env.fact <- length(env.factor)
# 
# # reduce dataset to selected number of samples
# rind.select <- sample(dim(data.input)[1], size=number.sample, replace = F)
# data.input <- data.input[rind.select,]
# 
# # replace NAs by 0 to simulate noise due to dispersal limitation
# if("noise.disp" %in% names(list.noise)){
#   data.input[is.na(data.input)] <- 0 # replace NAs by 0
# }
# 
# # replace 0/1 by absent/present
# for (taxon in taxa.colnames) {
#   data.input[which(data.input[,taxon] == 0),taxon] <- "Absent"
#   data.input[which(data.input[,taxon] == 1),taxon] <- "Present"
#   data.input[,taxon] = as.factor(data.input[,taxon])
# }
# 
# # correct coordinates from lv03 to lv95
# rind <- which(data.input$MonitoringProgram == "SyntheticPoint")
# data.input[rind, "X"] <- data.input[rind, "X"] - 2000000
# data.input[rind, "Y"] <- data.input[rind, "Y"] - 1000000
# 
# # sanity check of NAs
# sum(is.na(data.input[,c(vect.info,env.factor)])) # should be zero 
# sum(is.na(data.input[,c(taxa.colnames)])) # could have a lot, decide to replace by zero or not
# nb.rows.na.taxa <- sum(rowSums(is.na(data.input[,c(taxa.colnames)])) != 0)
# 
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Write metadata ####
# 
# # define experiment name
# experiment.name <- paste0(
#     number.sample, "Sites_",
#     no.taxa, "Taxa_",
#     no.models, "SDMs_",
#     no.env.fact, "EnvFact_",
#     name.list.noise,
#     ifelse(tune.ann, "TuneANN_", "")
# )
# experiment.name <- paste0(experiment.name,"testDownsampling_1500trees")
# print(experiment.name)
# 
# # create output directory for experiment results
# dir.experiment          <- paste0(dir.output, experiment.name, "/")
# dir.metadata            <- paste0(dir.experiment, "metadata.json")
# dir.create(dir.experiment)
# 
# # saving metadata to JSON file
# metadata.list           <- list("input_data"       = file.input.data,
#                                 "prev_taxa"        = file.prev.taxa,
#                                 "nb_rows_na"       = nb.rows.na.taxa,
#                                 "experiment.name"  = experiment.name,
#                                 "number_split"     = number.split,
#                                 "split.criterion"  = split.criterion,
#                                 "number.sample"    = number.sample,
#                                 "models"           = models,
#                                 "env_factor"       = env.factor,
#                                 "env_factor_full"  = env.factor.full,
#                                 "noise"            = list.noise) 
# 
# metadata.json           <- toJSON(metadata.list)
# write(metadata.json, dir.metadata)
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Prevalence analysis ####
# 
# # colnames(data.input)
# # plot.data <- data.input
# # taxon <- taxa.colnames[5]
# # plot.data[,taxon] <- as.factor(plot.data[,taxon])
# # p <- ggplot(data.input, aes(x = .data[[env.factor[1]]], y = .data[[env.factor[2]]], color = .data[[taxon]])) +
# #     geom_point()
# # p
# 
# # 
# # plot.suffix <- c(".png", ".pdf")
# # 
# # # plot of sites analysis (catchments and programs)
# # plot.data <- data.input %>%
# #     mutate(Site = ifelse(data.input$MonitoringProgram == "SyntheticPoint", "Additional", "Biomonitoring")) %>%
# #     mutate(Site = factor(Site, levels = c("Biomonitoring", "Additional")))
# # summary(plot.data$Site)
# # 
# # 
# # p <- ggplot() + 
# #     geom_sf(data = map.inputs$ch, fill="#E8E8E8", color="black") + 
# #     geom_sf(data = map.inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE) + 
# #     geom_sf(data = map.inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE) + 
# #     geom_point(data = plot.data, aes(x=X, y=Y, color = Watershed, shape = Site),  size= 1, alpha= 0.8) +
# #     scale_shape_manual(values = c(16,4)) +
# #     theme(legend.position = "right") +
# #     labs(color = "Catchment")
# # p
# # 
# # # save png and pdf
# # for(suffix in plot.suffix){
# #     file.name <- paste0(dir.experiment, "Plot_Swiss_map_Catch_Progr", suffix)
# #     ggsave(file.name, width = 2480*1.1,
# #            height = 3508*0.5,
# #            units = c("px"))
# # }
# # 
# # # plot of distribution of environmental predictors
# # plot.data <- data.input[,c(env.factor)]
# # name.env.fact <- names(env.factor)
# # colnames(plot.data) <- name.env.fact
# # plot.data <- gather(plot.data, key = factor, value = value) %>%
# #     mutate(factor = factor(factor, levels = name.env.fact))
# # 
# # q <- ggplot(data = plot.data, aes(x=value)) +
# #     geom_histogram(colour = "red", fill = alpha("red", 0.3))+#, 
# #                       # alpha = 0.3) +
# #     facet_wrap(~factor, scales = "free", ncol = 2) +
# #                     # labeller = env.fact_labeller, 
# #                     # strip.position="top") +
# #     labs(x = "Values of environmental input",
# #          y = "Count"#,
# #               # title = "Distribution of the environmental factors in splits"
# # )
# # q
# # 
# # # save png and pdf
# # for(suffix in plot.suffix){
# #     file.name <- paste0(dir.experiment, "Plot_Distribution_EnvPred", suffix)
# #     ggsave(file.name, width = 2480*1,
# #            height = 3508*0.7,
# #            units = c("px"))
# # }
# # 
# # # plot bars pres/abs/NA
# # # make long dataframe for analysis
# # long.df.input.taxa <- gather(data.input, key = Taxon, value = Observation, all_of(taxa.colnames))
# # long.df.input.taxa$Taxon <- gsub("Occurrence.", "", long.df.input.taxa$Taxon)
# # long.df.input.env.taxa <- gather(long.df.input.taxa, key = Factor, value = Value, all_of(env.factor))
# # 
# # plot.data <- long.df.input.taxa
# # plot.data$Observation <- factor(plot.data$Observation, levels=c("Present", NA, "Absent"))
# # 
# # temp.taxa <- unique(plot.data$Taxon)
# # temp.catch <- unique(plot.data$Watershed)
# # 
# # # analyze distribution of prevalence and NAs across catchments
# # p <- ggplot(plot.data, aes(x = Taxon, fill = Observation)) +
# #     geom_bar(stat="count", position='stack') +
# #     scale_x_discrete(limits=rev) +
# #     # facet_wrap(~ Watershed, scales = "free_y") +
# #     scale_fill_manual(values=vect.col.pres.abs) + 
# #     coord_flip(ylim= c(0, number.sample)) +
# #     labs(y = "Counts",
# #          fill = "")
# # p
# # 
# # # save png and pdf
# # for(suffix in plot.suffix){
# #     file.name <- paste0(dir.experiment, "Plot_Barplots_PresAbsNA_Taxa", suffix)
# #     ggsave(file.name, width = 2480*1,
# #            height = 3508*0.7,
# #            units = c("px"))
# # }
# # 
# # # analyze distribution of prevalence and NAs across taxa
# # q <- ggplot(plot.data, aes(x = Watershed, fill = Observation))
# # q <- q + geom_bar(stat = "count", position = "stack")
# # q <- q + facet_wrap(~ Taxon, scales = "free_y")
# # q <- q + scale_fill_manual(values=vect.col.pres.abs)
# # q <- q + theme_bw()
# # q <- q + scale_x_discrete(limits = temp.catch)
# # q <- q + theme(axis.text.x = element_text(angle=90))
# # q
# # 
# # # print plots in pdf
# # file.name <- "GeomBar_PresAbs_TaxaCatch"
# # print.pdf.plots(list.plots = list(q,p), width = 12, height = 9, dir.output = dir.experiment, info.file.name = "", file.name = file.name,
# #                 png = F)
# # 
# # plot.data <- filter(plot.data, Taxon %in% "Gammaridae")
# # q <- ggplot(plot.data, aes(x = Watershed, fill = Observation))
# # q <- q + geom_bar(stat = "count", position = "stack")
# # q <- q + scale_fill_manual(values=vect.col.pres.abs)
# # q <- q + theme_bw()
# # q <- q + scale_x_discrete(limits = temp.catch)
# # q <- q + theme(axis.text.x = element_text(angle=90))
# # q
# # file.name <- "GeomBar_PresAbs_GammaridaeCatch"
# # print.pdf.plots(list.plots = list(q), width = 7, height = 5, dir.output = dir.experiment, info.file.name = "", file.name = file.name,
# #                 png = F)
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Split and standardize data ####
# 
# preprocessed.data.cv  <- preprocess.data(data=data.input,
#                                          env.fact=env.factor,
#                                          dir=dir.experiment,
#                                          split.type="CV",
#                                          nb.split=number.split,
#                                          splitting.criterion=split.criterion)
# 
# preprocessed.data.fit <- preprocess.data(data=data.input,
#                                          env.fact=env.factor,
#                                          dir=dir.experiment,
#                                          split.type="FIT")
# 
# name.datasets <- c("FIT", names(preprocessed.data.cv))
# vect.add.col <- apply(expand.grid(c("prev", "na"), name.datasets), 1, paste, collapse=".")
# data.prev.na <- data.prev.taxa[which(data.prev.taxa$Occurrence.taxa %in% taxa.colnames),]
# data.prev.na[, vect.add.col] <- NA
# for (taxon in taxa.colnames) {
#     # taxon <- taxa.colnames[6]
#     for (name in name.datasets) {
#         # name <- name.datasets[1]
#         if(name == "FIT"){
#             temp.df <- preprocessed.data.fit[[1]]
#         } else {
#             temp.df <- preprocessed.data.cv[[name]][["Training data"]]
#         }
#         nb.na <- sum(is.na(temp.df[,taxon]))
#         # test <- temp.df[,taxon] == "Present"
#         nb.pres <- sum(temp.df[,taxon] == "Present", na.rm = T)
#         nb.abs <- sum(temp.df[,taxon] == "Absent", na.rm = T)
#         temp.prev <- round(nb.pres/(nb.pres + nb.abs), digits = 3)
#         
#         data.prev.na[which(data.prev.na$Occurrence.taxa == taxon), paste0("prev.", name)] <- temp.prev
#         data.prev.na[which(data.prev.na$Occurrence.taxa == taxon), paste0("na.", name)] <- nb.na
#     }
# }
# 
# # write summarized results in csv
# file.name <- paste0(experiment.name, "prevalence_NAs.csv")
# write.table(data.prev.na, file = paste0(dir.experiment, file.name), sep = ",", row.names = F)
# 
# ## Recover ICE data ####
# # recover Streambugs data for ICE
# data.base.ice            <- read.csv(paste0(dir.input.data, "WideDataInputs_ICE_50Sites.csv"), sep = ";")
# data.base.ice$observation_number <- 1:nrow(data.base.ice)
# data.base.ice$tempmaxC2  <- data.base.ice$tempmaxC^2
# data.base.ice$currentms2 <- data.base.ice$currentms^2
# 
# # preprocess input data for ICE plots
# preprocessed.data.ice <- preprocess.data(data=data.base.ice,
#                                              env.fact=env.factor,
#                                              dir=dir.experiment,
#                                              split.type="ICE")
# 
# name.ice.streambugs <- "8Catch_50Sites_ICE_50Steps_PolyInterp30_runC_100Yea_36501Steps_"
# ice.df.streambugs   <- read.csv(paste0(dir.input.data, name.ice.streambugs,
#                                        "WideData_ResultsProbObs.csv"), sep = ";")
# no.sites <- as.numeric(stringr::str_match(name.ice.streambugs, "Catch_(.+)Sites")[2])
# no.steps <- as.numeric(stringr::str_match(name.ice.streambugs, "ICE_(.+)Steps_Poly")[2])
# 
# pred.ice.streamb <- ice.df.streambugs[,c("ReachID", "Watershed", select.env.fact, taxa.colnames)] %>% #, paste0("Occurrence.", taxon))] %>%
#     mutate(#observation_number = rep(1:no.sites, each = no.steps),
#            column_label = 1,
#            model = "Streambugs") # %>%
# pred.ice.streamb$observation_number <- NA
# for(reach in data.base.ice$ReachID){
#     # reach <- data.base.ice$ReachID[1]
#     obs.num <- data.base.ice[which(data.base.ice$ReachID == reach), "observation_number"]
#     rind <- which(grepl(reach, pred.ice.streamb$ReachID))
#     pred.ice.streamb[rind, "observation_number"] <- obs.num
# }
# env.factor.sampled <- data.frame(variable = data.base.ice[,select.env.fact])
# 
# # # compute pdp for the selected taxon
# # mean.ice.streamb <- pred.ice.streamb %>%
# #     group_by_at(select.env.fact[1]) %>%
# #     summarize(across(all_of(taxa.colnames), mean, na.rm = TRUE)) %>%
# #     mutate(model = "Streambugs")
# # mean.ice.streamb.bounds  <- mean.ice.streamb %>%
# #     # group_by(model) %>%
# #     summarize(x.mean=max(Temperature), # ! here it's specific to factor
# #         y.mean.min=across(all_of(names.taxa), min, na.rm = TRUE),
# #         y.mean.max=across(all_of(names.taxa), max, na.rm = TRUE)) %>%
# #     mutate(# model = "Streambugs",
# #            x.mean = 21.94111) # ! here it's specific to factor to fit the other models
# 
# # points.ice.streamb <- data.base.ice[,c("ReachID", select.env.fact[1])]
# # points.ice.streamb$pred <- NA
# # for(site in points.ice.streamb$ReachID){
# #     
# # }
# 
# # list.plots.stream <- list()
# # names.taxa <- names(taxa.colnames)
# # 
# # for(taxon in names.taxa){
# # 
# #     # taxon <- names.taxa[1]
# #     # extract ice dataframe for taxon selected for analysis
# #     pred.ice.streamb <- ice.df.streambugs[,c("ReachID", "Watershed", select.env.fact, paste0("Occurrence.", taxon))] %>%
# #         mutate(observation_number = rep(1:no.sites, each = no.steps),
# #                column_label = 1,
# #                model = "Streambugs") %>%
# #         rename(pred = paste0("Occurrence.", taxon))
# # 
# #     env.factor.sampled <- data.frame(variable = data.base.ice[,select.env.fact])
# # 
# #     # compute pdp for the selected taxon
# #     mean.ice.streamb <- pred.ice.streamb %>%
# #         group_by_at(select.env.fact[1]) %>%
# #         summarize(avg = mean(na.omit(pred))) %>%
# #         mutate(model = "Streambugs")
# #     mean.ice.streamb.bounds  <- mean.ice.streamb %>%
# #         # group_by(model) %>%
# #         summarize(x.mean=max(Temperature), # ! here it's specific to factor
# #             y.mean.min=min(avg),
# #             y.mean.max=max(avg)) %>%
# #         mutate(# model = "Streambugs",
# #                x.mean = 21.94111) # ! here it's specific to factor to fit the other models
# # 
# #     # recover preference for temperature
# #     name.pref.temp <- "df.preferences_PolyInterp30__266Taxa_tempmaxtolval.csv"
# #     df.pref.temp <- read.csv(paste0(dir.input.data, name.pref.temp), sep = ",")
# #     df.pref.temp.taxon <- df.pref.temp %>%
# #         select(Values, all_of(taxon)) %>%
# #         filter(Values >= min(mean.ice.streamb[, name.select.env.fact]) & Values <= max(mean.ice.streamb[, name.select.env.fact])) %>%
# #         rename(PrefTrait = all_of(taxon)) %>%
# #         mutate(ExpTrans = exp.transform(PrefTrait, intercept = 0, curv = -10))  # Applying the log function to column2
# # 
# #     plot.data <- pred.ice.streamb
# #     plot.data.mean  <- mean.ice.streamb
# #     plot.data.mean.bounds  <- mean.ice.streamb.bounds
# #     plot.data.pref <- df.pref.temp.taxon
# #     # plot.data$model <- factor(plot.data$model, levels=c("GLM", "GAM", "RF"))
# # 
# #     fig3 <- ggplot(data=plot.data) +
# #         geom_line(aes(x=.data[[select.env.fact]],
# #                       y=pred,
# #                       group=observation_number,
# #                       color=as.character(observation_number)),
# #                   show.legend = FALSE, alpha = 0.7) +
# #         geom_line(data=plot.data.mean,
# #                   aes(x=.data[[name.select.env.fact]], y=avg),
# #                   size=1.5) +
# #         geom_rug(data = env.factor.sampled,
# #                  aes(x=variable),
# #                  color="grey20",
# #                  alpha=0.7,
# #                  inherit.aes=F) +
# #         geom_point(data = plot.data.pref,
# #                    aes(x = Values, y = ExpTrans),
# #                    shape = 4, color = "black", size = 2, stroke = 0.8) +
# #         scale_x_continuous(limits = c(min(env.factor.sampled), max(env.factor.sampled))) +
# #         geom_segment(data=plot.data.mean.bounds,
# #                      inherit.aes = FALSE,
# #                      lineend="round",
# #                      linejoin="round",
# #                      aes(x=x.mean,
# #                          y=y.mean.min,
# #                          xend=x.mean,
# #                          yend=y.mean.max),
# #                      arrow=arrow(length = unit(0.3, "cm"),
# #                                  ends = "both")) +
# #         facet_wrap(~factor(model, levels = c("Streambugs", names(sdm.models))), ncol = 5) +
# #         # theme_bw() +
# #         # theme(strip.background = element_rect(fill = "white"),
# #         #       legend.title = element_text(size=24),
# #         #       legend.text = element_text(size=20)) +
# #         labs(title = taxon,
# #              x = "Temperature",
# #              y = "Predicted probability of occurrence")
# # 
# #     list.plots.stream[[taxon]] <- fig3
# #     # fig3
# # }
# # 
# # # Figure SI A 5
# # # file.name <- "ice_streambugs_all_taxa_grid"
# # # pdf(paste0(dir.output, file.name, ".pdf"), width = width.a4*1.2, height = height.a4*1.2)
# # # grid.arrange(grobs=list.plots.stream, ncol = 5)
# # # dev.off()
# # 
# # p <- grid.arrange(grobs =list.plots.stream, ncol = 5)
# # 
# # # save png and pdf
# # for(suffix in plot.suffix){
# #     file.name <- paste0(dir.output, "ice_streambugs_all_taxa_grid", suffix)
# #     ggsave(file.name, plot = p, width = 2480*2.2,
# #            height = 3508*2.2,
# #            units = c("px"))
# # }
# 
# # list.plots.stream[[1]]
# # 
# # file.name <- "ice_streambugs"
# # print.pdf.plots(list.plots = list.plots.stream, width = width.a4*0.7, height = height.a4*0.5,
# #                 dir.output = dir.output, info.file.name = "",
# #                 file.name = file.name)
# 
# # pdf(paste0(dir.output, "ice_streambugs_", taxon.under.obs, ".pdf"), width = 5, height = 5)
# # print(fig3)
# # dev.off()
# # ggsave(paste0(dir.output, "ice_streambugs_", taxon.under.obs, ".png"), width = 1200,
# #        height = 1100,
# #        units = c("px"))
# 
# # # check balance presence/absence
# # for(i.split in 1:number.split){
# #     cat("Split number:", i.split, "\n")
# #     for (taxon in taxa.colnames) {
# #         cat(taxon)
# #         print(summary(preprocessed.data.cv[[i.split]][["Training data"]][,taxon]))
# #     }
# # }
# # print(summary(preprocessed.data.cv[[1]][["Training data"]][,"Occurrence.Baetisalpinus"]))
# # Absent Present 
# # 1149     836
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # APPLY MODELS ####
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## RF downsampling ####
# 
# # ## number of data points and simulations
# # ntr = 300 # number of training presence-absence records
# # ntst = 1000 # number of testing presence-absence records
# # # nsims = 50 # number of simulation/species
# # 
# # allsamp <- 1:dim(data.input)[1]
# # trainsamp <- sample(allsamp, ntr) # make sure train and test are independent
# # allsamp <- allsamp[-trainsamp]
# # testsamp <- sample(allsamp, ntst)
# # data_train <- data_full[trainsamp, ]
# # data_test <- data_full[testsamp, ]
# # data_train$occ <- as.factor(data_train$occ)
# # # fit the RF default
# # rf_def <- randomForest(formula = occ ~ ., 
# #                        data = data_train, 
# #                        ntree = 500)
# # rfpred_def <- predict(rf_def, data_test[,2:3], type = "prob") #predict
# # myAUC_def<-auc(data.frame(1:ntst, data_test[, 1], rfpred_def[, 2]))
# # 
# # # fit RF down-sample
# # prNum <- as.numeric(table(data_train$occ)["1"]) # nr presences
# # bgNum <- as.numeric(table(data_train$occ)["0"]) # nr absences
# # smpsize <- c("0" = min(prNum, bgNum), "1" = min(prNum, bgNum))  ##GGA: modified so works if num1>num0
# # rf_downsample <- randomForest(formula = occ ~ ., 
# #                               data = data_train, 
# #                               ntree = 1000, 
# #                               sampsize = smpsize,
# #                               replace = TRUE)
# # rfpred_downsample <- predict(rf_downsample, data_test[,2:3], type = "prob") #predict
# # myAUC_downsample <- auc(data.frame(1:ntst,data_test[,1],rfpred_downsample[,2]))
# 
# # training <- preprocessed.data.cv$Split1$`Training data`
# # testing <- preprocessed.data.cv$Split1$`Testing data`
# # taxon <- taxa.colnames[1]
# # table(training[,taxon])
# # nmin <- min(table(training[,taxon]))
# # # rows.pres <- which(training[,taxon] == "Present")
# # # rows.abs <- which(training[,taxon] == "Absent")
# # # samp.pres <- sample(rows.pres, nmin)
# # 
# # ctrl <- trainControl(method='cv',                  
# #                           number= 3, #nb.folds,                     
# #                           # index=folds.train,
# #                           classProbs=T,                 
# #                           summaryFunction=standard.deviance,
# #                           selectionFunction=lowest,
# #                           # verboseIter=TRUE)
# #                           verboseIter=F)
# # formul <- paste0(taxon,
# #                  " ~ ",
# #                  paste(env.factor, collapse = ' + '))
# # 
# # rfDownsampled <- train(formula(formul), data = training,
# #                        method = "rf",
# #                        ntree = 1500,
# #                        tuneLength = 5,
# #                        # metric = "ROC",
# #                        trControl = ctrl,
# #                        ## Tell randomForest to sample by strata. Here, 
# #                        ## that means within each class
# #                        strata = training[,taxon],
# #                        ## Now specify that the number of samples selected
# #                        ## within each class should be the same
# #                        sampsize = rep(nmin, 2))
# # 
# # rfUnbalanced <- train(formula(formul), data = training,
# #                       method = "rf",
# #                       ntree = 1500,
# #                       tuneLength = 5,
# #                       # metric = "ROC",
# #                       trControl = ctrl)
# # 
# # downPred <- predict(rfDownsampled, training)
# # downProbs <- predict(rfDownsampled, training, type = 'prob')
# # 
# # # Save model's performances
# # downPerf <- performance.wrapper(model=rfDownsampled,
# #                                           observation=training,
# #                                           prediction=downPred,
# #                                           prediction.probability=downProbs,
# #                                           taxa=taxon)
# # 
# # unbalPred <- predict(rfUnbalanced, training)
# # unbalProbs <- predict(rfUnbalanced, training, type = 'prob')
# # 
# # # Save model's performances
# # unbalPerf <- performance.wrapper(model=rfUnbalanced,
# #                                           observation=training,
# #                                           prediction=unbalPred,
# #                                           prediction.probability=unbalProbs,
# #                                           taxa=taxon)
# # downPerf$standard_deviance
# # unbalPerf$standard_deviance
# # 
# # downPerf$auc
# # unbalPerf$auc
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Tune ANN ####
# if(tune.ann){
#     tune.grid.ann          <- expand.grid(num.units                 = c(8),
#                                           num.layers                = c(1,2,3,5),
#                                           learning.rate             = c(0.01),
#                                           num.epochs                = c(100, 200),
#                                           batch.size                = c(64))
#     
#     list.trained.ann <- list()
#     for (i.row in seq(nrow(tune.grid.ann))) {
#         
#         # i.row <- 1
#         vect.hyperparam <- tune.grid.ann[i.row,]
#         name.model <- paste0(vect.hyperparam["num.units"], "uni_",
#                              vect.hyperparam["num.layers"], "lay_",
#                              vect.hyperparam["num.epochs"], "epo")
#         print(name.model)
#         
#         trained.ann <- apply.ann.model(data = preprocessed.data.cv, split.type = "CV", 
#                                        taxa.colnames = taxa.colnames, env.fact = env.factor, 
#                                        hyperparameters = vect.hyperparam)
#         
#         list.trained.ann[[name.model]] <- trained.ann
#     }
#     
#     # update models list
#     models <- names(list.trained.ann)
#     names(models) <- models
#     models.cv <- list.trained.ann
# }
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Cross-validation ####
# 
# tune.grid.ann          <- expand.grid(num.units                 = c(8),
#                                       num.layers                = c(2),
#                                       learning.rate             = c(0.01),
#                                       num.epochs                = c(100),
#                                       batch.size                = c(64))
# hyperparameters <- tune.grid.ann[1,]
# 
# models.cv  <- apply.ml.models(data=preprocessed.data.cv,
#                               models=models,
#                               split.type="CV", 
#                               taxa.colnames = taxa.colnames, 
#                               env.factor = env.factor, 
#                               env.factor.full = env.factor.full,
#                               hyperparameters = hyperparameters)
# 
# save.models(models=models.cv,
#             path=dir.experiment,
#             split.type="CV")
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Fit entire dataset ####
# 
# models.fit <- apply.ml.models(data=preprocessed.data.fit,
#                               models=models,
#                               split.type="FIT",
#                               taxa.colnames, 
#                               env.factor, env.factor.full,
#                               hyperparameters = hyperparameters)
# 
# save.models(models=models.fit,
#             path=dir.experiment,
#             split.type="FIT")
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # PLOTS EXPERIMENT ####
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# # models.fit      <- load.models(path=dir.experiment,
# #                                split.type="FIT")
# # models.cv       <- load.models(path=dir.experiment,
# #                                split.type="CV")
# 
# std.const.cv    <- readRDS(file=paste0(dir.experiment, "standardization_constant_CV.rds"))
# std.const.fit   <- readRDS(file=paste0(dir.experiment, "standardization_constant_FIT.rds"))
# prev.taxa <- data.prev.taxa
# prev.taxa$Prevalence.DispLim <- prev.taxa$Prevalence
# prev.taxa$Prevalence <- prev.taxa$Prevalence.NoDisp
# 
# ## Variable for plots ----
# all.results <- summarize.all.results(models.cv, prev.taxa)
# 
# # write summarized results in csv
# file.name <- paste0(experiment.name, "allresults.csv")
# write.table(all.results, file = paste0(dir.experiment, file.name), sep = ";", row.names = F)
# 
# ggplot.all.results <- restructure.all.results(all.results, models)
# ggplot.all.results$model <- factor(ggplot.all.results$model, levels= names(models))
# 
# # set colors for models
# if(tune.ann){
#     set.seed(13)
#     color.map <- palette(rainbow(length(models))) 
#     scales::show_col(color.map)
# } else {
#     # color.map <- c('GLM'           = 'deepskyblue',   # Generalized Linear Model
#     #                'GAM'           = 'green',         # Generalized Additive Model
#     #                'ANN'           = 'orange',        # Artificial Neural Network
#     #                'RF'            = 'red',           # Random Forest
#     #                'Null'          = 'black')         # Null Model
#     color.map <- c('downRF'           = 'deepskyblue',   # Generalized Linear Model
#                    # 'GAM'           = 'green',         # Generalized Additive Model
#                    # 'ANN'           = 'orange',        # Artificial Neural Network
#                    'RF'            = 'red',           # Random Forest
#                    'Null'          = 'black')         # Null Model
# }
# 
# ## Fig. 1 : boxplot perf ----
# gg.plot.dev.infos <-  ggplot.all.results %>%
#   group_by(fit_pred, model) %>%
#   summarise(mean=mean(dev),
#             max=max(dev),
#             min=min(dev),
#             median=median(dev))
# if(tune.ann){
#     order.ann <- gg.plot.dev.infos %>%
#         filter(fit_pred == "pred") %>%
#         arrange(median) %>%
#         pull(model)
#     ggplot.all.results$model <- factor(ggplot.all.results$model, levels = order.ann)
# }
# 
# median.null.pred <- min((gg.plot.dev.infos %>%
#                            filter(fit_pred=="pred" & model=="Null"))["median"])
# median.best.pred <- min((gg.plot.dev.infos %>%
#                            filter(fit_pred=="pred"))["median"])
# 
# fig1 <- ggplot(data=ggplot.all.results,
#                aes(x=model, y=dev)) +
#     geom_boxplot(aes(fill=fit_pred)) +
#     geom_hline(aes(yintercept=median.null.pred,
#                  color="median null model (pred)"),
#              linetype="dashed") +
#     geom_hline(aes(yintercept=median.best.pred,
#                  color="median best model (pred)"),
#              linetype="dashed") +
#     scale_fill_manual(values=c(fit = "#998ec3", pred = "#f1a340")) +
#     scale_colour_manual(values = c('green','black')) +
#     theme_minimal() +
#     theme(legend.title=element_blank()) +
#     labs(x = "Models",
#          y = "Standardized deviance\n ",
#          title = "") 
# if(tune.ann){
#     fig1 <- fig1 + theme(axis.text.x = element_text(angle = 90))
# }
# 
# pdf(paste0(dir.experiment, "boxplot_performance.pdf"), width = 8, height = 5)
# print(fig1)
# dev.off()
# 
# 
# 
# ## Fig. 2 : bellplot ----
# fig2 <- ggplot(data=subset(ggplot.all.results, model %in% "Null"),
#                aes(x=prevalence, y=dev, color=model)) +
#   geom_smooth(se = FALSE) + 
#   geom_point(data=ggplot.all.results, alpha = 0.6) +
#   facet_wrap(~fit_pred) +
#   scale_color_manual(values=color.map) +
#   theme_minimal() + 
#   theme(legend.title=element_blank()) +
#     labs(x = "Prevalence",
#          y = "Standardized deviance\n ",
#          title = "")
# 
# 
# pdf(paste0(dir.experiment, "bellplot_performance.pdf"), width = 8, height = 5)
# print(fig2)
# dev.off() 
# 
# ## Fig. 3 : ICE ----
# 
# 
# dir.exp.ice <- paste0(dir.experiment, "ICE/")
# dir.create(dir.exp.ice)
# preprocessed.data.ice <- preprocess.data(data=data.base.ice,
#                                          env.fact=env.factor,
#                                          dir=dir.exp.ice,
#                                          split.type="FIT")
# std.const.ice   <- readRDS(file=paste0(dir.exp.ice, "standardization_constant_FIT.rds"))
# 
# 
# # select.env.fact <- "temperature"
# select.env.fact <- "tempmaxC"
# 
# list.plots <- list()
# 
# for(taxon.under.obs in names(taxa.colnames)){
#   
#   cat("Plotting ICE for:", taxon.under.obs)
# 
#   # taxon.under.obs <- names(taxa.colnames)[1]
#   
#   temp.pred.ice.streamb <- pred.ice.streamb[,c("ReachID", "Watershed", "column_label", "observation_number", "model",
#                                                select.env.fact, paste0("Occurrence.", taxon.under.obs))] %>%
#     # mutate(observation_number = rep(1:no.sites, each = no.steps),
#     #        column_label = 1,
#     #        model = "Streambugs") %>%
#     rename(pred = paste0("Occurrence.", taxon.under.obs))
#   # 
#   # env.factor.sampled <- data.frame(variable = data.base.ice[,select.env.fact])
#   
#   mean.ice.streamb <- temp.pred.ice.streamb %>%
#     group_by_at(select.env.fact[1]) %>%
#     summarize(avg = mean(na.omit(pred))) %>%
#     mutate(model = "Streambugs")
#   mean.ice.streamb.bounds  <- mean.ice.streamb %>%
#     # group_by(model) %>%
#     summarize(# x.mean=max(tempmaxC), # ! here it's specific to factor
#               y.mean.min=min(avg),
#               y.mean.max=max(avg)) %>%
#     mutate(model = "Streambugs",
#            x.mean = 21.94111) # ! here it's specific to factor to fit the other models
# 
#   # input.env.factors <- get.input.env.factors(experiment.name)
#   input.env.factors <- env.factor
#   
#   ice.dfs <- plot.ice(models.performance = models.fit,
#                       select.env.fact=select.env.fact,
#                       taxa=taxon.under.obs,
#                       standardization.constant=std.const.ice[[1]],
#                       observations = preprocessed.data.ice[["Entire dataset"]],
#                       nb.sample=no.sites,
#                       resolution=no.steps,
#                       input.env.factors=input.env.factors)
#   
#   pred.models              <- ice.dfs[["observations"]] %>%
#       rename(pred = all_of(taxon.under.obs))
#   for(reach in data.base.ice$ReachID){
#       # reach <- data.base.ice$ReachID[1]
#       obs.num <- data.base.ice[which(data.base.ice$ReachID == reach), "observation_number"]
#       rind <- which(grepl(reach, pred.models$ReachID))
#       pred.models[rind, "observation_number"] <- obs.num
#   }
#   env.factor.sampled       <- ice.dfs[["env.factor.sampled"]]
#   pred.models.mean         <- pred.models %>%
#     group_by(across(all_of(select.env.fact)), model) %>%
#     summarise(avg = mean(pred))
#   pred.models.mean.bounds  <- pred.models.mean %>% group_by(model) %>%
#     summarise(x.mean=max(across(all_of(select.env.fact))),
#               y.mean.min=min(avg),
#               y.mean.max=max(avg))
#   # colnames(pred.ice.streamb)
#   # colnames(pred.models)
#   col.names <- c(select.env.fact, "pred", "observation_number", "column_label", "model")
#   plot.data <- rbind(temp.pred.ice.streamb[,col.names], pred.models[,col.names] )
#   plot.data.mean  <- rbind(mean.ice.streamb, pred.models.mean)       
#   plot.data.mean.bounds  <- rbind(mean.ice.streamb.bounds, pred.models.mean.bounds)
#   # plot.data$model <- factor(plot.data$model, levels=c("GLM", "GAM", "RF"))
#   
#   fig3 <- ggplot(data=plot.data) +
#     geom_line(aes(x=.data[[select.env.fact]],
#                   y=pred,
#                   group=observation_number, 
#                   color=as.character(observation_number)),
#               show.legend = FALSE) +
#     geom_line(data=plot.data.mean,
#               aes(x=.data[[select.env.fact]], y=avg),
#               size=1.5) +
#     geom_rug(data = env.factor.sampled,
#              aes(x=variable), 
#              color="grey20",
#              alpha=0.7,
#              inherit.aes=F) + 
#     scale_x_continuous(limits = c(min(env.factor.sampled), max(env.factor.sampled))) +
#     geom_segment(data=plot.data.mean.bounds,
#                  inherit.aes = FALSE,
#                  lineend="round",
#                  linejoin="round",
#                  aes(x=x.mean,
#                      y=y.mean.min,
#                      xend=x.mean,
#                      yend=y.mean.max),
#                  arrow=arrow(length = unit(0.3, "cm"),
#                              ends = "both")) +
#     facet_wrap(~factor(model, levels = c("Streambugs", names(sdm.models))), ncol = 5) +
#     theme_bw() +
#     theme(strip.background = element_rect(fill = "white"),
#           legend.title = element_text(size=24),
#           legend.text = element_text(size=20)) +
#     labs(title = taxon.under.obs,
#          x = "Temperature",
#          y = "Predicted probability of occurrence")
#   # fig3
#   
#   list.plots[[taxon.under.obs]] <- fig3
#   
#   # pdf(paste0(dir.experiment, "ice_", experiment.name, taxon.under.obs, ".pdf"))
#   # print(fig3)
#   # dev.off()
# }  
# 
# file.name <- paste0("ice_alltaxa_3x1")
# print.pdf.plots(list.plots = list.plots, width = 12, height = 6, 
#                 dir.output = dir.experiment, info.file.name = experiment.name, file.name = file.name, 
#                 png = F)
# 
# # dev.off()
# # # Fig 4
# # noised.dataset <- models.fit[[1]][[1]][[1]][[1]][["observation"]]
# # taxa.col <- colnames(noised.dataset)[which(grepl("Occurrence.", colnames(noised.dataset)))]
# # 
# # noised.prevalence <- lapply(taxa.col, FUN=function(taxon){
# #   occurence  <- noised.dataset[[taxon]]
# #   prevalence <- length(occurence[occurence=="present"])/length(occurence) 
# #   return(prevalence)
# # })  
# # 
# # noised.prev.taxa <- prev.taxa[which(prev.taxa$Occurrence.taxa %in% taxa.colnames),]  
# # noised.prev.taxa[["Prevalence"]] <- noised.prevalence
# # 
# # prev.taxa[["Status"]] <- "unnoised"
# # noised.prev.taxa[["Status"]] <- "noised"
# # 
# # full.prev.taxa <- rbind(prev.taxa, noised.prev.taxa)
# # full.prev.taxa[["Occurrence.taxa"]] <- as.factor(full.prev.taxa[["Occurrence.taxa"]])
# # full.prev.taxa$Prevalence <- unlist(full.prev.taxa$Prevalence)
# # full.prev.taxa$Occurrence.taxa <- gsub("Occurrence.", "", full.prev.taxa$Occurrence.taxa)
# # 
# # fig4 <- ggplot(data=full.prev.taxa,
# #                aes(x=Prevalence,
# #                    y=Occurrence.taxa,
# #                    color=Status)) + 
# #   geom_point() +
# #   scale_y_discrete(limits = rev(gsub("Occurrence.", "", prev.taxa$Occurrence.taxa))) +
# #   labs(x="Prevalence", y="Taxa")
# # theme_minimal() + 
# #   theme(legend.title=element_blank())
# # 
# # 
# # pdf(paste0(dir.experiment, "change_prevalence.pdf"))
# # print(fig4)
# # dev.off()
# 
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # PLOTS COMPARISON ####
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Options plots ----
# 
# # create directory for comparison plots
# dir.compar.plots <- paste0(dir.output, "comparison_plots/")
# dir.create(dir.compar.plots)
# 
# # select taxon and env. factor for later analysis
# taxon.under.obs <- names(taxa.colnames)[5]
# print(taxon.under.obs)
# 
# # set color map for models
# model.color.map <- c('GLM'     = "#619CFF",  # 'deepskyblue',   # Generalized Linear Model
#                      'GAM'     = "#00BA38", # 'green',         # Generalized Additive Model
#                      'ANN'     = 'orange',        # Artificial Neural Network
#                      'RF'      = "#F8766D") # 'red')#,           # Random Forest
#                      # 'Null'          = 'black')         # Null Model
# model.color.map <- c('downRF'           = 'deepskyblue',   # Generalized Linear Model
#                # 'GAM'           = 'green',         # Generalized Additive Model
#                # 'ANN'           = 'orange',        # Artificial Neural Network
#                'RF'            = 'red')           # Random Forest
# scales::show_col(model.color.map)
# 
# # select experiments to compare
# list.exp     <- list("1. Best-case scenario"                     = "3000Sites_35Taxa_4SDMs_8EnvFact_",
#                      "2. Reduced dataset size"                    = "300Sites_35Taxa_4SDMs_8EnvFact_",
#                      "3. Removed environmental predictors"        = "3000Sites_35Taxa_4SDMs_4EnvFact_",
#                      # "Remove flow velocity"      = "3000Sites_45Taxa_4SDMs_7EnvFact_",
#                      "4. Noise on temperature"                   = "3000Sites_35Taxa_4SDMs_8EnvFact_noise.temp",
#                      "5. Misdetection"                           = "3000Sites_35Taxa_4SDMs_8EnvFact_misdetection.all.taxa0.1",
#                      "6. Absences due to dispersal limitation"   = "3000Sites_35Taxa_4SDMs_8EnvFact_noise.disp")
# names.exp <- names(list.exp)
# 
# # test.colors <- RColorBrewer::brewer.pal(8, "Set1")
# # scales::show_col(test.colors)
# 
# # color.map <- RColorBrewer::brewer.pal(length(list.exp), "Dark2")
# # scales::show_col(color.map)
# 
# color.map <- c('1'            = '#4daf4a',
#                '2'            = '#377eb8',
#                '3'            = '#984ea3',
#                '4'            = '#e41a1c',
#                '5'            = '#ff7f00',
#                '6'            = '#ffd92f')
# # 
# # color.map <- c('1'            = '#4daf4a',
# #                '2'            = '#377eb8',
# #                '3'            = '#984ea3',
# #                '4'            = '#d95f02')
# #                 # '5'            = '#ff7f00',
# #                 # '6'            = '#ffd92f')
# 
# color.map <- color.map[1:length(list.exp)]
# names(color.map) <- names(list.exp)
# scales::show_col(color.map)
# 
# # create file names for saving plots
# file.name.exp <- paste0("comparison_", 
#                         length(list.exp), "exp_")
# print(file.name.exp)
# file.name.tax <- paste0(file.name.exp,
#                         taxon.under.obs, "_",
#                         select.env.fact)
# print(file.name.tax)
# 
# # create.comparison.plots("gauss_temp", list.exp.gauss, 
# #                         file.prev.taxa="8Catch_1416Sites_curve.curr.temp-10_interc.orgmic.sapro4_10Yea_3651Steps_PrevalenceTaxonomy_ThreshPresAbs.csv",
# #                         taxon="Gammaridae",
# #                         env.fact="tempmaxC")
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ## Load results ----
# 
# # load and process all results
# multi.all.results <- lapply(list.exp, FUN=function(name){
#     
#     # name <- list.exp[[2]]
#     dir.experiment <- paste0(dir.output, name, "/")
#     file.name <- paste0(dir.experiment, name, "allresults.csv")
#     
#     all.results <- read.table(file.name, header = T, sep = ";")
#     for(i in 1:ncol(all.results)){
#         # i <- 1
#         correct.name <- gsub(models[1], names(models)[1],
#                              gsub(models[2], names(models)[2],
#                                   gsub(models[3], names(models)[3],
#                                        gsub(models[4], names(models)[4],
#                                             colnames(all.results)[i]))))
#         colnames(all.results)[i] <- correct.name
#     }
#     
#     restr.all.results <- restructure.all.results(all.results, models)
#     restr.all.results["noise"] <- name
#     
#     return(list(all.results, restr.all.results))
# })
# 
# select.results <- c("dev_pred", "auc_pred", "likelihood_ratio")
# summary.metric <- "median"
# table.summary.metrics <- bind_rows(lapply(names(list.exp), FUN = function(name, multi.all.results, select.results, summary.metric){
#     
#     # name <- names(list.exp)[1]
#     temp.res <- multi.all.results[[name]][[1]]
#     temp.names <- apply(expand.grid(names(models), select.results), 1, paste, collapse="_")
#     temp.summary <- temp.res[, temp.names] %>%
#         summarise_if(where(is.numeric), summary.metric, na.rm = TRUE) %>%
#         mutate(scenario = name) %>% 
#         mutate(across(is.numeric, round, digits=2))
#     return(temp.summary)
#     
# }, multi.all.results, select.results, summary.metric), .id = "id")
# 
# list.results.taxa <- lapply(names(taxa.colnames), FUN = function(taxon){
#     bind_rows(lapply(names(list.exp), FUN = function(name, multi.all.results, select.results){
#         # name <- names(list.exp)[1]
#         temp.res <- multi.all.results[[name]][[1]]
#         temp.names <- apply(expand.grid(names(models), select.results), 1, paste, collapse="_")
#         temp.summary <- temp.res[which(temp.res$taxa == taxon), c("taxa", "prevalence", temp.names)] %>%
#             mutate(scenario = name) %>% 
#             mutate(across(is.numeric, round, digits=2)) %>%
#             select(last_col(), everything())
#         return(temp.summary)
#         
#     }, multi.all.results, select.results), .id = "id")
# })
# 
# df.results.all.taxa <- bind_rows(list.results.taxa) %>%
#     arrange(prevalence, taxa)
# 
# # write results per taxon in csv
# file.name <- paste0(file.name.exp, "per_taxon", "_results.csv")
# write.table(df.results.all.taxa, file = paste0(dir.compar.plots, file.name), sep = ",", row.names = F)
# 
# # extract dataframe with best median for each scenario
# df.best.median <- table.summary.metrics[,c(paste0(names(models), "_dev_pred"), "scenario")]
# colnames(df.best.median) <- gsub("_dev_pred", "", colnames(df.best.median))    
# df.best.median <- df.best.median  %>% 
#     pivot_longer(
#         cols = names(models), 
#         names_to = "model",
#         values_to = "dev") %>%
#     rename(column_label = scenario) %>%
#     group_by(column_label) %>% 
#     slice(which.min(dev))
# df.best.median$column_label <- factor(df.best.median$column_label, levels = names.exp)
# df.best.median <- as.data.frame(df.best.median[order(df.best.median$column_label),])
# 
# # write summarized results in csv
# file.name <- paste0(file.name.exp, summary.metric, "_results.csv")
# write.table(table.summary.metrics, file = paste0(dir.compar.plots, file.name), sep = ",", row.names = F)
# 
# final.multi.all.results <- bind_rows(lapply(multi.all.results, "[[", 2), .id = "column_label")
# final.multi.all.results$model <- factor(final.multi.all.results$model, levels= names(models))
# final.multi.all.results$column_label <- factor(final.multi.all.results$column_label, levels = names.exp)
# 
# # # filter results for selected taxon
# # filtered.multi.all.results <- final.multi.all.results %>%
# #     filter(taxa == taxon.under.obs)
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Fig 1: multi box plot ----
# 
# plot.data <- final.multi.all.results
# plot.data$pattern <- ifelse(plot.data$fit_pred == "fit", "Calibration", "Prediction")
# plot.data$pattern <- factor(plot.data$pattern, levels= c("Calibration", "Prediction"))
# 
# fig5 <- ggplot(data=plot.data, aes(x=model, y=dev)) +
#     geom_boxplot_pattern(aes(fill = model, pattern = pattern),
#                          pattern_colour = 'black',
#                          pattern_fill = 'black',
#                          pattern_density = 0.1,
#                          # pattern_spacing = 0.01,
#                          pattern_key_scale_factor = 1.5
#                          # pattern = 'stripe', 
#                          # pattern_size = 0.2
#     ) +
#     geom_hline(data = df.best.median, aes(yintercept = dev, color = model), 
#                linetype="longdash", size = 1, alpha = 0.8, show.legend = FALSE) +
#     scale_x_discrete(limits=names(models)) +
#     scale_y_continuous(limits = c(0, 1.7)) +
#     scale_color_manual(values = model.color.map) +
#     facet_wrap(~ column_label, ncol = 3, , labeller = label_wrap_gen()) +
#     scale_pattern_manual(values= c("Calibration" = "stripe", "Prediction" = "none")) + #,
#     scale_fill_manual(values = c("Null" = "grey", model.color.map)) +
#     labs(x = "Model",
#          y = "Standardized deviance",
#          # color = "Best median",
#          fill = "",
#          pattern = "")
# # fig5
# 
# pdf(paste0(dir.compar.plots, file.name.exp, "boxplot_colormodels_stripesfit.pdf"), width = width.a4*1, height = height.a4*0.6)
# print(fig5)
# dev.off() 
# 
# # rm(fig5)
# 
# # ICE and PDP data ----
# 
# dir.compar.ice <- paste0(dir.compar.plots, "ice/")
# dir.create(dir.compar.ice)
# 
# # taxa.for.ice <- names.taxa[c(5,9,11)]
# taxa.for.ice <- names.taxa
# select.env.fact <- env.factor[1]
# name.select.env.fact <-  names(select.env.fact)
# 
# multi.ice <- list()
# 
# for (i in 1:length(list.exp)) {
# 
#     # i <- 1
#   exp <- list.exp[[i]]
#   name.exp <- names(list.exp)[i]
#   dir.experiment          <- paste0(dir.output, exp, "/")
#   cat("\nLoading models and computing ICE for", length(taxa.for.ice),"taxa and predictor", select.env.fact , "for experiment:", exp, "\n")
#   
#   models.fit      <- load.models(path=dir.experiment, split.type="FIT")
#   std.const.fit   <- readRDS(file=paste0(dir.experiment, "standardization_constant_FIT.rds"))
#   std.const.ice   <- readRDS(file=paste0(dir.experiment, "standardization_constant_ICE.rds"))
#   
#   correct.names <- names(models)[which(names(models.fit) %in% models | names(models.fit) %in% names(models))]
#   names(models.fit) <- correct.names
#   file.metadata <- "metadata.json"
#   metadata.exp <- fromJSON(paste(readLines(paste0(dir.experiment, file.metadata)), collapse=""))
#   
#   input.env.factors <- metadata.exp$env_factor
#   
#   ice.dfs <- plot.ice(models.performance=models.fit,
#                       select.env.fact=select.env.fact,
#                       # taxa=taxon.under.obs,
#                       taxa=taxa.for.ice,
#                       
#                       standardization.constant=std.const.ice[[1]],
#                       observations = preprocessed.data.ice$`Entire dataset`,
#                       nb.sample=no.sites,
#                       resolution=no.steps,
#                       input.env.factors=input.env.factors)
#   
#   observations              <- ice.dfs[["observations"]]
#   
#   observations.mean         <- observations %>%
#       group_by(across(all_of(select.env.fact)), model) %>%
#       summarize(across(all_of(taxa.for.ice), mean, na.rm = TRUE))
#     # group_by(across(all_of(select.env.fact)), model) %>%
#     # summarise(avg = mean(all_of(taxa.for.ice)))
#   
#   observations.mean["noise"] <- name.exp
#   observations["noise"]      <- name.exp
#   
#   multi.ice[[name.exp]] <- list(observations, observations.mean)
#   # return(list(observations, observations.mean))
#   
# }
# 
# list.ice.bound <- lapply(lapply(multi.ice, "[[", 2), function(obs.mean){
#     observations.mean.bounds  <- as.data.frame(obs.mean %>% group_by(model, noise) %>%
#         summarise(x.mean=max(across(all_of(name.select.env.fact))),
#                   y.mean.min=across(all_of(taxa.for.ice), min, na.rm = TRUE),
#                   y.mean.max=across(all_of(taxa.for.ice), max, na.rm = TRUE)))
# })
# 
# long.multi.ice <- bind_rows(lapply(multi.ice, "[[", 1), .id = "column_label_noise")
# long.multi.ice$model <- factor(long.multi.ice$model, levels = names(sdm.models))
# long.multi.ice$column_label_noise <- factor(long.multi.ice$column_label_noise, levels = names.exp)
# max.temp <- max(long.multi.ice[,select.env.fact])
# min.temp <- min(long.multi.ice[,select.env.fact])
# 
# multi.mean <- lapply(multi.ice, "[[", 2)
# final.multi.ice <- bind_rows(multi.mean, .id = "column_label")
# final.multi.ice$model <- factor(final.multi.ice$model, levels = names(sdm.models))
# final.multi.ice$column_label_noise <- factor(final.multi.ice$column_label, levels = names.exp)
# 
# final.multi.bound <- bind_rows(list.ice.bound, .id = "column_label") 
# final.multi.bound$x.mean <- round(max.temp, digits = 1)
# final.multi.bound$model <- factor(final.multi.bound$model, levels = names(sdm.models))
# final.multi.bound$column_label_noise <- factor(final.multi.bound$column_label, levels = unlist(names(list.exp)))
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Fig 2: ICE comparison ----
# 
# plot.data <- long.multi.ice
# plot.data.mean <- final.multi.ice
# plot.data.bounds <- final.multi.bound
# 
# # # sub-select models for plotting
# # plot.data <- long.multi.ice %>%
# #     filter(model %in% c("GAM", "RF"))
# # plot.data.mean <- final.multi.ice %>%
# #     filter(model %in% c("GAM", "RF"))
# # plot.data.bounds <- final.multi.bound %>%
# #     filter(model %in% c("GAM", "RF"))
# 
# plot.data.stream <- mean.ice.streamb %>%
#     filter(Temperature < max.temp,
#            Temperature > min.temp) # ! predictor specific
# env.factor.sampled <- data.frame(variable = unique(data.base.ice[,select.env.fact]))
# 
#     
# # plot.data
# # colnames(plot.data)
# # colnames(plot.data.mean)
# # colnames(plot.data.bounds)
# 
# # unique(plot.data$column_label_noise)
# # unique(plot.data.bounds$column_label_noise)
# 
# list.plots <- list()
# for (taxon in taxa.for.ice) {
#     # taxon <- taxa.for.ice[5]
#     cat("\nProduce ICE plot for:", taxon)
#     occ.taxon <- paste0("Occurrence.", taxon)
# 
#     plot.data.bounds.taxon <- plot.data.bounds %>%
#         unnest(y.mean.min) %>%
#         rename(y.mean.min = all_of(taxon))  %>%
#         # select(-all_of(taxa.for.ice)) %>%
#         unnest(y.mean.max, names_sep = "") %>%
#         rename(y.mean.max = all_of(paste0("y.mean.max", taxon)))
# 
#         
#     fig3 <- ggplot(data=plot.data) +
#         geom_line(aes(x=.data[[select.env.fact]],
#                       y=.data[[taxon]],
#                       group=ReachID, 
#                       color=ReachID),
#                   show.legend = FALSE) +
#         geom_line(data=plot.data.mean,
#                   aes(x=.data[[name.select.env.fact]], y=.data[[taxon]]),
#                   linewidth=1.2) +
#         geom_line(data=plot.data.stream[,c( "avg", #taxon, 
#                                            name.select.env.fact)],
#                   aes(x=.data[[name.select.env.fact]], y=.data[["avg"]]), # y=.data[[taxon]]),
#                   linewidth=1.2, color = "grey50", linetype="dashed") + #,
#                   # alpha = 0.6, inherit.aes=F) +
#         geom_rug(data = env.factor.sampled,
#                  aes(x=variable), 
#                  color="grey20",
#                  alpha=0.7,
#                  inherit.aes=F) + 
#         # scale_x_continuous(limits = c(min(env.factor.sampled), max(env.factor.sampled))) +
#         geom_segment(data=plot.data.bounds.taxon,
#                      inherit.aes = FALSE,
#                      lineend="round",
#                      linejoin="round",
#                      aes(x=x.mean,
#                          y=y.mean.min,
#                          xend=x.mean,
#                          yend=y.mean.max),
#                      arrow=arrow(length = unit(0.3, "cm"),
#                                  ends = "both")) +
#         # facet_grid(column_label_noise ~ model, labeller = label_wrap_gen()) +
#         facet_grid(model ~ column_label_noise, labeller = label_wrap_gen()) +
#         
#         labs(title = "",
#             # title = taxon,
#              x = name.select.env.fact,
#              y = "Predicted probability of occurrence")
#     # fig3
#     
#     list.plots[[taxon]] <- fig3
#     file.name.tax <- paste0(file.name.exp,
#                             taxon, "_",
#                             select.env.fact)
#     cat("\nSaving:", file.name.tax)
#     
#     pdf(paste0(dir.compar.ice, file.name.tax, "_ice.pdf"), width = width.a4*1.2, height = height.a4*1.2)
#     print(fig3)
#     dev.off()
#     
#     ggsave(paste0(dir.compar.ice, file.name.tax, "_ice.png"), 
#            width = 2480*1.15,
#            height = 3508*1.15,
#            units = c("px"))
# }    
# # }
# 
# file.name <- "ice_all_taxa"
# print.pdf.plots(list.plots = list.plots, width = width.a4*1.2, height = height.a4*1.2,
#                 dir.output = dir.compar.plots, 
#                 info.file.name = paste0(file.name.exp, length(taxa.for.ice), "taxa_", select.env.fact, "_"),
#                 file.name = file.name)
# 
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Fig 1 & Fig 2: PDP comparison ----
# 
# dir.compar.pdp <- paste0(dir.compar.plots, "pdp/")
# dir.create(dir.compar.pdp)
# 
# # compute boundaries
# min.boundaries <- lapply(multi.mean, FUN=function(ice){
#   min(ice[[name.select.env.fact]])
# })
# 
# max.boundaries <- lapply(multi.mean, FUN=function(ice){
#   max(ice[[name.select.env.fact]])
# })
# 
# lb <- max(unlist(min.boundaries))
# hb <- min(unlist(max.boundaries))
# 
# plot.data <- final.multi.ice
# plot.data.stream <- mean.ice.streamb %>%
#     filter(Temperature < max.temp,
#            Temperature > min.temp) # ! predictor specific
# 
# list.plots.pdp <- list()
# for (taxon in taxa.for.ice) {
#     # taxon <- taxa.for.ice[1]
#     cat("\nProduce PDP plot for:", taxon)
#     occ.taxon <- paste0("Occurrence.", taxon)
#         
#     fig1 <- ggplot(data=plot.data) +
#       geom_line(aes(x=.data[[name.select.env.fact]],
#                     y=.data[[taxon]],
#                     group=model,
#                     colour=model),
#                 size=1, alpha = 0.8) +
#       geom_line(data=plot.data.stream,
#                 aes(x=.data[[name.select.env.fact]], y=.data[[taxon]]),
#                 size=1, color = "darkgrey", alpha = 0.8) +
#       facet_wrap(~factor(column_label, levels=unlist(names(list.exp))), ncol = 3,
#                  labeller = label_wrap_gen())+
#       #facet_wrap(~column_label) +
#       xlim(lb, hb) +
#       scale_y_continuous(limits = c(0,1)) +
#       scale_color_manual(values=model.color.map) +
#         # theme(axis.text.x = element_text(hjust = -1)) +
#       labs(x = name.select.env.fact,
#            y = "Predicted probability of occurrence",
#            colour="Models",
#            title = taxon)
#     # fig1
#     
#     list.plots.pdp[[taxon]] <- fig1
#     file.name.tax <- paste0(file.name.exp,
#                             taxon, "_",
#                             select.env.fact)
#     cat("\nSaving:", file.name.tax)
#     
#     pdf(paste0(dir.compar.pdp, file.name.tax, "_pdp.pdf"), width = width.a4*1.1, height = height.a4*0.6)
#     print(fig1)
#     dev.off()
#     
#     ggsave(paste0(dir.compar.pdp, file.name.tax, "_pdp.png"), width = 2480*1.1,
#            height = 3508*0.6,
#            units = c("px"))
# }    
# 
# file.name <- "pdp_all_taxa"
# print.pdf.plots(list.plots = list.plots.pdp, width = width.a4*1.1, height = height.a4*0.6,
#                 dir.output = dir.compar.plots, 
#                 info.file.name = paste0(file.name.exp, length(taxa.for.ice), "taxa_", select.env.fact, "_"),
#                 file.name = file.name)
# 
# # rm(fig1)
# 
# # # comparison pdp coloured scenarios
# # fig2 <- ggplot(data=plot.data) +
# #   geom_line(aes(x=.data[[name.select.env.fact]],
# #                 y=.data[[taxon]],
# #                 group=column_label,
# #                 colour=column_label),
# #             size=1) +
# #   geom_line(data=plot.data.stream,
# #             aes(x=.data[[name.select.env.fact]], y=.data[[taxon]]),
# #             size=1, color = "darkgrey", alpha = 0.8) +
# #   scale_color_manual(values=color.map) +
# #   xlim(lb, hb) +
# #   # scale_y_continuous(limits = c(0,1)) +
# #   facet_wrap(~model) +
# #   # facet_wrap(~factor(model, levels = c("GLM", "GAM", "RF"))) +
# #   # theme_bw() +
# #   # theme(strip.background = element_rect(fill = "white"),
# #   #       legend.position = "bottom") +#,
# #   # legend.title = element_text(size=24),
# #   # legend.text = element_text(size=20)) +
# #   labs(x = name.select.env.fact,
# #        y = "Predicted probability of occurrence",
# #        colour="Scenario", 
# #        title = taxon)
# # # fig2
# # 
# # pdf(paste0(dir.compar.plots, file.name.tax, "_pdp_per_model.pdf"), width = width.a4, height = height.a4*2/3)
# # print(fig2)
# # dev.off()
# 
# # rm(fig2)
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Fig 3: plotting the score of given taxon ####
# 
# dir.compar.perf <- paste0(dir.compar.plots, "perf/")
# dir.create(dir.compar.perf)
# 
# list.plots.score <- list()
# for (taxon in names(taxa.colnames)) {
#     # taxon <- names(taxa.colnames)[1]
#     # filter results for selected taxon
#     plot.data <- final.multi.all.results %>%
#         filter(taxa == taxon)
#     plot.data$shape <- ifelse(plot.data$fit_pred == "fit", "Calibration", "Prediction")
#     plot.data$shape <- factor(plot.data$shape, levels= c("Calibration", "Prediction"))
#     
#     fig3 <- ggplot(data=plot.data) +
#         geom_point(aes(x=column_label,
#                        y=dev,
#                        shape=shape, 
#                        color=model),
#                    size=3,
#                    alpha=0.7) +
#         scale_x_discrete(limits=names(list.exp)) +
#         scale_color_manual(values=model.color.map) +
#         labs(x="Scenario",
#              y="Standardized deviance",
#              title = taxon,
#              color = "Model",
#              shape = "") +
#         theme(axis.text.x = element_text(angle = 45, #vjust = 0.5, 
#                                          hjust=1))
#     list.plots.score[[taxon]] <- fig3
#     file.name.tax <- paste0(file.name.exp,
#                             taxon)
#     cat("\nSaving:", file.name.tax)
#     
#     pdf(paste0(dir.compar.perf, file.name.tax, "_perf.pdf"), width = width.a4*1.1, height = height.a4*0.6)
#     print(fig3)
#     dev.off()
#     
#     ggsave(paste0(dir.compar.perf, file.name.tax, "_perf.png"), width = 2480*1.1,
#            height = 3508*0.6,
#            units = c("px"))
# }    
# 
# file.name <- "perf_all_taxa"
# print.pdf.plots(list.plots = list.plots.score, width = width.a4*1.1, height = height.a4*0.6,
#                 dir.output = dir.compar.plots, 
#                 info.file.name = paste0(file.name.exp, length(taxa.for.ice), "taxa_"),
#                 file.name = file.name)
# 
# # print pdf with pdp and performance plots combined in one page
# pdf(paste0(dir.compar.plots, paste0(file.name.exp, length(taxa.for.ice), "taxa"), "_perf_pdp.pdf"), width = width.a4*1.1, height = height.a4*1.1)
# for(taxon in taxa.for.ice){
#     # taxon <- taxa.for.ice[1]
#     # plot <- ggarrange(list.plots.pdp[[taxon]], ggarrange(list.plots.score[[taxon]], list.plots.stream[[taxon]], ncol = 2), ncol = 1)
#     plot <- ggarrange(list.plots.pdp[[taxon]], list.plots.score[[taxon]], ncol = 1)
#     annotate_figure(plot, top = text_grob(taxon, 
#                                           color = "black", face = "bold", size = 16))
#     print(plot)
# }
# dev.off()
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Fig 6: multi-scenario bellplots ----
# 
# plot.data <- final.multi.all.results %>%
#     filter(fit_pred == "pred")
# 
# fig6 <- ggplot(data=plot.data,
#                aes(x=prevalence, y=dev, color=column_label)) + 
#     geom_point() +
#     facet_wrap(~model)  +
#     scale_color_manual(values=color.map) + 
#     scale_y_continuous(limits = c(0,2)) +
#     labs(x = "Prevalence",
#          y = "Standardized deviance") +
#     # scale_x_discrete(names(taxa.colnames)) +
#     theme_minimal() + 
#     theme(legend.title=element_blank())
# # fig6
# 
# pdf(paste0(dir.compar.plots, file.name.exp, "_bellplot_pred_per_scenario.pdf"), height = 6, width = 10)
# print(fig6)
# dev.off() 
# 
# fig7 <- ggplot(data=plot.data,
#                aes(x=prevalence, y=dev, color=model)) + 
#     geom_point(alpha = 0.6) +
#     facet_wrap(~column_label, labeller = label_wrap_gen())  +
#     scale_color_manual(values=model.color.map) + 
#     scale_y_continuous(limits = c(0,1.5)) +
#     labs(x = "Prevalence",
#          y = "Standardized deviance",
#          color = "Model")# +
#     # scale_x_discrete(names(taxa.colnames)) +
#     # theme_minimal() + 
#     # theme(legend.title=element_blank())
# 
# # fig7
# 
# pdf(paste0(dir.compar.plots, file.name.exp, "_bellplot_pred_per_model.pdf"), height = height.a4*2/3, width = width.a4)
# print(fig7)
# dev.off() 
# 
# ggsave(paste0(dir.compar.plots, file.name.exp, "bellplot_pred_per_model.png"), width = 2700,
#        height = 2300,
#        units = c("px"))
# 
# # Fig X: box plot per model ####
# 
# fig6 <- ggplot(data=plot.data) +
#     geom_boxplot(aes(x=model,
#                      y=dev,
#                      fill=column_label)) +
#     scale_x_discrete(limits=names(models))+
#     scale_y_continuous(limits = c(0, 2)) +
#     facet_wrap(~fit_pred, nrow = 2) + 
#     theme_minimal() +
#     scale_fill_manual(values = color.map) +
#     theme(legend.title=element_blank()) +
#     labs(x="Model",
#          y="Standardized deviance")
# # fig6
# 
# 
# pdf(paste0(dir.compar.plots, file.name.general, "_boxplot_per_model.pdf"), width = 12, height = 8)
# print(fig6)
# dev.off() 
# 
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Fig 4: dual ICE ----
# # exp.baseline <- list.exp[[1]]
# # # exp.extreme  <- list.exp[[length(list.exp)]]
# # no.select.exp <- 3
# # exp.extreme  <- list.exp[[no.select.exp]]
# 
# list.dual.exp <- list.exp[c(1,3)]
# 
# # taxon.under.obs <- paste0("Occurrence.", taxon.under.obs)
# 
# # list.dual.exp <- list("Baseline" = exp.baseline,
# #                       names(list.exp)[no.select.exp] = exp.extreme)
# 
# dual.ice <- lapply(list.dual.exp, FUN=function(name){
#   
#     # name <- list.dual.exp[[1]]
#   dir.experiment          <- paste0(dir.output, name, "/")
#   print(dir.experiment)
#   models.fit      <- load.models(path=dir.experiment, split.type="FIT")
#   std.const.fit   <- readRDS(file=paste0(dir.experiment, "standardization_constant_FIT.rds"))
#   correct.names <- names(models)[which(names(models.fit) %in% models | names(models.fit) %in% names(models))]
#   names(models.fit) <- correct.names
#   
#   # input.env.factors <- get.input.env.factors(exp.name=name)
#   input.env.factors <- input.env.factors
#   
#   ice.dfs <- plot.ice(models.performance=models.fit,
#                       select.env.fact=select.env.fact,
#                       taxa=taxon.under.obs,
#                       standardization.constant=std.const.fit[[1]],
#                       observations = preprocessed.data.fit$`Entire dataset`,
#                       nb.sample=no.sites,
#                       resolution=no.steps,
#                       input.env.factors=input.env.factors)
#   
#   return(ice.dfs)
# })
# 
# obs1      <- dual.ice[[1]][["observations"]]
# env.fact1 <- dual.ice[[1]][["env.factor.sampled"]]
# obs2      <- dual.ice[[2]][["observations"]]
# env.fact2 <- dual.ice[[2]][["env.factor.sampled"]]
# 
# merged.obs       <- bind_rows(list(obs1, obs2),
#                               .id = "column_label")
# merged.env.fact  <- bind_rows(list(env.fact1, env.fact2),
#                               .id = "column_label")
# 
# observations.mean         <- merged.obs %>%
#   group_by(across(all_of(select.env.fact)), model, column_label) %>%
#   summarise(avg = mean(pred))
# observations.mean.bounds  <- observations.mean %>% group_by(model, column_label) %>%
#   summarise(x.mean=max(across(all_of(name.select.env.fact))),
#             y.mean.min=min(avg),
#             y.mean.max=max(avg))
# 
# list.obs.mean <- list(observations.mean %>% filter(column_label==1),
#                       observations.mean %>% filter(column_label==2))
# 
# min.boundaries <- lapply(list.obs.mean, FUN=function(ice){
#   min(ice[[name.select.env.fact]])
# })
# 
# max.boundaries <- lapply(list.obs.mean, FUN=function(ice){
#   max(ice[[name.select.env.fact]])
# })
# 
# lb <- max(unlist(min.boundaries))
# hb <- min(unlist(max.boundaries))
# 
# names.scenarios <- names(list.dual.exp)
# 
# merged.obs["column_label"] <- ifelse(merged.obs[["column_label"]]==1, names.scenarios[1], names.scenarios[2]) 
# observations.mean["column_label"] <- ifelse(observations.mean[["column_label"]]==1, names.scenarios[1], names.scenarios[2])
# merged.env.fact["column_label"] <- ifelse(merged.env.fact[["column_label"]]==1, names.scenarios[1], names.scenarios[2])
# observations.mean.bounds["column_label"] <- ifelse(observations.mean.bounds[["column_label"]]==1, names.scenarios[1], names.scenarios[2])
# 
# plot.data <- merged.obs
# plot.data$model <- factor(plot.data$model, levels= names(models))
# 
# fig4 <- ggplot(data=plot.data) +
#   geom_line(aes(x=.data[[select.env.fact]],
#                 y=pred,
#                 group=observation_number, 
#                 color=as.character(observation_number)),
#             alpha=0.4,
#             show.legend = FALSE) +
#   geom_line(data=observations.mean,
#             aes(x=.data[[name.select.env.fact]], y=avg),
#             size=1.5) +
#   geom_rug(data = merged.env.fact,
#            aes(x=variable),
#            alpha=0.7,
#            inherit.aes=F) + 
#   geom_segment(data=observations.mean.bounds,
#                inherit.aes = FALSE,
#                lineend="round",
#                linejoin="round",
#                aes(x=hb,
#                    y=y.mean.min,
#                    xend=hb,
#                    yend=y.mean.max),
#                arrow=arrow(length = unit(0.3, "cm"),
#                            ends = "both"),
#                alpha=0.9) +
#   xlim(lb, hb) +
#   facet_wrap(~column_label+model) +
#   theme_bw() +
#   theme(strip.background = element_rect(fill = "white")) +#,
#   # legend.title = element_text(size=24),
#   # legend.text = element_text(size=20)) +
#   labs(x =name.select.env.fact,
#        y = "Predicted probability of occurrence", 
#        title = taxon.under.obs)
# # fig4
# 
# pdf(paste0(dir.compar.plots, filename, "_ice.pdf"))
# print(fig4)
# dev.off()
# 
# rm(fig4)
# 
# # repeat fig4 for each model separately
# for (model.name in c("GLM", "GAM", "RF")){ # c("glm", "gamloess", "rf", "ann")){
#   
#   filtered.merged.obs               <- merged.obs               %>% filter(model == model.name)
#   filtered.observations.mean        <- observations.mean        %>% filter(model == model.name)
#   filtered.observations.mean.bounds <- observations.mean.bounds %>% filter(model == model.name)
#   
#   fig <- ggplot(data=filtered.merged.obs) +
#     geom_line(aes(x=.data[[select.env.fact]],
#                   y=pred,
#                   group=observation_number, 
#                   color=as.character(observation_number)),
#               show.legend = FALSE) +
#     geom_line(data=filtered.observations.mean,
#               aes(x=.data[[name.select.env.fact]], y=avg),
#               size=1.5) +
#     geom_rug(data = merged.env.fact,
#              aes(x=variable), 
#              color="grey20",
#              alpha=0.7,
#              inherit.aes=F) + 
#     geom_segment(data=filtered.observations.mean.bounds,
#                  inherit.aes = FALSE,
#                  lineend="round",
#                  linejoin="round",
#                  aes(x=hb,
#                      y=y.mean.min,
#                      xend=hb,
#                      yend=y.mean.max),
#                  arrow=arrow(length = unit(0.3, "cm"),
#                              ends = "both")) +
#     xlim(lb, hb) +
#     facet_wrap(~column_label) +
#     theme_bw() +
#     theme(strip.background = element_rect(fill = "white")) +#,
#     # legend.title = element_text(size=24),
#     # legend.text = element_text(size=20)) +
#     labs(x =name.select.env.fact,
#          y = "Predicted probability of occurrence")
#   
#   # dir <- "../output_data/comparison_plots/"
#   
#   pdf(paste0(dir.compar.plots, filename, "_ice_", model.name, ".pdf"), width = 8, height = 4)
#   print(fig)
#   dev.off()
#   
#   rm(fig)
#   
# }
# 
# # # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # # Fig 5: multi box plot ---
# # multi.all.results <- lapply(list.exp, FUN=function(name){
# #   
# #   dir.experiment <- paste0(dir.output, name, "/")
# #   models.cv      <- load.models(path=dir.experiment, split.type="CV")
# #   
# #   all.results    <- summarize.all.results(models.cv, data.prev.taxa)
# #   
# #   rm(models.cv) 
# #   
# #   all.results <- restructure.all.results(all.results, models)
# #   
# #   all.results["noise"] <- name
# #   
# #   return(all.results)
# # })
# # 
# # final.multi.all.results <- bind_rows(multi.all.results, .id = "column_label")
# # 
# # fig5 <- ggplot(data=final.multi.all.results) +
# #   geom_boxplot(aes(x=column_label,
# #                    y=dev,
# #                    fill=fit_pred)) +
# #   scale_x_discrete(limits=names(list.exp)) +
# #   facet_wrap(~model) + 
# #   theme_minimal() +
# #   theme(legend.title=element_blank())
# # 
# # 
# # pdf(paste0(dir.compar.plots, filename, "_boxplot.pdf"))
# # print(fig5)
# # dev.off() 
# # 
# # rm(fig5)
# 
# 
# # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# # Fig 6: multi-scenario bellplots ----
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
# 
# fig6 <- ggplot(data=final.multi.all.results,
#                aes(x=prevalence, y=dev, color=column_label)) + 
#   geom_point(data=final.multi.all.results) +
#   facet_wrap(~model + fit_pred)  +
#   scale_color_manual(values=color.map) + 
#   theme_minimal() + 
#   theme(legend.title=element_blank())
# 
# pdf(paste0(dir.compar.plots, filename, "_bellplot.pdf"))
# print(fig6)
# dev.off() 
# 
# rm(fig6)
# 
