


create.plots <- function(experiment.name="experiment_name",
                         file.prev.taxa="All_2729samples_9envfact_lme.area.elev_PrevalenceTaxa.csv"){
  
  
  # The goal of this function is to create the different plots assessing 
  # individually the performance of an experiment (i.e. one model trained under
  # specific noise conditions)
  #
  # arguments:
  #   - experiment.name:        The name of the folder (located in /output_data)
  #                             containing all the data of the experiment we 
  #                             wish to plot
  #   - file.prev.taxa:         The name of the file containing the information 
  #                             about the prevalence of taxa for instance
  #
  # returns:
  #   - nothing is returned but plots are saved as .pdf files
  
  
  # Libraries ----
  if (!require("dplyr")){install.packages("dplyr"); library("dplyr")}
  if (!require("readr")){install.packages("readr"); library("readr")}
  
  # to compute AUC
  if (!require("pROC")){install.packages("pROC"); library("pROC")}                       
  
  # For Neural Network  
  if (!require("reticulate")){install.packages("reticulate"); library("reticulate")}
  #install_miniconda()              # run this the very first time reticulate is installed
  #install.packages("tensorflow")
  library("tensorflow")
  #install_tensorflow()             # run this line only when opening new R session
  #install.packages("keras")
  library("keras")
  #install_keras()                  # run this line only when opening new R session
  #use_condaenv()
  
  # caret has to be loaded at the end to not cache function 'train'
  if (!require("caret")){install.packages("caret"); library("caret")}
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Custom R scripts ---- 
  source("global_variables.r")
  source("ml_models.r")
  source("performances_assessment.r")
  source("dataset_preprocessing.r")
  
  # load the the files required for plotting ----
  dir.input.data  <- "../input_data/"
  prev.taxa       <- read.csv(paste0(dir.input.data, file.prev.taxa),
                              header = T, 
                              sep = ";", 
                              stringsAsFactors = F)
  prev.taxa <- filter(prev.taxa, Occurrence.taxa %in% TAXA.COLNAMES)
  
  
  dir.output              <- "../output_data/"
  dir.experiment          <- paste0(dir.output, experiment.name, "/")
  
  
  models.fit      <- load.models(path=dir.experiment,
                                 split.type="FIT")
  models.cv       <- load.models(path=dir.experiment, 
                                 split.type="CV")
  std.const.cv    <- readRDS(file=paste0(dir.experiment, "standardization_constant_CV.rds"))
  std.const.fit   <- readRDS(file=paste0(dir.experiment, "standardization_constant_FIT.rds"))
  
  # Variable for plots ----
  all.results <- summarize.all.results(models.cv, prev.taxa)
  
  ggplot.all.results <- restructure.all.results(all.results)
  
  color.map <- c('glm'           = 'deepskyblue',   # Generalized Linear Model
                 'gamloess'      = 'green',         # Generalized Additive Model
                 'ann'           = 'orange',        # Artificial Neural Network
                 'rf'            = 'red',           # Random Forest
                 'null'          = 'black')         # Null Model
  
  
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
    geom_hline(aes(yintercept=median.null.pred,
                   color="median null model (pred)"),
               linetype="dashed") + 
    geom_hline(aes(yintercept=median.best.pred,
                   color="median best model (pred)"),
               linetype="dashed") +
    scale_colour_manual(values = c('green','black')) +
    theme_minimal() +
    theme(legend.title=element_blank())
  
  
  # Fig. 2 ----
  fig2 <- ggplot(data=subset(ggplot.all.results, model %in% "null"),
                 aes(x=prevalence, y=dev, color=model)) +
    geom_smooth(se = FALSE) + 
    geom_point(data=ggplot.all.results) +
    facet_wrap(~fit_pred) +
    scale_color_manual(values=color.map) +
    theme_minimal() + 
    theme(legend.title=element_blank())
  
  
  # Fig. 3 ----
  # env.fact.under.obs <- "temperature"
  env.fact.under.obs <- "tempmaxC"
  
  taxon.under.obs <- "Occurrence.Gammaridae"
  
  input.env.factors <- get.input.env.factors(experiment.name)
  
  ice.dfs <- plot.ice(models.performance=models.fit,
                      env.factor=env.fact.under.obs,
                      taxa=taxon.under.obs,
                      standardization.constant=std.const.fit[[1]],
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
  
  
  # Fig 4
  noised.dataset <- models.fit[[1]][[1]][[1]][[1]][["observation"]]
  taxa.col <- colnames(noised.dataset)[which(grepl("Occurrence.", colnames(noised.dataset)))]
  
  noised.prevalence <- lapply(taxa.col, FUN=function(taxon){
    occurence  <- noised.dataset[[taxon]]
    prevalence <- length(occurence[occurence=="present"])/length(occurence) 
    return(prevalence)
  })  
  
  noised.prev.taxa <- prev.taxa[which(prev.taxa$Occurrence.taxa %in% TAXA.COLNAMES),]  
  noised.prev.taxa[["Prevalence"]] <- noised.prevalence
  
  prev.taxa[["Status"]] <- "unnoised"
  noised.prev.taxa[["Status"]] <- "noised"
  
  full.prev.taxa <- rbind(prev.taxa, noised.prev.taxa)
  full.prev.taxa[["Occurrence.taxa"]] <- as.factor(full.prev.taxa[["Occurrence.taxa"]])
  full.prev.taxa$Prevalence <- unlist(full.prev.taxa$Prevalence)
  full.prev.taxa$Occurrence.taxa <- gsub("Occurrence.", "", full.prev.taxa$Occurrence.taxa)
  
  fig4 <- ggplot(data=full.prev.taxa,
                 aes(x=Prevalence,
                     y=Occurrence.taxa,
                     color=Status)) + 
          geom_point() +
          scale_y_discrete(limits = rev(gsub("Occurrence.", "", prev.taxa$Occurrence.taxa))) +
          labs(x="Prevalence", y="Taxa")
          theme_minimal() + 
          theme(legend.title=element_blank())
  
  
  # Saving figures to experiment folder ----
  pdf(paste0(dir.experiment, "fig1.pdf"))
  print(fig1)
  dev.off()
  
  pdf(paste0(dir.experiment, "fig2.pdf"))
  print(fig2)
  dev.off() 
  
  pdf(paste0(dir.experiment, "fig3.pdf"))
  print(fig3)
  dev.off()
  
  pdf(paste0(dir.experiment, "fig4.pdf"))
  print(fig4)
  dev.off()
  
}

create.comparison.plots <- function(filename,
                                    list.exp,
                                    file.prev.taxa="All_2729samples_9envfact_lme.area.elev_PrevalenceTaxa.csv",
                                    taxon="Gammaridae",
                                    env.fact="temperature"){
  
  # The goal of this function is to create plots comparing the performance of a 
  # different experiment (i.e. models trained under different noise conditions)
  #
  # arguments:
  #   - filename:               This string will be taken as prefix when naming 
  #                             and saving the different plot in the |plots 
  #                             folder
  #   - list.exp:               List of names of the folder (located in
  #                             /output_data) containing all the data of the
  #                             experiment we wish to compare
  #   - file.prev.taxa:         The name of the file containing the information 
  #                             about the prevalence of taxa for instance
  #   - taxon:                  Which taxon will be considered when plotting the
  #                             ICE
  #   - env.fact                Which environment factor will be considered when
  #                             plotting the ICE
  #
  # returns:
  #   - nothing is returned but plots are saved as .pdf files

  
  env.fact.string <- env.fact
  
  if (length(names(ENV.FACT.MAP)[which(env.fact %in% ENV.FACT.MAP)])!=0){
    env.fact.string <- names(ENV.FACT.MAP)[which(env.fact %in% ENV.FACT.MAP)]
  }
  
  
  dir <- "../plots/comparison_plots/"
  dir.create(dir)
  dir.input       <- "../input_data/"
  dir.output      <- "../output_data/"
  prev.taxa       <- read.csv(paste0(dir.input, file.prev.taxa),
                              header = T, 
                              sep = ";", 
                              stringsAsFactors = F)
  
  color.map <- c('1'            = 'deepskyblue',
                 '2'            = 'green',
                 '3'            = 'orange',
                 '4'            = 'red')
  
  if (length(color.map)==length(list.exp)){
    names(color.map) <- names(list.exp)
  }
  
  model.color.map <- c('glm'     = 'deepskyblue',   # Generalized Linear Model
                 'gamloess'      = 'green',         # Generalized Additive Model
                 'ann'           = 'orange',        # Artificial Neural Network
                 'rf'            = 'red',           # Random Forest
                 'null'          = 'black')         # Null Model
  
  
  
  # Fig 1 & Fig 2: PDP comparison 
  taxon.under.obs <- paste0("Occurrence.", taxon)
  
  multi.ice <- lapply(list.exp, FUN=function(name){
    
    dir.experiment          <- paste0(dir.output, name, "/")
    
    models.fit      <- load.models(path=dir.experiment, split.type="FIT")
    std.const.fit   <- readRDS(file=paste0(dir.experiment, "standardization_constant_FIT.rds"))
    
    input.env.factors <- get.input.env.factors(exp.name=name)
    
    ice.dfs <- plot.ice(models.performance=models.fit,
                        env.factor=env.fact,
                        taxa=taxon.under.obs,
                        standardization.constant=std.const.fit[[1]],
                        nb.sample=100,
                        resolution=200,
                        input.env.factors=input.env.factors)
    
    observations              <- ice.dfs[["observations"]]
    observations.mean         <- observations %>%
      group_by(across(all_of(env.fact)), model) %>%
      summarise(avg = mean(pred))
    
    observations.mean["noise"] <- name
    
    return(observations.mean)
  })
  
  final.multi.ice <- bind_rows(multi.ice, .id = "column_label")
  
  # compute boundaries
  min.boundaries <- lapply(multi.ice, FUN=function(ice){
    min(ice[[env.fact]])
  })
  
  max.boundaries <- lapply(multi.ice, FUN=function(ice){
    max(ice[[env.fact]])
  })
  
  lb <- max(unlist(min.boundaries))
  hb <- min(unlist(max.boundaries))
  
  fig1 <- ggplot(data=final.multi.ice) +
          geom_line(aes(x=.data[[env.fact]],
                        y=avg,
                        group=model,
                        colour=model),
                    size=1) +
          facet_wrap(~factor(column_label, levels=unlist(names(list.exp))))+
          #facet_wrap(~column_label) +
          xlim(lb, hb) +
          labs(x =env.fact.string,
               y = "Predicted probability of occurrence",
               colour="Models")
  
  
  pdf(paste0(dir, filename, "_pdp_per_scenario.pdf"))
  print(fig1)
  dev.off()
  
  rm(fig1)
  
  fig2 <- ggplot(data=final.multi.ice) +
          geom_line(aes(x=.data[[env.fact]],
                        y=avg,
                        group=column_label,
                        colour=column_label),
                    size=1) +
          scale_color_manual(values=color.map) +
          xlim(lb, hb) +
          facet_wrap(~model) +
          labs(x =env.fact.string,
               y = "Predicted probability of occurrence",
               colour="Scenario")
  
  pdf(paste0(dir, filename, "_pdp_per_model.pdf"))
  print(fig2)
  dev.off()
  
  rm(fig2)
  
  
  # Fig 3: plotting the score of given taxon
  multi.all.results <- lapply(list.exp, FUN=function(name){
    
    dir.experiment <- paste0(dir.output, name, "/")
    models.cv      <- load.models(path=dir.experiment, split.type="CV")
    
    all.results    <- summarize.all.results(models.cv, prev.taxa)
    
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
  
  
  pdf(paste0(dir, filename, "_taxa_score.pdf"))
  print(fig3)
  dev.off()
  
  rm(fig3)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Fig 4: dual ICE ----
  exp.baseline <- list.exp[[1]]
  exp.extreme  <- list.exp[[length(list.exp)]]
  
  taxon.under.obs <- paste0("Occurrence.", taxon)
  
  list.dual.exp <- list("Baseline" = exp.baseline,
                        "High noise"  = exp.extreme)
  
  dual.ice <- lapply(list.dual.exp, FUN=function(name){
    
    dir.experiment          <- paste0(dir.output, name, "/")
    
    models.fit      <- load.models(path=dir.experiment, split.type="FIT")
    std.const.fit   <- readRDS(file=paste0(dir.experiment, "standardization_constant_FIT.rds"))
    
    input.env.factors <- get.input.env.factors(exp.name=name)
    
    ice.dfs <- plot.ice(models.performance=models.fit,
                        env.factor=env.fact,
                        taxa=taxon.under.obs,
                        standardization.constant=std.const.fit[[1]],
                        nb.sample=100,
                        resolution=200,
                        input.env.factors=input.env.factors)
    
    return(ice.dfs)
  })
  
  obs1      <- dual.ice[["Baseline"]][["observations"]]
  env.fact1 <- dual.ice[["Baseline"]][["env.factor.sampled"]]
  obs2      <- dual.ice[["High noise"]][["observations"]]
  env.fact2 <- dual.ice[["High noise"]][["env.factor.sampled"]]
  
  merged.obs       <- bind_rows(list(obs1, obs2),
                         .id = "column_label")
  merged.env.fact  <- bind_rows(list(env.fact1, env.fact2),
                         .id = "column_label")
  
  observations.mean         <- merged.obs %>%
    group_by(across(all_of(env.fact)), model, column_label) %>%
    summarise(avg = mean(pred))
  observations.mean.bounds  <- observations.mean %>% group_by(model, column_label) %>%
    summarise(x.mean=max(across(all_of(env.fact))),
              y.mean.min=min(avg),
              y.mean.max=max(avg))
  
  list.obs.mean <- list(observations.mean %>% filter(column_label==1),
                        observations.mean %>% filter(column_label==2))
  
  min.boundaries <- lapply(list.obs.mean, FUN=function(ice){
    min(ice[[env.fact]])
  })
  
  max.boundaries <- lapply(list.obs.mean, FUN=function(ice){
    max(ice[[env.fact]])
  })
  
  lb <- max(unlist(min.boundaries))
  hb <- min(unlist(max.boundaries))
  
  merged.obs["column_label"] <- ifelse(merged.obs[["column_label"]]==1, "Baseline", "High noise") 
  observations.mean["column_label"] <- ifelse(observations.mean[["column_label"]]==1, "Baseline", "High noise")
  merged.env.fact["column_label"] <- ifelse(merged.env.fact[["column_label"]]==1, "Baseline", "High noise")
  observations.mean.bounds["column_label"] <- ifelse(observations.mean.bounds[["column_label"]]==1, "Baseline", "High noise")
  
  
  fig4 <- ggplot(data=merged.obs) +
    geom_line(aes(x=.data[[env.fact]],
                  y=pred,
                  group=observation_number, 
                  color=as.character(observation_number)),
              alpha=0.4,
              show.legend = FALSE) +
    geom_line(data=observations.mean,
              aes(x=.data[[env.fact]], y=avg),
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
    labs(x =env.fact.string,
         y = "Predicted probability of occurrence")
  
  pdf(paste0(dir, filename, "_ice.pdf"))
  print(fig4)
  dev.off()
  
  rm(fig4)
  
  # repeat fig4 for each model separately
  for (model.name in c("glm", "gamloess", "rf", "ann")){
    
    filtered.merged.obs               <- merged.obs               %>% filter(model == model.name)
    filtered.observations.mean        <- observations.mean        %>% filter(model == model.name)
    filtered.observations.mean.bounds <- observations.mean.bounds %>% filter(model == model.name)
    
    fig <- ggplot(data=filtered.merged.obs) +
      geom_line(aes(x=.data[[env.fact]],
                    y=pred,
                    group=observation_number, 
                    color=as.character(observation_number)),
                show.legend = FALSE) +
      geom_line(data=filtered.observations.mean,
                aes(x=.data[[env.fact]], y=avg),
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
      labs(x =env.fact.string,
           y = "Predicted probability of occurrence")
    
    # dir <- "../output_data/comparison_plots/"
    
    pdf(paste0(dir, filename, "_ice_", model.name, ".pdf"))
    print(fig)
    dev.off()
    
    rm(fig)
    
  }
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Fig 5: multi box plot ----
  multi.all.results <- lapply(list.exp, FUN=function(name){
  
    dir.experiment <- paste0(dir.output, name, "/")
    models.cv      <- load.models(path=dir.experiment, split.type="CV")
  
    all.results    <- summarize.all.results(models.cv, prev.taxa)
  
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
  
  
  pdf(paste0(dir, filename, "_boxplot.pdf"))
  print(fig5)
  dev.off() 
  
  rm(fig5)
  
  
  # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  # Fig 6: multi-scenario bellplots ----
  multi.all.results <- lapply(list.exp, FUN=function(name){
  
    dir.experiment <- paste0(dir.output, name, "/")
    models.cv      <- load.models(path=dir.experiment, split.type="CV")
  
    all.results    <- summarize.all.results(models.cv, prev.taxa)
  
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
  
  pdf(paste0(dir, filename, "_bellplot.pdf"))
  print(fig6)
  dev.off() 
  
  rm(fig6)
}


get.input.env.factors <- function(exp.name){
  
  # this function takes as input as experiment and returns the env.factor that 
  # were used when training the model. To do so, it gets the env. factors found
  # in the metadata of the experiment and modify them according to the noises
  # that were used in this experiment
  
  # arguments:
  #   - exp.name:               The experiment for which it is desired to know 
  #                             what environmental factors were used when 
  #                             training the models
  #
  # returns:
  #   - env.factor:             The name of the environmental factors used to
  #                             train the models in the experiment given as the 
  #                             input
  
  dir.metadata <- paste0("../output_data/", exp.name, "metadata.json")
  metadata <- read_json(dir.metadata)
  
  # 1 get input.env.factor from metadata
  env.factors <- unlist(metadata[["env_factor"]])
  # 2 get noise from metadata
  noise <- metadata[["noise"]]
  
  # 3 apply all noises to input.env.factors
  for (n in noise){
    
    if (n[["type"]]=="add_factor"){
      env.factors <- c(env.factors, n[["target"]])
    }
    
    if (n[["type"]]=="remove_factor"){
      env.factors <- env.factors[!grepl(n[["target"]], env.factors)]
    }
  }
  
  return(unlist(env.factors))
}
