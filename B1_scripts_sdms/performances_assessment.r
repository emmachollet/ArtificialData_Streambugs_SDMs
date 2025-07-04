summarize.all.results <- function(models.data, prev.taxa, taxa.colnames){
  
  # This function takes as input the models data generated by the function 
  # apply.models and restructure it to create the data needed for plotting 
  # Fig. 1. The data structure of models data is as follows:
  #     -> model name
  #       -> split
  #         -> training, testing
  #           -> taxa 
  #             -> model, observation, prediction factors, prediction 
  #                probabilities, likelihood, performance
  # The data is extracted to form a dataframe with the following column: taxa,
  # prevalence, taxonomic_level, model, dev, auc, fit_pred. Meaning each column 
  # describe for a taxa, a model and either training (fit) or testing (pred)
  # the obtained scores (i.e. dev and auc). The scores are averaged over the 
  # different splits.
  #
  # arguments:
  #   - models.data: Where the scores of the different algorithm is obtained 
  #   - prev.taxa: Where the information about the taxa is obtained (i.e taxa
  #                name, taxa prevalance, taxa family level)
  #
  # returns:
  #   - all.results: A dataframe whose columns are taxa, prevalence, 
  #                  taxonomic_level, model, dev, auc, fit_pred.
  
  # to debug
    # models.data <- models.cv
    # models.data <- models_CV_15
    
  # get dataframe 
  all.results <- data.frame(prev.taxa) 
  
  # # remove useless columns
  # all.results <- all.results %>% select(-Missing.values) # remove useless col
  # 
  # remove useless rows
  all.results <- all.results %>% filter(all.results$Occurrence.taxa %in% taxa.colnames)
  
  # keep only relevant columns and give them appropriate names
  all.results <- all.results[c("Taxon", "Prevalence", "Taxonomic.level")]
  colnames(all.results) <- c('taxa','prevalence','taxonomic_level')
  
  # Change Occurence.taxaname to only taxaname
  all.results['taxa'] <- apply(all.results['taxa'],
                               MARGIN=1,
                               FUN=function(x){gsub('Occurrence.', '', x)})
  
  
  # then for each model : dev.fit, dev.pred, auc.fit, auc.pred, likelihood.ratio, auc.ratio,
  
  # models.data structure:model -> split -> training, testing -> taxa -> list performance
  
  model.names <- names(models.data)
  for (model in model.names){
      # model <- model.names[2]
    performance <- models.data[[model]]
    
    nb.split <- length(performance)
    
    dev.fit.sum <- list(double(nrow(all.results)))
    auc.fit.sum <- list(double(nrow(all.results)))
    dev.pred.sum <- list(double(nrow(all.results)))
    auc.pred.sum <- list(double(nrow(all.results)))
    
    
    for (split.index in seq(nb.split)){ # ! nb.split not declared as function input
      # split.index <- 2
      training.perf <- performance[[split.index]][["training"]]
      
      dev.fit <- lapply(taxa.colnames,
                        FUN=get.metric,
                        training.perf,
                        "standard_deviance")
      
      auc.fit <- lapply(taxa.colnames,
                        FUN=get.metric,
                        training.perf,
                        "auc")
      
      dev.fit.sum <- mapply(sum, dev.fit.sum, dev.fit, SIMPLIFY=FALSE)
      auc.fit.sum <- mapply(sum, auc.fit.sum, auc.fit, SIMPLIFY=FALSE)
      
       
      if (length(performance[[split.index]])>1){
        testing.perf <- performance[[split.index]][["testing"]]
        
        dev.pred <- lapply(taxa.colnames,
                           FUN=get.metric,
                           testing.perf,
                           "standard_deviance")
        auc.pred <- lapply(taxa.colnames,
                           FUN=get.metric,
                           testing.perf,
                           "auc")
        
        dev.pred.sum <- mapply(sum, dev.pred.sum, dev.pred, SIMPLIFY=FALSE)
        auc.pred.sum <- mapply(sum, auc.pred.sum, auc.pred, SIMPLIFY=FALSE)
      } 
    }
    
    dev.fit.avg <- lapply(dev.fit.sum, FUN=function(x){x/nb.split})
    auc.fit.avg <- lapply(auc.fit.sum, FUN=function(x){x/nb.split})
    
    all.results[paste0(model, "_dev_fit")] <- unlist(dev.fit.avg)
    all.results[paste0(model, "_auc_fit")] <- unlist(auc.fit.avg)
    
    if (length(performance[[split.index]])>1){
   
      dev.pred.avg <- lapply(dev.pred.sum, FUN=function(x){x/nb.split})
      auc.pred.avg <- lapply(auc.pred.sum, FUN=function(x){x/nb.split})
      
      dev.diff <- mapply('-', dev.fit.avg, dev.pred.avg, SIMPLIFY=FALSE)
      auc.diff <- mapply('-', auc.fit.avg, auc.pred.avg, SIMPLIFY=FALSE)  
      
      likelihood.ratio <- lapply(dev.diff, FUN=function(x){exp(-x/2)})
      auc.ratio        <- lapply(auc.diff, FUN=function(x){exp(-x/2)})
      
      all.results[paste0(model, "_dev_pred")] <- unlist(dev.pred.avg)
      all.results[paste0(model, "_auc_pred")] <- unlist(auc.pred.avg)
      
      all.results[paste0(model, "_likelihood_ratio")] <- unlist(likelihood.ratio)
      all.results[paste0(model, "_auc_ratio")]        <- unlist(auc.ratio)
    }
  }
  
  
  
  return(all.results)
}

get.metric <- function(taxa, performance, metric.name){
  
  # This function is a helper function for the summarize all results function.
  # In a list of performance of a model, it fetches the metric associated to
  # its name and a taxa.
  
  # arguments:
  #   - taxa: the name of the taxa for which the score needs to be extracted 
  #   - performance: the data structure storing all the score metric of a model
  #                  for different taxa
  #   - metric.name: the name of the metric for which we want the value
  #
  # returns:
  #   - metric: the score associated to a taxa and with name metric.name 
  
  taxa.short <- gsub("Occurrence.", "", taxa)
  metric <- performance[[taxa.short]][[metric.name]]
  
  return(metric)
}

restructure.all.results <- function(all.results, models){
  
  # This function takes as input a dataframe where there is one different column
  # for every model, for every metric and for fitting and predicting, e.g. there
  # are column ann.auc.fit or rf.dev.pred. This method reorder the those columns
  # such that (#model*#metric*#fit.pred = 5*2*2= 50) columns get reorderd in  
  # only three columns model model, metric, fit.pred, 
  #
  # arguments:
  #   - all.results: the dataframe to restructure
  #
  # returns:
  #   - fit.pred.all.results: the restructured dataframe 
  
  name.models <- names(models)
    
  model.all.results <- data.frame()
  for (name.model in name.models) {
      # name.model <- name.models[1]
      model.all.results <- rbind(model.all.results,
          extract.model.results(name.model, all.results))
  }
  # model.all.results <- rbind(extract.model.results("Null", all.results),
  #                            extract.model.results("GLM", all.results),
  #                            extract.model.results("GAM", all.results),
  #                            extract.model.results("RF", all.results))#,
  #                           # extract.model.results("ann", all.results))
  
  model.all.results <- model.all.results[c("taxa",
                                           "prevalence",
                                           "taxonomic_level",
                                           "model",
                                           "dev_fit",
                                           "auc_fit",
                                           "dev_pred",
                                           "auc_pred",
                                           "likelihood_ratio",
                                           "auc_ratio")]
  
  fit.pred.all.results <- rbind(extract.fit.pred.results("fit", model.all.results),
                                extract.fit.pred.results("pred", model.all.results))
  
  return(fit.pred.all.results)
}

extract.model.results <- function(model.name, all.results){
  
  # This function is a helper function for the restructure results function. It
  # iterates over the all.results data.frame and only retains the metric columns 
  # which corresponds to the given model.name. Moreover it adds one column 
  # "model" filled with the given model name
  #
  # arguments:
  #   - model.name: The name model of the model whose metrics are to be kept.  
  #   - all.results: dataframe storing the metrics for the different models
  #
  # returns:
  #   - model.results: dataframe storing the metrics only for the model with 
  #                    name given as argument.
  
  patterns <- c("taxa", "prevalence", "taxonomic_level",
                  paste0("^", model.name))
    
  model.results <- data.frame(all.results %>%
                                  select(matches(paste(patterns, collapse = "|"))))
 
  model.results["model"] <- model.name
  
  new.colnames <- c("taxa",
                    "prevalence",
                    "taxonomic_level", 
                    "dev_fit",
                    "auc_fit",
                    "dev_pred",
                    "auc_pred",
                    "likelihood_ratio",
                    "auc_ratio",
                    "model")
  
  # if the desired model name is not find in the columns of all.results, an 
  # empty dataframe is returned
  if (length(model.results)!=length(new.colnames)){
    
    empty.df = data.frame(matrix(nrow=0,
                                 ncol=length(new.colnames))) 
    colnames(empty.df) = new.colnames
    
    return(empty.df)
  }
  
  
  colnames(model.results) <- new.colnames
  
  return(model.results)
}

extract.fit.pred.results <- function(fit.pred.name, all.results){
  
  # This function is a helper function for the restructure results function. It
  # iterates over the all.results data.frame and only retains the metric columns 
  # which corresponds to the either "fitting" or "predicting". Moreover it adds one 
  #  column "fit.pred" filled with either "fitting" or "predicting".
  #
  # arguments:
  #   - fit.pred.name: Whether the metrics to be kept should be the one of the 
  #                    of the fitting or the predicting.
  #   - all.results: dataframe storing the metrics for either fitting or 
  #                  predicting. 
  #
  # returns:
  #   - model.results: dataframe storing the metrics only for either fitting or
  #                    predicting, depending on "fit.pred.name" argument.
  
  fit.pred.results <- data.frame(all.results %>% select(contains(c("taxa",
                                                                   "prevalence",
                                                                   "taxonomic_level",
                                                                   "model",
                                                                   fit.pred.name))))
  
  
  fit.pred.results["fit_pred"] <- fit.pred.name
  
  colnames(fit.pred.results) <- c("taxa",
                                  "prevalence",
                                  "taxonomic_level",
                                  "model",
                                  "dev",
                                  "auc",
                                  "fit_pred")
  
  return(fit.pred.results)
}

plot.ice <- function(models.performance,
                     select.env.fact,
                     taxa,
                     standardization.constant,
                     observations,
                     nb.sample=100,
                     resolution=200,
                     input.env.factors=c("Temperature"                      = "temperature",
                                         "Flow velocity"                    = "velocity",          # FV
                                         "Riparian agriculture"             = "A10m",              # A10m
                                         "Livestock unit density"           = "cow.density",       # LUD
                                         "Insecticide application rate"     = "IAR",               # IAR
                                         "Urban area"                       = "urban.area",        # Urban
                                         "Forest-river intersection"        = "FRI",               # FRI
                                         "Forest-river intersection buffer" = "bFRI",              # bFRI
                                         "Width variability"                = "width.variability")){
  
  
  # This function is used to create the dataframe to plot Individual Conditional
  # Expectation (ICE). To do so, first the observation used to train the models 
  # are extracted from the argument "model.performance". Then from these
  # observations, nb.sample of them are sampled. Then each of this observation 
  # are duplicated the number of time stored in argument resolution. For the 
  # env.factor given as argument, the value of each each different sample (the
  # "resolution" time) are modified such that they go from the minimum value to^
  # the maximal value that the env.factor was taking in the observation. Once
  # this new dataframe of observation has been created, prediction probability
  # are performed for every model. The final dataframe with new prediction and 
  # corresponding dataframe is returned, alongside, the env.factor value of the 
  # the original sampled observation.
  #
  # arguments:
  #   - model.performance: data structure containing the trained models as well
  #                        all the observations used to train the models.
  #   - select.env.fact: The environmental factor to variate on training sample
  #   - taxa: the taxa to variate values for plotting the ICE
  #   - standardization.constant: used to scale back the standardized
  #                               environmental factors
  #   - nb.sample: the number of observation to sample
  #   - resolution: the number of duplicate of every sample
  #
  # returns:
  #   - observations: the observations sampled with their duplicates, as well as
  #                   prediction for every model. 
  #   - env.factor.sampled: the original values of the sampled observation
  
    # # to debug
    # models.performance=models.fit
    # # taxa=taxon.under.obs
    # taxa=names.taxa[c(1:2)]
    # standardization.constant=std.const.ice[[1]]
    # observations = preprocessed.data.ice$`Entire dataset`
    # nb.sample=no.sites
    # resolution=no.steps
    
  models.names <- names(models.performance)
  models.sdm <- models.performance[which(models.names!="Null" & models.names!="null")] # remove null model from list of models
  
  models.trained = lapply(models.sdm,
                          FUN=function(x){
                            lapply(x[["entire_dataset"]][["training"]],
                                   FUN=function(y){y["model"]})
                          })
  
  # 1 generate df_ext ==========================================================
  
  # get range of env.factor
  min.env.factor <- min(observations[, select.env.fact])
  max.env.factor <- max(observations[, select.env.fact])
  
  # create range of all env.factor
  range.size <- max.env.factor-min.env.factor
  step.size  <- range.size/(resolution-1)
  
  env.factor.seq <- seq(from=min.env.factor,
                        to=max.env.factor,
                        by=step.size)
  
  # filter observation to keep only nb.sample observation
  sampled.observations <- observations[sample(nrow(observations), nb.sample),]
  
  env.factor.sampled <- sampled.observations[, select.env.fact]
  #rownames(sampled.observation) <- 1:nb.sample
  
  # add column with observation number
  sampled.observations['observation_number'] <- 1:nb.sample
  
  # repeat each observation "resolution" times
  sampled.observations <- data.frame(lapply(sampled.observations,
                                            rep,
                                            resolution))
  
  # order rows by observation
  sampled.observations <- sampled.observations[order(sampled.observations$observation),]
  
  # for the factor under consideration, change the repeated value by the range
  sampled.observations[select.env.fact] <- rep(env.factor.seq, times = nb.sample)
  
  # reset observation index
  rownames(sampled.observations) <- 1:nrow(sampled.observations)
  
  # 2 make prediction for every model ==========================================
  
  models.names <- names(models.trained)
  
  predictions <- lapply(models.names,
                        FUN=make.prediction,
                        models.trained,
                        sampled.observations,
                        taxa,
                        input.env.factors)
  
  names(predictions) <- models.names
  
  
  list.obs.df <- lapply(models.names,
                        # name <- models.names[1]
                        function(name, obs, pred){
                          
                            # obs <- sampled.observations
                            # pred <- predictions
                          obs  <- pred[[name]][,c("ReachID", input.env.factors, taxa)]
                          obs["model"] <- name
                          
                          return(obs)
                        },
                        sampled.observations,
                        predictions)
  
  final.observation <- bind_rows(list.obs.df, .id = "column_label")
  
  # unstandardize
  env.fact.std.const <- standardization.constant %>% filter(env_fact==select.env.fact)
  
  mean      <- as.numeric(env.fact.std.const[,"mean"])
  std.dev   <- as.numeric(env.fact.std.const[,"standard_deviation"])
  
  env.factor.sampled            <- env.factor.sampled*std.dev+mean
  final.observation[select.env.fact] <- final.observation[select.env.fact]*std.dev+mean
  
  return(list("observations"=final.observation, 
              "env.factor.sampled"=data.frame(variable=env.factor.sampled)))
  
}

make.prediction <- function(model.name, models, obs, taxa, input.env.factors){
  
  # model.name <- models.names[1]
  # obs <- sampled.observations
  # models <- models.trained
  # models.orig <- models
    
  # This is a helper function used by the plot.ice function. It uses a model
  # passed as argument, to make the prediction corresponding to the observations
  # passed as argument for a given taxa.
  #
  # arguments:
  #   - model.name: name of the model to use
  #   - models: list containing all the models.
  #   - obs: the observation to use to make prediction 
  #   - taxa: which taxa to predicts 
  #
  # returns:
  #   - single.taxa.prediction: the prediction for a the given taxa performed 
  #                             with the given model.
  
  
  taxa.short <- gsub("Occurrence.", "", taxa)
  df.pred <- obs
  
  for(taxon in taxa.short){
      
      # taxon <- taxa.short[1]
      model <- models[[model.name]][[taxon]]
      # input.env.factors <- model$model$coefnames
      
      if(model.name == "ANN"){
          
          obs.matrix  <- as.matrix(obs[ ,input.env.factors])
          predictions <- model[[1]] %>% predict(obs.matrix)
          
          colnames(predictions) <- taxa.colnames # ! Global variable !!
          
          single.taxa.prediction <- predictions[,paste0("Occurrence.", taxon)]
          
      } else {
          
          single.taxa.prediction.df <- predict(model, obs, type='prob')
          single.taxa.prediction    <- single.taxa.prediction.df[[1]][, "Present"]
          
      }
      
      df.pred[, taxon] <- single.taxa.prediction
  }
  
  
  return(df.pred)
}



