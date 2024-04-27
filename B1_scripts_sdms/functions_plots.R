plot.ice <- function(models.performance,
                     env.factor,
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
  #   - env.factor: The environmental factor to variate on training sample
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
  
  # observations <- models.performance[["null"]][["entire_dataset"]][["training"]][["Simuliidae"]][["observation"]] 
  models = lapply(models.performance,
                  FUN=function(x){
                    lapply(x[["entire_dataset"]][["training"]],
                           FUN=function(y){y["model"]})
                  })
  
  
  # 1 generate df_ext ==========================================================
  
  # get range of env.factor
  min.env.factor <- min(observations[[env.factor]])
  max.env.factor <- max(observations[[env.factor]])
  
  # create range of all env.factor
  range.size <- max.env.factor-min.env.factor
  step.size  <- range.size/(resolution-1)
  
  env.factor.seq <- seq(from=min.env.factor,
                        to=max.env.factor,
                        by=step.size)
  
  # filter observation to keep only nb.sample observation
  sampled.observations <- observations[sample(nrow(observations), nb.sample),]
  
  env.factor.sampled <- sampled.observations[, env.factor]
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
  sampled.observations[env.factor] <- rep(env.factor.seq, times = nb.sample)
  
  # reset observation index
  rownames(sampled.observations) <- 1:nrow(sampled.observations)
  
  # 2 make prediction for every model ==========================================
  
  models.names <- names(models)
  models.names <- models.names[models.names!="null"] # remove null model from list of models 
  #models.names <- models.names[models.names!="gamloess"]
  
  predictions <- lapply(models.names,
                        FUN=make.prediction,
                        models,
                        sampled.observations,
                        taxa,
                        input.env.factors)
  
  names(predictions) <- models.names
  
  
  list.obs.df <- lapply(models.names,
                        function(name, obs, pred){
                          
                          obs["pred"]  <- predictions[[name]]
                          obs["model"] <- name
                          
                          return(obs)
                        },
                        sampled.observations,
                        predictions)
  
  final.observation <- bind_rows(list.obs.df, .id = "column_label")
  
  # unstandardize
  env.fact.std.const <- standardization.constant %>% filter(env_fact==env.factor)
  
  mean      <- as.numeric(env.fact.std.const[,"mean"])
  std.dev   <- as.numeric(env.fact.std.const[,"standard_deviation"])
  
  env.factor.sampled            <- env.factor.sampled*std.dev+mean
  final.observation[env.factor] <- final.observation[env.factor]*std.dev+mean
  
  return(list("observations"=final.observation, 
              "env.factor.sampled"=data.frame(variable=env.factor.sampled)))
  
}
