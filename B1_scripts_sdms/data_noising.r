
add.gaussian.noise <- function(data, col.name, standard.deviation, min_value=-Inf, max_value=Inf){
  
  # This function adds gaussian noise to a specific column of a dataframe
  #
  # arguments:
  #   - data:                   the data.frame to which noise should be added
  #   - col.name:               to which column of the dataframe the noise 
  #                             should be added
  #   - standard.deviation:     the standard deviation of the gaussian noise to 
  #                             add
  #   - min_value:              minimum value that can take a value of the column
  #                             col.name
  #   - max_value:              maximum value that can take a value of the column
  #                             col.name
  #
  # returns:
  #   - the dataframe with noise added
  
    # col.name <- "tempmaxC"
    # rnorm(nrow(data), mean=0, sd=1.5)
    
  # Adding noise
  data[,col.name] <- data[,col.name] + rnorm(nrow(data), mean=0, sd=standard.deviation)
  # Clipping
  data[,col.name] <- clipper(data[,col.name], min_value, max_value)
  
  return(data)
  
}


binomial.flipper <- function(data, col.name, chance.flipping){
  
  # This function flips (0->1, 1->0) the presence of taxa of a specific column 
  # of a dataframe.
  #
  # arguments:
  #   - data:                   the dataframe where the missclassifcation 
  #                             should be happen.
  #   - col.name:               which taxa column of the data frame should be 
  #                             missclassified.
  #   - chance.missclassifying: the likelihood of flipping for each observation 
  #                             of the data set.
  #
  # returns:
  #   - the dataframe with missclassified taxa
  
  data[,col.name] <- (data[,col.name] + rbinom(n=nrow(data), size=1, prob=chance.flipping)) %% 2
  
  return(data)
  
}

missdetection <- function(data, col.name, chance.missdetection){
  
  # This function missdetect entries. That means that if a for an observation
  # the taxa of col.name has been observed as present, it is changed to observed
  # as absent, with a probability chance.missdetection.
  #
  # arguments:
  #   - data:                   the dataframe where the missdetection should 
  #                             take place.
  #   - col.name:               In which taxa column of the dataframe the miss-
  #                             detection should be performed.
  #   - chance.missdetection:   the likelihood of not detecting a taxa that was 
  #                             there. 
  #
  # returns:
  #   - the dataframe with missclassified taxa
  
  data[,col.name] <- ifelse(data[,col.name] == 0, 0, rbinom(dim(data)[1], 1, 1-chance.missdetection))
  return(data)
}

missclassification <- function(data){
  
  # TODO: write a function that missclassify species based on a meaningful set of
  #rule
  
  return(data)
  
}


add.noise <- function(data, noise, env.fact, env.fact.full){
  
    # data            = data.input.orig
    # noise           = list.noise
    # env.fact        = env.factor
    # env.fact.full   = env.factor.full
  noised.data <- data
  
  for (n in noise){
    # n <- noise[[1]]
    if (n[["type"]]=="gaussian"){
      noised.data <- noised.data %>% 
        add.gaussian.noise(col.name=n[["target"]],
                           standard.deviation=n[["amount"]],
                           min_value=n[["parameters"]][["min"]],
                           max_value=n[["parameters"]][["max"]])
    }
    if (n[["type"]]=="missdetection"){
      noised.data <- noised.data %>%
        missdetection(col.name=n[["target"]],
                      chance.missdetection=n[["amount"]])
    }
    if (n[["type"]]=="missclassification"){
      noised.data <- noised.data %>%
        missmissclassification(data)
    }
    if (n[["type"]]=="add_factor"){
    
      new.factor.name                <- n[["target"]]
      noised.data[[new.factor.name]] <- runif(nrow(data), -1, 1)
      env.fact                       <- c(env.fact,
                                          new.factor.name = new.factor.name)
      env.fact.full                  <- c(env.fact.full,
                                          new.factor.name = new.factor.name)
    }
    if (n[["type"]]=="remove_factor"){
      target_factor <- n[["target"]]
      
      env.fact <- env.fact[!grepl(target_factor, env.fact)]
      env.fact.full <- env.fact.full[!grepl(target_factor, env.fact.full)]
      cind.fact.remove <- which(grepl(target_factor, colnames(data)))
      noised.data <- noised.data[,-cind.fact.remove]
      
    }
  }

  return(list("noised data"   = noised.data,
              "env fact"      = env.fact,
              "env fact full" = env.fact.full))
}
  