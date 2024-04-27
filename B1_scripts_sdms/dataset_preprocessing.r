# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Select taxa with intermediate prevalence ----

select.taxa.interm.prev <- function(data, 
                                    prev.taxa, 
                                    prev.threshold = 0.05){
    # This function returns the data with only the taxa having an intermediate prevalence.
    #
    # arguments:
    #   - data: the dataframe where the taxa list should be subselcted by prevalence
    #   - prev.taxa: the dataframe with the list of all taxa and their prevalence
    #   - prev.threshold: a numeric between 0 and 1 used as threshold above/below 
    #                     the taxa prevalence is selected
    #
    # returns:  
    #   - data: the dataframe with only taxa with selected prevalence
    
    taxa.col.index = which(grepl("Occurrence.", colnames(data)))
    data.info <- data[,-taxa.col.index]
    prev.taxa.interm.prev <- prev.taxa[which(prev.taxa$Prevalence > prev.threshold &
                                                prev.taxa$Prevalence < (1 - prev.threshold)), ]
    vect.taxa.interm.prev <- prev.taxa.interm.prev[order(prev.taxa.interm.prev$Prevalence), "Occurrence.taxa"]
    data <- cbind(data.info, data[,vect.taxa.interm.prev])
    
    return(data)
} 

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Reorganizing data ----

categorize.taxa.occurence <- function(data){
  # This function takes a dataframe as input and respectively replace the 
  # numerical 0s and 1s describing taxa's occurrences by categorical "absent"  
  # and "present".
  #
  # arguments:
  #   - data: the dataframe where the taxa's occurrence need to be categorized
  #
  # returns:  
  #   - data: the dataframe with taxa's occurrences categorized
  
  taxa.col.index = which(grepl("Occurrence.",colnames(data)))
  
  for (i in taxa.col.index) {
    data[which(data[,i] == 0),i] <- "absent"
    data[which(data[,i] == 1),i] <- "present"
    data[,i] = as.factor(data[,i])
  }
  return(data)
} 

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Preprocessing data: Splitting and standardizing data  ----

preprocess.data <- function(data,
                            dir,
                            split.type="FIT",
                            nb.split=3,
                            splitting.criterion="ReachID"){
  
  # This function split the data according to a split.type. If split.type is
  # "CV", then the data get split in a 3-Fold for cross validation. If 
  # split.type is "ODG", data gets split for Out-of-Domain Generalization
  #
  # arguments:
  #   - data:       the dataframe to split
  #   - dir:        the directory where to save the data when split
  #   - split.type: the type of split to be perform. It can be either "CV" for 
  #                 a 3-fold cross-validation, "ODG" for Out-of-Domain Generali-
  #                 zation and "FIT" for no split.
  #
  # returns:
  #   - splits:     the splits performed according to split.type and  standard-
  #                 ized
  
  
  # split.type must be either "CV", "ODG", or "FIT". If not, data is not split
  if (!((split.type == "CV")|(split.type == "ODG")|(split.type == "FIT"))){
    split.type <- "FIT"
  }
  
  # get file.name
  split.file <- paste0(dir, "preprocessed_dataset_", split.type, ".rds")
  const.file <-paste0(dir, "standardization_constant_", split.type, ".rds")
  
  # Add other temperatures to env fact list
  all.temp <- colnames(data)[which(grepl("temperature.lm", colnames(data)))]
  list.env.fact <- c(ENV.FACT.FULL.COLNAMES, all.temp)
  
  # Either read splits from file or create them
  if (file.exists(split.file)){  
    
    cat("> File with data splits already exists\n")
    cat("> Read splits from file: \n", split.file,"\n")
    splits <- readRDS(file=split.file)
 
  } else {
    
    cat("> No data splits exist yet, they will be created.\n")
    
    if (split.type=="CV"){
      splits.const <- cv.preprocess.data(data,
                                         list.env.fact,
                                         nb.split,
                                         splitting.criterion)
    } else if (split.type=="ODG"){
      splits.const <- odg.preprocess.data(data, list.env.fact)
    } else { # FIT
      splits.const <- fit.preprocess.data(data, list.env.fact, splitting.criterion)
    }
    
    splits <- splits.const[["splits"]]
    const  <- splits.const[["const"]]
    folds.train <- splits.const[["folds.train"]]
    
    saveRDS(splits, file=split.file)
    saveRDS(const, file=const.file)
    cat("> Splits have been created and saved in:\n   ", split.file, ".\n")
    cat("> Standardization constants have been saved in:\n   ", const.file, ".\n")
    
  }
  
  return(splits)
}

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Subfunctions used by the above functions ----


cv.preprocess.data <- function(data, env.fact, nb.split, splitting.criterion){
  # This function splits the data set for cross-validation with 3 folds. Moreover
  # it directly standardizes the data (without data-leakage). This function is a 
  # subfunction of split.data
  
  # Creating folds
  # No magic number
  
  #folds <- groupKFold(data$SiteId, nb.split) # Keep same sites in same split to avoid data leakage
  folds <- groupKFold(data[,splitting.criterion], nb.split)
  
  splits <- list()
  const  <- list()
  splits.names <- character()
  
  # standardizing the folds
  for (fold.index in seq(nb.split)){
    
    fold.name <- paste("Fold", fold.index, sep="")
    
    splits.names <- c(splits.names, paste("Split", fold.index, sep=""))
    
    train.test <- standardize.data(env.fact,
                                   train=data[folds[[fold.name]],],
                                   test=data[-folds[[fold.name]],])
    
    train     <- train.test[["train"]]
    test      <- train.test[["test"]]
    std.const <- train.test[["const"]]
    folds.train <- groupKFold(train[,splitting.criterion], 3)
    
    splits <- c(splits,
                list(list("Training data" = train, 
                          "Testing data"  = test,
                          "Folds train"   = folds.train)))
    const <- c(const, list(std.const))
    
  } 
  
  names(splits) <- splits.names
  names(const) <- splits.names
  
  return(list("splits"=splits, "const"=const))
}

odg.preprocess.data <- function(data, env.fact){
  # This function split the data set for Out-of-Domain generalization. Then it 
  # standardize the data without data leakage. This function is a subfunction of 
  # preprocess.data
  
  # TODO: write data split function for Out-of-Domain generalization. Currently
  # does nothing, i.e. returns entire dataset without standardization  
  splits <- list("Entire dataset"=data)
  return(list("splits"=splits, "const"=NULL))
}

fit.preprocess.data <- function(data, env.fact, splitting.criterion = "ReachID"){
  # This function standardizes the data when there is no split happening, i.e.
  # split.type = "FIT". This function is a subfunction of split.data
  
  # Standardize
  train.test <- standardize.data(env.fact,
                                train=data,
                                test=NULL)
  data.standardized <- train.test[["train"]]
  std.const         <- train.test[["const"]]
  folds       <- groupKFold(data[,splitting.criterion], 3)
  # TODO: save standardization const to file
  
  splits <- list("Entire dataset" = data.standardized,
                 "Folds train"    = folds)
  const  <- list("Entire dataset" = std.const)

  return(list("splits"=splits, "const"=const))
}

standardize.data <- function(env.fact, train, test=NULL){
  
  # dataframe to save standardization constants (empty with three columns for
  # the environmental factors and their mean and standard deviation)
  std.const.df = data.frame(matrix(nrow=0, ncol=3))
  colnames(std.const.df) = c("env_fact", "mean", "standard_deviation")
  
  for (fact in env.fact){
    
    mean <- mean(train[[fact]])
    standard.deviation <- sd(train[[fact]]) 
    
    # save standardization constant to dataframe
    std.const.df[nrow(std.const.df) + 1,] = c(fact,mean, standard.deviation)
    
    # standardizing training set
    train[[fact]] <- standardize.function(train[[fact]], mean, standard.deviation)
    
    # standardizing testing set if it exists
    if (!(is.null(test))){
      test[[fact]] <- standardize.function(test[[fact]], mean, standard.deviation)
    }
  }
  
  
  return (list("train"=train, "test"=test, "const"=std.const.df))
}


standardize.function <- function(list, mean, standard.deviation){
  return <- unlist(lapply(list, function(x){return((x-mean)/standard.deviation)}))
}

