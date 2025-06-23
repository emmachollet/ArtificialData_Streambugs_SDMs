library(caret)

standard.deviance <- function(data, lev=c("Present","Absent"), model=NULL){
    
    
    # This function computes the standard deviance metric for given data, factors
    # and model. This function is used as loss when training the caret models 
    #
    # arguments:
    #   - data: the data to use to compute the standard deviance.
    #   - lev: the factors of the classification of the taxa
    #   - model: the model to use to compute the standard deviance
    #
    # returns:
    #   - standard.dev: the loss to return
    
    
    no.obs <- dim(data)[1]
    likelihood <- 1:no.obs
    
    likelihood[which(data$obs==lev[1])] <- data[which(data$obs==lev[1]), lev[1]]
    likelihood[which(data$obs==lev[2])] <- data[which(data$obs==lev[2]), lev[2]]
    
    threshold <- 1e-4
    likelihood[which(likelihood < threshold)] <- threshold
    
    standard.dev <- -2*sum(log(likelihood)) / no.obs
    names(standard.dev) <- "StandardizedDeviance"
    
    return(standard.dev)
}

lowest <- function (x, metric, maximize = F){
    
    
    # returns the lowest value of a specific metric 
    #
    # arguments:
    #   - x: dataframe containing different metric for each hyperparameters
    #   - metric: metric in datadframe x to get the lowest
    #   - maximize: flag not being used
    #
    # returns:
    #   - best: the lowest value of the metric  
    
    
    best <- which.min(x[, metric])
    return(best)
}


simulate_binomial_data <- function(n = 1000,
                                   intercept = -1,
                                   beta = c(x1 = 0.8, x2 = -0.5, x3 = 0.3),
                                   seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    
    # Simulate predictors
    x1 <- rnorm(n)                        # continuous normal
    x2 <- sample(0:1, n, replace = TRUE) # binary
    x3 <- runif(n, 0, 10)                # uniform
    
    # Compute linear predictor
    eta <- intercept + beta["x1"] * x1 + beta["x2"] * x2 + beta["x3"] * x3
    
    # Inverse logit to get probability
    p <- 1 / (1 + exp(-eta))
    
    # Simulate binary response
    y <- rbinom(n, size = 1, prob = p)
    
    # Return data frame
    data <- data.frame(
        y = factor(y),
        x1 = x1,
        x2 = factor(x2),
        x3 = x3
    )
    
    data$y <- ifelse(data$y == 1, "present", "absent")
    return(data)
}



simulate_binomial_data2 <- function(n = 1000,seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    x <- c(rep(0, n / 2), rep(1, n / 2))
    y <- c(rep(0, n / 2), rep(1, n / 2))  # perfectly separated
    
    data = data.frame(y=y, x1=x)
    
    data$y2 <- ifelse(data$y == 1, "present", "absent")
    
    data
}

# --------------------------------------------------------

data <- simulate_binomial_data()



trainctrl <- trainControl(method='cv',                  
                          number= 3, #nb.folds,                     
                          # index=folds.train,
                          classProbs=T,                 
                          summaryFunction=standard.deviance,
                          selectionFunction=lowest,
                          # verboseIter=TRUE)
                          verboseIter=F)

caret.model <- train(form= y~x1+x2+x3,
                     data=data,
                     method="glm",
                     trControl=trainctrl)
coef <- caret.model$finalModel$coefficients

caret.model2 <- train(form= y~x1+x2+x3,
                     data=data,
                     method="glm",
                     trControl=trainctrl)
coef2 <- caret.model2$finalModel$coefficients




# --------------------------------------------------------

data.x <- simulate_binomial_data2()


mod=glm(y~x1, family = binomial, data=data.x)
mod$converged


trainctrl <- trainControl(method='cv',                  
                          number= 3, #nb.folds,                     
                          # index=folds.train,
                          classProbs=T,                 
                          summaryFunction=standard.deviance,
                          selectionFunction=lowest,
                          # verboseIter=TRUE)
                          verboseIter=F)

caret.model <- train(form= y2~x1,
                     data=data.x,
                     method="glm",
                     trControl=trainctrl)
coef3 <- caret.model$finalModel$coefficients

caret.model$finalModel$converged

