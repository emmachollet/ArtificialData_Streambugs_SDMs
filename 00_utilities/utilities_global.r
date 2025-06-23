## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## --- utilities needed to preprocess data and run Streambugs  ---
## --- to create synthetic data ---
##
## --- December 2024 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GENERAL ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ~~~~~~~~~~~~~~~~~~~~~~~~~~
## Utilities ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~

list.depth <- function(x) {
    if (!is.list(x)) {
        return(0)
    } else if (length(x) == 0) {
        return(1)
    } else {
        return(1 + max(sapply(x, list.depth)))
    }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plot related functions ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~

print.pdf.plots <- function(list.plots, width = 12, height = width*3/4, 
                            dir.output, info.file.name = "", file.name = "list.plots.pdf", 
                            png = FALSE, png.square = F, png.vertical = F, png.ratio = 1){
    
    # This function process files (gdb) to "sf" objects to plot maps
    #
    # arguments:
    #   - list plots:             list of (gg)plots
    #   - dir.output:             directory to output the pdf
    #   - info.file.name, file.name: information and name of the file
    #   - png:                    if TRUE print PNG (options for square, vertical and modify ratio for quality)
    #
    # returns:
    #   - list.swiss.map.inputs:  list with "sf" inputs for plotting Swiss map 
    #                             with major lakes and rivers 
    
    pdf(paste0(dir.output, info.file.name, file.name, ".pdf"), paper = 'special', width = width, height = height, onefile = TRUE)
    cat("\nPDF:")
    for (n in 1:length(list.plots)) {
        plot(list.plots[[n]])
        cat(" ",n, " ")
    }
    dev.off()
    
    if(png){
        cat("\nPNG:")
        # Print a png file
        if(png.square){
            w = 1280
            h = 960
        } else if(png.vertical){
            w = 960
            h = 1280
        } else {
            w = 1280
            h = 720
        }
        
        for (n in 1:length(list.plots)) {
            # n = 1
            cat(" ",n, " ")
            png(paste0(dir.output, info.file.name, file.name, n, ".png"), width = png.ratio*w, height = png.ratio*h) # size in pixels (common 16:9 --> 1920�1080 or 1280�720)
            plot(list.plots[[n]])
            dev.off()
            
        }    
    }
}

map.inputs <- function(directory){
    
    # This function process files (gdb) to "sf" objects to plot maps
    #
    # library:
    #   - library("sf")           needs to be loaded for specific functions
    # 
    # arguments:
    #   - directory:              the path to the folder containing the gdb files 
    #                             with geographic layer
    #
    # returns:
    #   - list.swiss.map.inputs:  list with "sf" inputs for plotting Swiss map 
    #                             with major lakes and rivers 
    
    list.swiss.map.inputs = list()
    
    # Obtain simple feature for borders of Switzerland, and major lakes and rivers
    list.swiss.map.inputs$ch <- st_read(directory, layer="switzerland", stringsAsFactors = F)
    list.swiss.map.inputs$ch <- filter(list.swiss.map.inputs$ch, NAME=="Schweiz")
    
    list.swiss.map.inputs$rivers.major <- st_read(directory, layer = "major_rivers", stringsAsFactors = F)
    list.swiss.map.inputs$lakes.major <- st_read(directory, layer = "major_lakes", stringsAsFactors = F)
    
    return(list.swiss.map.inputs)
}

maps.env.fact <- function(map.inputs, vect.env.fact, vect.info, data){
    
    # inputs:
    # map.inputs : list with "sf" inputs for plotting Swiss map 
    #              with major lakes and rivers
    # vect.env.fact : named vector of environmental factors (columns of data)
    #                 to be plotted on Swiss map
    # data : dataframe that as X, Y, ReachID, and vect.env.fact as columns
    # data <- na.omit(data)
    
    list.plots <- list()
    
    for(variable in vect.env.fact){
        # variable <- vect.analysis.env.fact[1]
        print(variable)
        name.variable <- names(variable)
        print(name.variable)
        
        p <- ggplot() + 
            geom_sf(data = map.inputs$ch, fill="#E8E8E8", color="black") + 
            geom_sf(data = map.inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE) + 
            geom_sf(data = map.inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE) + 
            geom_point(data = data, aes(x=X, y=Y, color = data[, variable]), size= 5, alpha= 0.8) + 
            # geom_text(aes(x=X, y=Y,  label = data[, "ReachID"]), size = 3, vjust = 0, nudge_y = 0.5) +
            scale_colour_gradient2(name = name.variable,
                                   low = "lightskyblue",
                                   high = "firebrick3") +
            theme_void(base_size = 18) +
            # theme_void(base_size = 18) + 
            theme(panel.grid.major = element_line(colour="transparent"),
                  plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
                  legend.title = element_text(size=24),
                  legend.text = element_text(size=20))
        
        list.plots[variable] <- p
    }
    # plot_maps_variable = function (pair.fact.name, map.inputs, data) {
    #     variable <- pair.fact.name[1]
    #     print(variable)
    #     name.variable <- pair.fact.name[2]
    #     if(is.na(name.variable)){name.variable <- variable}
    #     print(name.variable)
    #     # data[,"ReachID"] <- gsub("CSCF.CH.", "", data[,"ReachID"])
    #     ggplot(data = data[,c("X","Y", "ReachID", variable)]) + 
    #         geom_sf(data = map.inputs$ch, fill="#E8E8E8", color="black") + 
    #         geom_sf(data = map.inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE) + 
    #         geom_sf(data = map.inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE) + 
    #         geom_point(aes(x=X, y=Y, color = data[, variable]), size= 5, alpha= 0.8) + 
    #         # geom_text(aes(x=X, y=Y,  label = data[, "ReachID"]), size = 3, vjust = 0, nudge_y = 0.5) +
    #         scale_colour_gradient2(name = name.variable,
    #                                low = "lightskyblue",
    #                                high = "firebrick3") +
    #         theme_void(base_size = 18) +
    #         # theme_void(base_size = 18) + 
    #         theme(panel.grid.major = element_line(colour="transparent"),
    #               plot.margin = unit(c(0.1,0.1,0.1,0.1), "lines"),
    #               legend.title = element_text(size=24),
    #               legend.text = element_text(size=20))
    # }
    
    # !! almost there, just missing labels 
    # list.plots <- lapply(vect.env.fact, plot_maps_variable, map.inputs  = map.inputs, data = data)
    
    return(list.plots)
}

# small basic ggplot colors function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## ~~~~~~~~~~~~~~~~~~~~~~~
## Taxonomy ####
## ~~~~~~~~~~~~~~~~~~~~~~~

get.df.prev.pres.taxonomy <- function(data.inv, data.taxonomy, catch.variable, vect.catch.select, dir.plots){
    # required input: 
    # data.inv = the wide-format monitoring data of the invertebrates, with taxa in columnnames starting with Occurrence_
    # data.tax = a taxonomic dictionary with all taxa in Switzerland and their taxonomic resolution
    # catch.variable = column name of data.inv of a (geographical/catchment) factor for which we want to get the prevalence
    
    ## add a column with the full name of taxa, and their tax resolution to the data.taxonomy
    df.taxonomy <- data.taxonomy
    df.taxonomy$taxon <- NA
    df.taxonomy$lowest.level <- NA
    for(i in 1:nrow(df.taxonomy)){
        if(df.taxonomy[i,"Species"]!=""){
            df.taxonomy$taxon[i] <- paste(df.taxonomy[i,"Genus"],df.taxonomy[i,"Species"],sep="_" )
            df.taxonomy$lowest.level[i] <- "species"
        } else {
            if(df.taxonomy[i,"Genus"]!=""){
                df.taxonomy$taxon[i] <- paste(df.taxonomy[i,"Genus"])
                df.taxonomy$lowest.level[i] <- "genus"
            } else {
                if(df.taxonomy[i,"Family"]!=""){
                    df.taxonomy$taxon[i] <- paste(df.taxonomy[i,"Family"])
                    df.taxonomy$lowest.level[i] <- "family"
                } else {
                    if(df.taxonomy[i,"Order"]!=""){
                        df.taxonomy$taxon[i] <- paste(df.taxonomy[i,"Order"])
                        df.taxonomy$lowest.level[i] <- "order"
                    } else{
                        if(df.taxonomy[i,"Class"]!=""){
                            df.taxonomy$taxon[i] <- paste(df.taxonomy[i,"Class"])
                            df.taxonomy$lowest.level[i] <- "class"
                        } else {
                            df.taxonomy$taxon[i] <- paste(df.taxonomy[i,"Phylum"])
                            df.taxonomy$lowest.level[i] <- "phylum"
                        }
                    }
                }
            }
        }
    }
    # change species name to match streambugs workflow
    df.taxonomy$taxon.streambugs <- gsub("_", "", df.taxonomy$taxon)
    
    # get taxa vector
    cind.taxa <- which(grepl("Occurrence.", colnames(data.inv)))
    vect.occ.taxa <- colnames(data.inv)[cind.taxa]
    
    # get catchment vector
    # catch.variable <- "Watershed"
    vect.catch.all <- unique(data.inv[, catch.variable])
    
    # prepare dataframe prevalence
    data.prev <- data.frame("Taxon" = vect.occ.taxa)
    data.prev$Taxonomic.level <- NA
    data.prev$Prevalence <- NA
    data.prev[, vect.catch.all] <- NA
    
    # prepare dataframe prevalence
    data.nb.pres <- data.frame("Taxon" = vect.occ.taxa)
    data.nb.pres$Taxonomic.level <- NA
    data.nb.pres$NbPresPoints <- NA
    data.nb.pres[, vect.catch.all] <- NA
    
    for(i in cind.taxa){
        # i <- cind.taxa[6]
        taxon <- colnames(data.inv)[i]
        
        # fill prevalence column
        n.obs <- sum(!is.na(data.inv[,i]))
        n.pres <- sum(na.omit(data.inv[,i]))
        prev <- n.pres/n.obs # calculate the of prevalence of this taxon
        data.nb.pres[which(data.nb.pres[,1] == taxon), "NbPresPoints"] <- n.pres
        data.prev[which(data.prev[,1] == taxon), "Prevalence"] <- round(prev, digits = 3)
        
        for (catch in vect.catch.all) {
            # catch <- vect.catch.all[3]
            temp.df.catch.taxon <- na.omit(data.inv[which(data.inv[,catch.variable] == catch), taxon])
            n.pres.catch <- sum(temp.df.catch.taxon)
            prev.tax.catch <- n.pres.catch/length(temp.df.catch.taxon)
            data.nb.pres[which(data.nb.pres$Taxon == taxon), which(colnames(data.nb.pres) == catch)] <- n.pres.catch
            data.prev[which(data.prev$Taxon == taxon), which(colnames(data.prev) == catch)] <- prev.tax.catch
            # }
        }
        
        # fill taxonomic level column
        rind <- which(df.taxonomy[,"taxon.streambugs"] == sub("Occurrence.","", taxon)) # look for taxon in df.taxonomy
        if (length(rind) != 0) {
            
            data.nb.pres[which(data.nb.pres[,1] == taxon), "Taxonomic.level"] <- df.taxonomy[rind,"lowest.level"]
            data.prev[which(data.prev[,1] == taxon), "Taxonomic.level"] <- df.taxonomy[rind,"lowest.level"]
            
        } else {
            cat(sub("Occurrence.","", taxon), "is not in taxa homogenization \n")
        }
    }
    
    data.prev[,vect.catch.all] <- round(data.prev[,vect.catch.all], 3) # round numbers
    data.prev$Taxon <- gsub("Occurrence.", "", data.prev$Taxon)
    data.nb.pres$Taxon <- gsub("Occurrence.", "", data.nb.pres$Taxon)
    
    list.df.prev.pres <- list("Prevalence" = data.prev, "NbPresPoints" = data.nb.pres)
    list.plots <- list()
    for (n in 1:length(list.df.prev.pres)) {
        # plot distribution prev/nbprespoints across selected catchments
        # n <- 1
        var <- names(list.df.prev.pres)[n]
        plot.data <- list.df.prev.pres[[n]][,c("Taxon", "Taxonomic.level", var, vect.catch.select)]
        plot.data <- gather(plot.data, key = "Catchment", value = Value, -c("Taxon", "Taxonomic.level"))
        plot.data <- filter(plot.data, Taxonomic.level %in% c("species", "genus", "family"))
        
        p <- ggplot(plot.data, aes(x = Catchment, y = Value, fill = Taxonomic.level))
        p <- p + geom_boxplot()
        p <- p + theme_bw()
        p <- p + labs(y=var)
        p
        
        list.plots[[n]] <- p
    }
    
    # add close up on smaller nb of presence points
    q <- p + ylim(0,50)
    list.plots[[length(list.plots) + 1]] <- q
    
    file.name <- "Distribution_Prev_NbPresPoints_SelectCatch"
    print.pdf.plots(list.plots = list.plots, width = 12, height = 8, dir.output = dir.plots, info.file.name = "", file.name = file.name, 
                    png = F)
    
    
    return(list.df.prev.pres)
    
}




plot.histogram.prev.catch <- function(df.taxa.catch, vect.taxa){
    
    list.plots <- list()
    plot.data <- gather(df.taxa.catch, key = Catchment, value = Prevalence, -Taxa)
    
    for (taxon in vect.taxa) {
        # taxon <- vect.taxa[5]
        name.taxon <- gsub("Occurrence.", "", taxon)
        temp.plot.data <- filter(plot.data, Taxa == taxon)
        p <- ggplot(data = temp.plot.data, aes(x=Catchment, y = Prevalence)) 
        p <- p + geom_bar(stat = 'identity', fill = 'lightblue', color = 'black', width = 0.75)
        # p <- p + geom_hline(yintercept = 0.5)
        p <- p + coord_cartesian(ylim = c(0,1))
        p <- p + theme_bw(base_size = 15)
        p <- p + labs(title = paste("Distribution of prevalence for:", name.taxon))
        # p
        list.plots[[taxon]] <- p
    }
    
    return(list.plots)
}


## ~~~~~~~~~~~~~~~~~~~~~~~
## Preprocessing data ####
## ~~~~~~~~~~~~~~~~~~~~~~~

preprocess.data <- function(data,
                            env.fact,
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
  #                 zation, "FIT" for no split or "ICE" for no split but input
  #                 data for ICE plots.
  #
  # returns:
  #   - splits:     the splits performed according to split.type and  standard-
  #                 ized
  
  
  # split.type must be either "CV", "ODG", "FIT" or "ICE". If not, data is not split
  if (!((split.type == "CV")|(split.type == "ODG")|
        (split.type == "FIT")|(split.type == "ICE"))){
    split.type <- "FIT"
  }
  
  # get file.name
  split.file <- paste0(dir, "preprocessed_dataset_", split.type, ".rds")
  const.file <- paste0(dir, "standardization_constant_", split.type, ".rds")
  
  # Either read splits from file or create them
  if (file.exists(split.file)){  
    
    cat("> File with data splits already exists\n")
    cat("> Read splits from file: \n", split.file,"\n")
    splits <- readRDS(file=split.file)
    
  } else {
    
    cat("> No data splits exist yet, they will be created.\n")
    
    if (split.type=="CV"){
      splits.const <- cv.preprocess.data(data,
                                         env.fact,
                                         nb.split,
                                         splitting.criterion)
    } else if (split.type=="ODG"){
      splits.const <- odg.preprocess.data(data, env.fact)
    } else { # FIT or ICE
      splits.const <- fit.preprocess.data(data, env.fact, splitting.criterion)
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


clipper <- function(vec, min_value=-Inf, max_value=Inf){
  
  # This function clips the values of received vector vec to lie in [min_value, 
  # max_value].
  #
  # arguments:
  #   - vec:                    vector of values to clip.
  #   - min_value:              lower bound for clipping vec.
  #   - max_value:              upper bound for clipping vec
  #
  # returns:
  #   - the clipped version of vec
  
  return(pmax(rep(min_value, length(vec)), pmin(vec, rep(max_value, length(vec)))))
}

# Done with ChatGPT May 2025
generate.scenario.names <- function(list.scenarios, na.to.absence, no.taxa, no.models, vect.seeds) {
    # Helper to get either full range (if flag is TRUE) or just the best value
    get_values <- function(scenario) {
        if (scenario$flag) {
            return(as.list(scenario$range))
        } else {
            return(list(best = scenario$range["best"]))
        }
    }
    
    # Get list of value sets for each parameter
    value_sets <- lapply(list.scenarios, get_values)
    
    # Add seeds as a scenario dimension
    value_sets$seed <- as.list(vect.seeds)
    
    # Create all combinations
    scenario_grid <- expand.grid(value_sets, stringsAsFactors = FALSE)
    
    # Generate scenario names
    scenario_names <- apply(scenario_grid, 1, function(row) {
        parts <- c()
        if ("dataset.size" %in% names(row)) {
            parts <- c(parts, paste0(row[["dataset.size"]], "sites"))
        }
        if ("seed" %in% names(row)) {
            parts <- c(parts, paste0(row[["seed"]], "seed"))
        }
        if ("nb.predictors" %in% names(row)) {
            parts <- c(parts, paste0(row[["nb.predictors"]], "pred"))
        }
        if ("noise.temperature" %in% names(row)) {
            parts <- c(parts, paste0(row[["noise.temperature"]], "noisetemp"))
        }
        if ("misdetection" %in% names(row)) {
            parts <- c(parts, paste0(row[["misdetection"]], "misdet"))
        }
        
        # Final name with NA handling, taxa, and model count
        paste0(
            paste(parts, collapse = "_"), "_",
            ifelse(na.to.absence, "NAtoabs", "NAtoNA"), "_",
            no.taxa, "taxa_", no.models, "models_"
        )
    })
    
    return(scenario_names)
}

# generate.scenario.names <- function(list.scenarios, na.to.absence, no.taxa, no.models) {
#     # Helper to get either full range (if flag is TRUE) or just the best value
#     get_values <- function(scenario) {
#         if (scenario$flag) {
#             print(scenario$range)
#             return(as.list(scenario$range))
#         } else {
#             return(list(best = scenario$range["best"]))
#         }
#     }
#     
#     # Get list of value sets for each parameter
#     value_sets <- lapply(list.scenarios, get_values)
#     
#     # Create all combinations
#     scenario_grid <- expand.grid(value_sets, stringsAsFactors = FALSE)
#     
#     # Generate scenario names
#     scenario_names <- apply(scenario_grid, 1, function(row) {
#         parts <- c()
#         if ("dataset.size" %in% names(row)) {
#             parts <- c(parts, paste0(row[["dataset.size"]], "sites"))
#         }
#         if ("nb.predictors" %in% names(row)) {
#             parts <- c(parts, paste0(row[["nb.predictors"]], "pred"))
#         }
#         if ("noise.temperature" %in% names(row)) {
#             parts <- c(parts, paste0(row[["noise.temperature"]], "noisetemp"))
#         }
#         if ("misdetection" %in% names(row)) {
#             parts <- c(parts, paste0(row[["misdetection"]], "misdet"))
#         }
#         
#         
#         # Final name with taxa and model count
#         paste0(paste(parts, collapse = "_"), "_", ifelse(na.to.absence, "NAtoabs", "NAtoNA"), "_", no.taxa, "taxa_", no.models, "models_")
#     })
#     
#     return(scenario_names)
# }

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
## Subfunctions used by the above functions ----


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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# STREAMBUGS ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

exp.transform <-   function(x,intercept=0,curv=0)
{
    #!if curv > 0 and intercept <1: function is curved to the right, if curv < 0 and intercept <1 function is curved to the left
    #!if curv > 0 and intercept >1: function is curved to the left,  if curv < 0 and intercept >1 function is curved to the right
    
    if(curv == 0)
    { 
        y = intercept-(intercept-1)*x
    } else {
        y = intercept - (intercept -1) * (1 - exp(-curv * x)) / (1-exp(-curv)) 
    }
    return(y)
}

write.df.preferences <- function(system.def, Invertebrates, dir.outputs){
    
    # make list with dataframes with two columns: name of env. traits and their classes
    
    # initiate list with feeding types
    feedingtypes <- c("gra", "min", "xyl", "shr", "gat", "aff", "pff", "pre", "par", "other")
    list.preferences <- list(data.frame( "class.names" = paste0("feedingtype_", feedingtypes), "class.values" = rep(1, length(feedingtypes))))
    names(list.preferences) <- "feedingtype"
    
    # add other env. traits to the list
    env.traits <- names(system.def$par.taxaprop.traits)
    short.traits <- c("sapro", "orgmicro", "current", "tempmax")
    for (t in short.traits) {
        # t <- short.traits[1]
        trait <- env.traits[which(grepl(t, env.traits))]
        class.names <- colnames(system.def$par.taxaprop.traits[[trait]][[2]])
        
        par.global <- system.def$par.global.envtraits
        ind <- which(grepl(t, names(par.global)))
        class.values <- par.global[[ind]][["parvals"]]
        if(grepl("temp", t)){ class.values <- class.values - 273.15}
        
        temp.trait.list <- list(data.frame(class.names, class.values))
        names(temp.trait.list) <- trait
        list.preferences <- append(list.preferences, temp.trait.list)
    }
    
    # extract parameters
    parameters <- system.def$par
    
    
    list.df.preferences <- list()
    
    for (n in 1:length(list.preferences)) {
        env.pref <- names(list.preferences)[n]
        cat(env.pref, "classes :", list.preferences[[n]][,1], "\n",
            env.pref, "values :", list.preferences[[n]][,2], "\n")
        classes.pref <- list.preferences[[n]][,1]
        classes.val <- list.preferences[[n]][,2]
        taxa.pref <- parameters[which(grepl(env.pref, names(parameters)))]
        df.pref <- data.frame(Classes = classes.pref, Values = classes.val)
        for (taxon in Invertebrates) {
            df.pref[,taxon] <- NA
            for (c in classes.pref) {
                df.pref[which(df.pref$Classes == c), taxon] <-
                    taxa.pref[which(grepl(taxon, names(taxa.pref)) & grepl(c, names(taxa.pref)))]
            }
        }
        n.taxa <- length(Invertebrates)
        file.name <- paste0(dir.outputs, "df_taxapref_", n.taxa, "Taxa_", env.pref, ".csv")
        cat("Writing:", file.name, "\n")
        write.csv(df.pref, file = file.name)
        
        list.df.preferences[[env.pref]] <- df.pref
    }
    
    return(list.df.preferences)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Probability of occurrence ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


get.prob.obs <- function(ss, M.taxon){
    
    # definitions of parameters from Schuwirth et al., 2016, Func. Ecol. (eqn 2 and 3)
    p.abs   <- 0.1  # y-intercept
    K.abs   <- 100  # predicted abundance at which the probability pabs is 50% of p.abs.
    # parameters above used to compute final pabs:
    # probability that a taxon is absent at sampling date k despite the fact that it has a stable population at sampling site j
    # (e.g. due to short term disturbance by floods or emergence just before the sampling event),
    
    D.drift <- 0.04 # typical number of individuals per m2 that are introduced by drift or misclassification
    p.obs   <- 0.75 # probability of observing (i.e. catching and correctly identifying) a taxon at sampling event k 
    # that is present at site j at average environmental conditions at steady state with one individual per sampling area Aj,k
    
    A.jk <- 0.25^2*8 # area of observation, in MSK methods eight times the area of a foot (squared) (randomly select size 40 EU)
    
    if(ss<0){
        if(ss < -1e-7) warning("res.ss of ",names(ss)," <0: ",ss," set to zero")
        ss <- 0
    }
    
    p.absent <-  p.abs*exp(-log(2)*(ss/M.taxon)/K.abs)
    
    p.ijk.0 <- p.absent + (1-p.absent) * (1-p.obs)^(A.jk*((ss/M.taxon)+D.drift))
    p.ijk.1 <- 1 - p.ijk.0
    
    return(p.ijk.1)
}


calc.prob.occ.obs <- function(res.ss, y.names, M.taxa, obs.error = F){
# calc.prob.occ <- function(res.ss, y.names, M.taxa, observed.abund, p.obs, p.abs, K.abs, D.drift){
        # res.ss <- temp.res.ss
    y.names <- temp.y.names
    if(obs.error == T){
        # definitions of parameters from Schuwirth et al., 2016, Func. Ecol. (eqn 2 and 3)
        p.abs   <- 0.1  # y-intercept
        K.abs   <- 100  # predicted abundance at which the probability pabs is 50% of p.abs.
        # parameters above used to compute final pabs:
        # probability that a taxon is absent at sampling date k despite the fact that it has a stable population at sampling site j
        # (e.g. due to short term disturbance by floods or emergence just before the sampling event),
        
        D.drift <- 0.04 # typical number of individuals per m2 that are introduced by drift or misclassification
        p.obs   <- 0.75 # probability of observing (i.e. catching and correctly identifying) a taxon at sampling event k 
        # that is present at site j at average environmental conditions at steady state with one individual per sampling area Aj,k

        A.jk <- 0.25^2*8 # area of observation, in MSK methods eight times the area of a foot (squared) (randomly select size 40 EU)
        
        p.1 <- NULL

        for(i in 1:length(res.ss)){
            # i <- 5
            ind.y <- which(y.names$y.names==names(res.ss)[i])
            
            group <- y.names$y.groups[ind.y]
            # print(group)
            if(group=="Invertebrates"){
                # site <- y.names$y.reaches[ind.y]
                # hab  <- y.names$y.habitats[ind.y]
                # tax  <- y.names$y.taxa[ind.y]
                
                # ind.j <- which(observed.abund[,"ReachID"] == site & observed.abund[,"Habitat"]==hab)
                
                # dates <- observed.abund[ind.j,"Date"]
                
                M.i <- M.taxa[names(res.ss)[i]] # gDM
                
                # if(length(dates)>0){
                    # for(k in ind.j){
                        
                        # A.jk <- observed.abund[k,"area_m2"]
                        
                        # if(sum(A.jk)==0) warning("area not given, A.jk:",A.jk)
                        
                        if(res.ss[i]<0){
                            if(res.ss[i] < -1e-7) warning("res.ss of ",names(res.ss)[i]," <0: ",res.ss[i]," set to zero")
                            res.ss[i] <- 0
                        }
                        
                        p.absent <-  p.abs*exp(-log(2)*(res.ss[i]/M.i)/K.abs)
                        
                        p.ijk.0 <- p.absent + (1-p.absent) * (1-p.obs)^(A.jk*((res.ss[i]/M.i)+D.drift))
                        p.ijk.1 <- 1 - p.ijk.0
                        
                        # p.0 <- c(p.0, p.ijk.0)
                        # names(p.0)[length(p.0)] <- paste(y.names$y.names[ind.y],observed.abund[k,"Date"],sep="_")
                        
                        p.1 <- c(p.1, p.ijk.1)
                        # names(p.1)[length(p.1)] <- paste(y.names$y.names[ind.y],observed.abund[k,"Date"],sep="_")
                        
                        # if(observed.abund[k,tax]>0) p.ijk <- p.ijk.1 else {
                        #     if(observed.abund[k,tax]==0) p.ijk <- p.ijk.0 else {warning(y.names$y.names,"not observed at",observed.abund[k,"Date"] )}}
                        # 
                        # p.y <- c(p.y,p.ijk)
                        # names(p.y)[length(p.y)] <- paste(y.names$y.names[ind.y],observed.abund[k,"Date"],sep="_")
                        
                        #           cat(log(p.ijk),names(M.i),i,k,"\n")
                        
                        # loglikeli <- loglikeli + log(p.ijk)
                        
                    # }        
                # }
            }
        }

    } else {
        p.abs <- 0 # no probability of being absent if there are some indiv
        K.abs <- 0 # no need to define K.abs if p.abs = 0
        D.drift <- 0 # no indiv introduced by drift
        p.obs   <- 1 # probability of observing (i.e. catching and correctly identifying) a taxon at sampling event k 
        
    }
    
    return(p.1)    
    # return(list( #"loglikeli"=loglikeli,
    #             "p.0"=p.0,"p.1"=p.1))
    #             # ,"p.y"=p.y))
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Additional results ####
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get.plot.data.add.res <- function(catch.results, list.factors.type){
    
    # catch.results <- list.results$RheinabBS
    
    df.res.add <- as.data.frame(catch.results$res$res.add)
    col.df.res.add <- colnames(df.res.add)
    temp.vect.taxa <- catch.results$y.names$taxa[-c(1:2)]
    temp.vect.sites <- catch.results$y.names$reaches
    y.names <- catch.results$y.names$y.names
    streambugs.results <- as.matrix(catch.results$res$res)
    env.data <- catch.results$env.data
    
    # make list of data frames with additional results per taxon per site
    list.df.rates.limfact.taxa <- list()
    for (taxon in temp.vect.taxa) {
        # taxon <- temp.vect.taxa[3]
        c.ind <- which(grepl(taxon, col.df.res.add) & !grepl("__", col.df.res.add))
        temp.df <- df.res.add[,c(1, c.ind)]
        temp.list <- list()
        for (site in temp.vect.sites) {
            # site <- temp.vect.sites[1]
            c.ind2 <- which(grepl(site, colnames(temp.df)))
            temp.df2 <- temp.df[,c(1,c.ind2)]
            c.ind3 <- which(grepl(taxon, colnames(streambugs.results)) & grepl(site, colnames(streambugs.results)))
            temp.df2$Biomass <- streambugs.results[,c.ind3]
            temp.list[[site]] <- temp.df2
        }
        list.df.rates.limfact.taxa[[taxon]] <- temp.list
    }
    
    # make long dataframe for plotting additional result
    plot.data.all.taxa.sites <- data.frame()
    for (taxon in temp.vect.taxa) {
        # taxon <- temp.vect.taxa[3]
        temp.list.df.sites <- list.df.rates.limfact.taxa[[taxon]]
        # list.plot.data.biom.rates.fact[[taxon]] <- list()
        
        for (site in temp.vect.sites) {
            # site <- temp.vect.sites[1]
            temp.df.rates.limfact <- temp.list.df.sites[[site]]
            ind.y.names <- which(grepl(taxon, y.names) & grepl(site, y.names))
            state.variable <- y.names[ind.y.names]
            colnames(temp.df.rates.limfact) <- gsub(paste0(".", state.variable), "", colnames(temp.df.rates.limfact))
            if(grepl("Algae", taxon)){ # if limfactor of algae. add 1 - shade
                temp.df.rates.limfact$fshade <- 1 - env.data[which(env.data$ReachID == site), "shade"]
            }
            
            temp.plot.data <- gather(temp.df.rates.limfact, key = Factor, value = Value, -c("time"))
            
            temp.plot.data$Taxon <- taxon
            temp.plot.data$Site <- site
            
            # add type of limitation factor for easier plotting later
            temp.plot.data$FactorType <- NA
            for (type in fact.types) {
                # type <- fact.types[1]
                r.ind <- which(temp.plot.data$Factor %in% list.factors.type[[type]])
                temp.plot.data[r.ind, "FactorType"] <- type
            }
            plot.data.all.taxa.sites <- rbind(plot.data.all.taxa.sites, temp.plot.data)
            # list.plot.data.biom.rates.fact[[taxon]][[site]]
        }
    }
    
    return(list(list.df.add.res.taxa = list.df.rates.limfact.taxa,
                plot.data = plot.data.all.taxa.sites))
}

