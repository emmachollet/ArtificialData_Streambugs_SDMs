## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## --- utilities needed to preprocess data and run Streambugs  ---
## --- to create synthetic data ---
##
## --- October 2023 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




print.pdf.plots <- function(list.plots, width = 12, height = width*3/4, dir.output, info.file.name = "", file.name, png = FALSE, png.square = F, png.vertical = F, png.ratio = 1){
    
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

# exp.transform <-   function(x,intercept=0,curv=0)
# {
#   #!if curv > 0 and intercept <1: function is curved to the right, if curv < 0 and intercept <1 function is curved to the left
#   #!if curv > 0 and intercept >1: function is curved to the left,  if curv < 0 and intercept >1 function is curved to the right
#   
#   if(curv == 0)
#   { 
#     y = intercept-(intercept-1)*x
#   } else {
#     y = intercept - (intercept -1) * (1 - exp(-curv * x)) / (1-exp(-curv)) 
#   }
#   return(y)
# }


# small basic ggplot colors function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

## Get taxonomic level ####

get.df.prev.pres.taxonomy <- function(data.inv, data.taxonomy, catch.variable, vect.catch.select, dir.outputs){
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
    print.pdf.plots(list.plots = list.plots, width = 12, height = 8, dir.output = dir.outputs, info.file.name = "", file.name = file.name, 
                    png = F)
    
    
    return(list.df.prev.pres)
    
}


write.df.preferences <- function(system.def, Invertebrates, file.name, dir.outputs){
    
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
      # n <- 2
        env.pref <- names(list.preferences)[n]
        cat(env.pref, "classes :", list.preferences[[n]][,1], "\n",
            env.pref, "values :", list.preferences[[n]][,2], "\n")
        classes.pref <- list.preferences[[n]][,1]
        classes.val <- list.preferences[[n]][,2]
        taxa.pref <- parameters[which(grepl(env.pref, names(parameters)))]
        df.pref <- data.frame(Classes = classes.pref, Values = classes.val)
        for (taxon in Invertebrates) {
          # taxon <- Invertebrates[1]
            df.pref[,taxon] <- NA
            for (c in classes.pref) {
              # c <- classes.pref[1]
                df.pref[which(df.pref$Classes == c), taxon] <-
                    taxa.pref[which(grepl(taxon, names(taxa.pref)) & grepl(paste0(c, "$"), names(taxa.pref)))] # $ Asserts that we are at the end.
            }
        }
        n.taxa <- length(Invertebrates)
        csv.name <- paste0(dir.outputs, file.name, n.taxa, "Taxa_", env.pref, ".csv")
        cat("Writing:", csv.name, "\n")
        write.csv(df.pref, file = csv.name)
        
        list.df.preferences[[env.pref]] <- df.pref
    }
    
    return(list.df.preferences)
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


# Figure SI A 7-15: spatial distribution env. fact. ####
maps.env.fact <- function(map.inputs, vect.env.fact, data){
    
    # data <- na.omit(data)
    
    plot_maps_variable = function (pair.fact.name, map.inputs, data) {
        variable <- pair.fact.name[1]
        print(variable)
        name.variable <- pair.fact.name[2]
        if(is.na(name.variable)){name.variable <- variable}
        print(name.variable)
        # data[,"ReachID"] <- gsub("CSCF.CH.", "", data[,"ReachID"])
        ggplot(data = data[,c("X","Y", "ReachID", variable)]) + 
            geom_sf(data = map.inputs$ch, fill="#E8E8E8", color="black") + 
            geom_sf(data = map.inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE) + 
            geom_sf(data = map.inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE) + 
            geom_point(aes(x=X, y=Y, color = data[, variable]), size= 5, alpha= 0.8) + 
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
    }
    
    # !! almost there, just missing labels 
    list.plots <- lapply(vect.env.fact, plot_maps_variable, map.inputs  = map.inputs, data = data)
    
    return(list.plots)
}



# process additional results ####
get.plot.data.add.res <- function(catch.results, list.factors.type){
    
    # catch.results <- list.results$RheinabBS
    
    df.res.add <- as.data.frame(catch.results$res.add)
    col.df.res.add <- colnames(df.res.add)
    temp.vect.taxa <- catch.results$y.names$taxa[-c(1:2)]
    temp.vect.sites <- catch.results$y.names$reaches
    y.names <- catch.results$y.names$y.names
    streambugs.results <- as.matrix(catch.results$res)
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
            c.ind2 <- which(grepl(paste0(site, "_"), colnames(temp.df)))
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

# 
# # if additional results calculated
# if(res.add){
#     
#     # extract results
#     df.res.add <- as.data.frame(res.catch$res$res.add)
#     colnames(df.res.add)[which(colnames(df.res.add) == "time")] <- "Time"
#     # str(df.res.add)
#     
#     # extract rates and limitation factors value over time (only if res.add = T)
#     col.df.res.add <- colnames(df.res.add)
#     vect.factor <- col.df.res.add
#     for (variable in y.names) {
#         vect.factor <- gsub(variable, "", vect.factor)
#     }
#     vect.factor <- unique(vect.factor)
#     print(vect.factor)
#     vect.factor <- vect.factor[2:21]
#     
#     # [1] "time"                       "r_basal."                   "r_miner."
#     # [4] "r_resp."                    "fsapro."                    "forgmicropoll."
#     # [7] "r_death."                   "flimI."                     "flimP."
#     # [10] "flimN."                     "flimnutrients."             "fselfshade."
#     # [13] "r_prod."                    "fcurrent."                  "ftempmax."
#     # [16] "fmicrohab."                 "sum_food."                  "sum_food_pref."
#     # [19] "fselfinh."                  "ffoodlim."                  "r_cons_tot."
#     
#     # based on result above, define vectors of rates and factors per group
#     # factor types
#     list.factors.type <- list(
#         "Biomass" = c("Biomass"),
#         "Rates" = c("r_basal", "r_resp", "r_prod", "r_death", "r_cons_tot", "r_miner"),
#         "DeathFactors" = c("fsapro", "forgmicropoll"),
#         "AlgaeFactors" = c("flimI", "flimP", "flimN",
#                            "fshade", "flimnutrients", "fselfshade"),
#         "InhibFactors" = c("fmicrohab", "fcurrent",
#                            "ftempmax",  "fselfinh", "ffoodlim"),
#         "FoodFactors" = c("sum_food", "sum_food_pref"))
#     fact.types <- names(list.factors.type)
#     
#     # assign colors to rates for better plotting
#     mypalette.8.col <- scales::hue_pal()(8)
#     scales::show_col(mypalette.8.col)
#     
#     vect.colors <- c()
#     set.seed(3)
#     for (n in 1:length(list.factors.type)) {
#         # n <- 2
#         # names(list.factors.type[[n]]) <- gg_color_hue(length(list.factors.type[[n]]))
#         names(list.factors.type[[n]]) <- sample(mypalette.8.col, length(list.factors.type[[n]]))
#         vect.colors <- c(vect.colors, list.factors.type[[n]])
#     }
#     mycolors = names(vect.colors)
#     scales::show_col(mycolors)
# }
# 
# temp.vect.taxa <- vect.taxa[3:length(vect.taxa)]
# print(temp.vect.taxa)
# temp.vect.sites <- vect.sites
# 
# if(res.add){
#     
#     
#     # make list of data frames with additional results per taxon per site
#     list.df.rates.limfact.taxa <- list()
#     for (taxon in temp.vect.taxa) {
#         # taxon <- vect.taxa[3]
#         c.ind <- which(grepl(taxon, col.df.res.add) & !grepl("__", col.df.res.add))
#         temp.df <- df.res.add[,c(1, c.ind)]
#         temp.list <- list()
#         for (site in temp.vect.sites) {
#             # site <- temp.vect.sites[1]
#             c.ind2 <- which(grepl(site, colnames(temp.df)))
#             temp.df2 <- temp.df[,c(1,c.ind2)]
#             c.ind3 <- which(grepl(taxon, colnames(streambugs.results)) & grepl(site, colnames(streambugs.results)))
#             temp.df2$Biomass <- streambugs.results[,c.ind3]
#             temp.list[[site]] <- temp.df2
#         }
#         list.df.rates.limfact.taxa[[taxon]] <- temp.list
#     }
#     
#     # make long dataframe for plotting additional result
#     plot.data.all.taxa.sites <- data.frame()
#     for (taxon in temp.vect.taxa) {
#         # taxon <- vect.taxa[3]
#         temp.list.df.sites <- list.df.rates.limfact.taxa[[taxon]]
#         # list.plot.data.biom.rates.fact[[taxon]] <- list()
#         
#         for (site in temp.vect.sites) {
#             # site <- temp.vect.sites[1]
#             temp.df.rates.limfact <- temp.list.df.sites[[site]]
#             ind.y.names <- which(grepl(taxon, y.names) & grepl(site, y.names))
#             state.variable <- y.names[ind.y.names]
#             colnames(temp.df.rates.limfact) <- gsub(paste0(".", state.variable), "", colnames(temp.df.rates.limfact))
#             if(grepl("Algae", taxon)){ # if limfactor of algae. add 1 - shade
#                 temp.df.rates.limfact$fshade <- 1 - env.data[which(env.data$ReachID == site), "shade"]
#             }
#             
#             temp.plot.data <- gather(temp.df.rates.limfact, key = Factor, value = Value, -c("Time"))
#             
#             temp.plot.data$Taxon <- taxon
#             temp.plot.data$Site <- site
#             
#             # add type of limitation factor for easier plotting later
#             temp.plot.data$FactorType <- NA
#             for (type in fact.types) {
#                 # type <- fact.types[1]
#                 r.ind <- which(temp.plot.data$Factor %in% list.factors.type[[type]])
#                 temp.plot.data[r.ind, "FactorType"] <- type
#             }
#             plot.data.all.taxa.sites <- rbind(plot.data.all.taxa.sites, temp.plot.data)
#             # list.plot.data.biom.rates.fact[[taxon]][[site]]
#         }
#     }
# } 
# 
# # select one taxon for analysis
# taxon <- vect.taxa[3]
# Invertebrates[7]
# 
# if(res.add){
#     
#     # assign colors to rates for better plotting
#     mypalette.8.col <- scales::hue_pal()(8)
#     scales::show_col(mypalette.8.col)
#     
#     vect.colors <- c()
#     set.seed(3)
#     for (n in 1:length(list.factors.type)) {
#         # n <- 2
#         # names(list.factors.type[[n]]) <- gg_color_hue(length(list.factors.type[[n]]))
#         names(list.factors.type[[n]]) <- sample(mypalette.8.col, length(list.factors.type[[n]]))
#         vect.colors <- c(vect.colors, list.factors.type[[n]])
#     }
#     mycolors = names(vect.colors)
#     scales::show_col(mycolors)
#     
#     # plot all rates and limitation factors
#     list.plots <- list()
#     for (taxon in temp.vect.taxa) {
#         plot.data <- plot.data.all.taxa.sites %>%
#             filter(Taxon == taxon) # %>%
#         # filter(Site != 107)
#         p <- ggplot(plot.data, aes(x = Time, y = Value, color = Factor))
#         p <- p + geom_line(size=1)
#         p <- p + facet_grid(FactorType ~ Site, scales = "free")
#         if(!grepl("Algae", taxon)){ p <- p + scale_color_manual(values = mycolors)}
#         p <- p + theme_bw()
#         p <- p + labs(title = taxon)
#         ggplotly(p)
#         
#         list.plots[[taxon]] <- p
#     }
#     
#     file.name <- "_AddRes_RatesLimFact"
#     print.pdf.plots(list.plots = list.plots, width = 23, height = 8, dir.output = dir.outputs, info.file.name = name.run, file.name = file.name, 
#                     png = F)
# }
