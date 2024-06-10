## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##
## --- utilities needed to preprocess data and run Streambugs  ---
## --- to create synthetic data ---
##
## --- October 2023 -- Emma Chollet ---
##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get the state variables for a specific catchment
construct.variables.par.catch <- function(catch, data.env.inputs, 
                                          data.par.update, par.adjust,
                                          data.taxa.selection, taxa.selection, selected.taxa.analysis,
                                          sites.selection, select.taxonomy, 
                                          catch.variable,
                                          plot.foodweb = F, name.run,
                                          dir.inputs, dir.outputs){
    # catch = vect.catch.select[1]
    cat("\nConstructing variables and parameters for catchment:\n", 
        catch, "\n")

    list.metadata.catch <- list()
    list.metadata.catch[["Catchment"]] <- catch
    
    # environmental data for the catchment
    if(catch == taxa.selection["metric"]){
        env.data <- data.env.inputs[1,]
        threshold <- 0
    } else {
        env.data <- data.env.inputs[which(data.env.inputs[,catch.variable] == catch),] # filter sites from a specific catchment
        if(sites.selection["select.sites"] == "random"){
            if(sites.selection["n.sites"] != "all" & as.numeric(sites.selection["n.sites"]) < dim(env.data)[1]){
                env.data <- env.data[1:sites.selection["n.sites"],] # temporary take only few sites for trial
            }
        } else {
            rind.sites.selected <- which(grepl(sites.selection["select.sites"], env.data$ReachID))
            if(sites.selection["n.sites"] != length(rind.sites.selected)){
                rind.sites.selected <- c(1:(n.sites - length(rind.sites.selected)), rind.sites.selected) 
            }
            env.data <- env.data[rind.sites.selected,] # temporary take only few sites for trial
        }
        threshold <- taxa.selection["threshold"]
    }

    Reaches   <- env.data$ReachID
    Habitats <- env.data$Habitat
    list.metadata.catch[["Reaches"]] <- Reaches
    list.metadata.catch[["Habitats"]] <- Habitats
    
    # taxa pool of the catchment
    Invertebrates <- data.taxa.selection[which(data.taxa.selection$Taxonomic.level %in% select.taxonomy &
                                                    data.taxa.selection[, catch] > threshold), "Taxon"]
    list.metadata.catch[["Invertebrates"]] <- Invertebrates
    
    # Pom and algae
    POM   <- c("FPOM","CPOM")
    Algae <- "Algae"

    # construct variable name
    y.names <- construct.statevariables(POM=POM, Algae=Algae, Invertebrates=Invertebrates,
                                        Reaches=unique(Reaches), Habitats=unique(Habitats))
    y.names <- decode.statevarnames(y.names)
    
    # construct parameters
    
    # reach and habitat specific parameters - environmental inputs
    # w L T I0 fshade P N
    par.env <- construct.envpars(Reaches, Habitats, env.data) #xxx nis: this functions wants some microhabitats, we could change it to avoid this
    
    # Invertebrates trait parameters
    
    # databases to be used for deriving trait parameters
    file.db.feeding    = paste0(dir.inputs, "databases/db_freshwaterecology_20131206_AP.dat")
    file.db.sapro      = paste0(dir.inputs, "databases/db_freshwaterecology_20131206_AP.dat")
    file.db.spear      = paste0(dir.inputs, "databases/db_traits_rspear_Eurasia_20130114_AP.dat")
    #file.db.microhab   = paste0(dir.inputs, "databases/db_traits_tachet_AP.dat") 	#xxx nis: not needed for our study
    file.db.current    = paste0(dir.inputs, "databases/db_traits_tachet_AP.dat")
    # file.db.temp       = paste0(dir.inputs, "databases/db_traits_tachet_AP.dat")
    file.db.temp      = paste0(dir.inputs, "databases/db_freshwaterecology_20131206_AP.dat")
    
    file.db.bodymass   = paste0(dir.inputs, "databases/db_bodymass_meas_aquaplus.dat")
    file.db.taxonomy   = paste0(dir.inputs, "databases/db_taxonomy.dat")
    
    # mid values of classes of environmental conditions used for traits 
    par.env.global <- derive.par.env.traits(file.db.temp=file.db.temp,
                                            file.db.current=file.db.current,
                                            file.db.sapro=file.db.sapro,
                                            file.db.spear=file.db.spear)  
    
    par.invtraits <- construct.invpars.traits(Invertebrates,
                                              file.db.taxonomy=file.db.taxonomy,
                                              file.db.feeding=file.db.feeding,
                                              #file.db.microhab=file.db.microhab, #xxx nis: not needed for our study
                                              file.db.current=file.db.current,
                                              file.db.temp=file.db.temp,
                                              file.db.sapro=file.db.sapro,
                                              file.db.spear=file.db.spear,
                                              file.db.bodymass=file.db.bodymass)
    par.invtraits.orig <- par.invtraits
    # # extract name of traits from streambugs function
    # name.traits <- names(par.invtraits)
    # for (taxon in Invertebrates) {
    #     name.traits <- gsub(paste0(taxon, "_"), "", name.traits)
    # }
    # name.traits <- unique(name.traits)
    
    if(par.adjust[["update.traits"]][["Flag"]]){
      cat("Using updated (Vermeiren's) invertebrate ecological preferences only for selected taxa.")
      # ecr/nis 18.12.23: overwrite here some preference traits from Vermeiren et al. 2021 for selected taxa
      # see table: T1d_BDM_CH_Tvsappestsubst_maxpost_trait_pars_2023-12-15.dat
      # colnames(data.par.update)[2:9]
      names.selected.taxa <- gsub("Occurrence.", "", selected.taxa.analysis)

      for (taxon in Invertebrates) {
        # taxon <- data.par.update$Taxon[2]
        # taxon <- Invertebrates[1]
          if(taxon %in% data.par.update$Taxon && taxon %in% names.selected.taxa){
            cat("\nUpdating temperature and current preference for taxon:", taxon)
              for (trait in colnames(data.par.update)[2:9]) { # update only temp and current trait
                  # trait <- colnames(data.par.update)[2]
                  ind <- which(grepl(taxon, names(par.invtraits)) & grepl(trait, names(par.invtraits)))
                  # ind <- 1073
                  if(length(ind) > 0){
                      par.invtraits[ind] <- data.par.update[which(data.par.update$Taxon == taxon), trait]
                  }
              }
          } else {
            # cat("\nFollowing taxa in our taxa list but not in the Vermeiren taxa list: ")
            #   cat(taxon, " ")
          }
      }
    } else {
      cat("Using original invertebrate ecological preferences.")
    }
    
    # try to round low preferences to 0 to get stronger responses
    ind.current <- which(grepl("currenttolval", names(par.invtraits)))
    ind.temp <- which(grepl("tempmaxtolval", names(par.invtraits)))
    if(par.adjust[["round.inv.traits"]][["Flag"]]){
      for(i in c(ind.current, ind.temp)){
        # print(i)
        # i <- 780
        # pref <- par.invtraits[i]
        # print(par.invtraits[i])
        thresh.round <- par.adjust[["round.inv.traits"]][["Value"]]
        if(par.invtraits[i] < thresh.round){
          print("COUCOU")
          par.invtraits[i] <- 0
          print(par.invtraits[i])
        }
      }
    }

    # set dummy values with habitat suitability 1 for all invertebrates and all types
    par.invtraits.mh <- rep(1,4*length(Invertebrates))
    names(par.invtraits.mh) <- c(paste0(Invertebrates,"_microhabtolval_type",1),
                                 paste0(Invertebrates,"_microhabtolval_type",2),
                                 paste0(Invertebrates,"_microhabtolval_type",3),
                                 paste0(Invertebrates,"_microhabtolval_type",4))
    
    par.invtraits <- c(par.invtraits,par.invtraits.mh)
    
    # other parameters
    
    # read in definitions for parameter distributions from file
    par <- sysanal.read.distdef(paste0(dir.inputs, "parameter_input_fix.dat"))
    
    # adjust curve and intercept of limitation factor function
    if(par.adjust[["curve.curr.temp"]][["Flag"]]){
        par$ftempmax_curv[[2]] <- par.adjust[["curve.curr.temp"]][["Value"]]
        par$fcurrent_curv[[2]] <- par.adjust[["curve.curr.temp"]][["Value"]]
    }
    if(par.adjust[["curve.orgmic.sapro"]][["Flag"]]){
      par$forgmicropoll_curv[[2]] <- par.adjust[["curve.orgmic.sapro"]][["Value"]]
      par$fsapro_curv[[2]] <- par.adjust[["curve.orgmic.sapro"]][["Value"]]
    }
    if(par.adjust[["interc.orgmic.sapro"]][["Flag"]]){
      par$forgmicropoll_intercept[[2]] <- par.adjust[["interc.orgmic.sapro"]][["Value"]]
      par$fsapro_intercept[[2]] <- par.adjust[["interc.orgmic.sapro"]][["Value"]]
    }
    
    if(par.adjust[["hdens"]][["Flag"]]){
      par$Invertebrates_hdens[[2]] <- par.adjust[["hdens"]][["Value"]] * as.numeric(par$Invertebrates_hdens[[2]])
    }
    
    # automatic generated parameters for organisms:
    for ( i in c(Invertebrates,"Algae") ){
        par[[paste(i,"_fgrotax",sep="")]]   <- par[["Organism_fgrotax"]]    # -
        par[[paste(i,"_fbasaltax",sep="")]] <- par[["Organism_fbasaltax"]]  # -
    }
    
    # assign taxon dependent hdens and Kfood parameters: 
    par.feedtype <- par.invtraits[grepl("feedingtype",names(par.invtraits))]
    
    for ( i in c(Invertebrates) ){
        par[[paste(i,"_hdens",sep="")]] <- par[["Invertebrates_hdens"]]   
        
        # test if pure filter feeders:
        # because food is per m3, need to be converted to density and change preference factor
        if(file.db.feeding==paste0(dir.inputs, "databases/db_freshwaterecology_20131206_AP.dat")) {
            if( sum(par.feedtype[paste(i,"feedingtype_aff",sep="_")]+
                    par.feedtype[paste(i,"feedingtype_pff",sep="_")])>0 & 
                sum(par.feedtype[paste(i,"feedingtype_shr",sep="_")],
                    par.feedtype[paste(i,"feedingtype_min",sep="_")],
                    par.feedtype[paste(i,"feedingtype_xyl",sep="_")],
                    par.feedtype[paste(i,"feedingtype_gra",sep="_")],
                    par.feedtype[paste(i,"feedingtype_gat",sep="_")],
                    par.feedtype[paste(i,"feedingtype_pre",sep="_")])==0)
            {
                par[[paste(i,"_Kfood",sep="")]] <- par[["Filt_Kfood"]] 
                cat(i ,"is pure filterfeeder, Kfood is set to",par[["Filt_Kfood"]],"\n" )
            } else {par[[paste(i,"_Kfood",sep="")]] <- par[["Invertebrates_Kfood"]]  }
        } else {
            if(file.db.feeding==paste0(dir.inputs, "databases/db_traits_tachet.dat")) {
                if( sum(par.feedtype[paste(i,"feedingtype_filter.feeder",sep="_")]) >0 & 
                    sum(par.feedtype[paste(i,"feedingtype_deposit.feeder",sep="_")],
                        par.feedtype[paste(i,"feedingtype_shredder",sep="_")],
                        par.feedtype[paste(i,"feedingtype_scraper",sep="_")],
                        par.feedtype[paste(i,"feedingtype_piercer",sep="_")],
                        par.feedtype[paste(i,"feedingtype_predator",sep="_")])==0)
                {
                    par[[paste(i,"_Kfood",sep="")]] <- par[["Filt_Kfood"]] 
                    cat(i ,"is pure filterfeeder, Kfood is set to",par[["Filt_Kfood"]],"\n" )
                } else { par[[paste(i,"_Kfood",sep="")]] <- par[["Invertebrates_Kfood"]]  }
                
            } else { warning("file.db.feeding not recognized:", file.db.feeding)}
        }
    }
    
    # remove parameters that are no longer needed
    ind_rm <- which(names(par)=="Organism_fgrotax"|names(par)=="Organism_fbasaltax"|
                        names(par)=="Invertebrates_Kfood"|names(par)=="Invertebrates_hdens"|
                        names(par)=="Filt_Kfood")
    par <- par[-ind_rm]
    
    # fixed stoichiometric parameters
    par.stoich.fix <- numeric(0)
    
    # generate stoichiometric parameters that are fixed (mineralization, respiration, production):
    taxa.POM <- unique(y.names$y.taxa[y.names$y.groups=="POM"])                         # POM:
    par.stoich.fix[paste("Miner",taxa.POM,taxa.POM,sep="_")] <- -1                      # mineralization
    taxa.Algae <- unique(y.names$y.taxa[y.names$y.groups=="Algae"])                     # Algae:
    par.stoich.fix[paste("Resp",taxa.Algae,taxa.Algae,sep="_")]  <- -1                  # respiration
    par.stoich.fix[paste("Prod",taxa.Algae,taxa.Algae,sep="_")]  <- 1                   # production
    taxa.Invertebrates <- unique(y.names$y.taxa[y.names$y.groups=="Invertebrates"])     # Invertebrates:
    par.stoich.fix[paste("Resp",taxa.Invertebrates,taxa.Invertebrates,sep="_")]  <- -1  # respiration
    
    # xxx ecr: could be simplified to fix par directly
    # combine parameters
    par.unc <- c(par,
                 par.env,
                 par.env.global,
                 par.invtraits,
                 par.stoich.fix)
    
    
    # Generate parameter sample
    par.samp <- generate.par.samp.matrix(n.samp=1, par.unc=par.unc)
    class(par.samp) <- "numeric"
    
    # write.table(par.samp,paste(dir.catch,"/par_samp_",name.run,".dat",sep=""),
    #             sep="\t", row.names = T, col.names=NA)
    
    # Set up time dependant input (to NA for our application case)
    # ecr: problem if this is only in the function
    # !! now "inp" = NA globally, but should be defined in this function, but there is an error
    # Error in streambugs:::get.inpind.parval.taxaprop.traits(trait.names = "feedingtype", : 
    #                                                             object 'inp' not found    
    # inp <- NA
    
    par.fix <- par.samp[,1]
    names(par.fix) <- rownames(par.samp)
    
    # calc stoichiometric parameters
    par.stoich.out <- calc.stoich(par=as.list(par.fix),returns="parout")
    
    # assign calculated group stoichiometric parameters to the taxa
    par.stoich.taxa <- assign.par.stoich(par.invtraits,par.stoich.out,y.names)
    
    # combine parameters
    par.fixc <- c(par.fix, par.stoich.taxa)
    
    # convert CSusPOM to DSusPOM
    # new: we convert CSusPOM to DSusPOM just before we call the model,
    # in this way Filt_scope and CSusPOM can both be certain or uncertain 
    # parameters independently of each other.
    par.fixc <- convert.CSusPOM(par.fixc)
    
    
    # save results in catch folder for other use (e.g., Bayesian inference script)
    dir.catch <- paste0(dir.outputs, catch, "_", sites.selection["n.sites"], "Sites", "/")
    dir.create(dir.catch)
    
    if(plot.foodweb == T){
        # plot foodweb
      cat("Printing foodweb in pdf.")
        pdf(paste(dir.catch, catch, "_", length(Invertebrates), "Taxa_", "Foodweb",
                  # "feedtypesFWB_R_",
                  #par.fix["ratio_pred_prey"],
                  ".pdf",sep=""),
            width=9,height=5,onefile=T)
        
        foodweb.plot(y.names,par=par.fixc,cex=0.7,font=1,title="complete foodweb r3",ncrit=8,
                     lcrit=15,lwd=0,bg=NA,lcol=colors()[555]) #,texts=F,pointcol=T)
        
        foodweb.plot(y.names,par=par.fixc,cex=0.9,font=2,title="complete foodweb",ncrit=8,
                     lcrit=15,lwd=0,lcol=1,bg=NA,texts=F,pointcol=T)
        dev.off()
    }
    
    list.results <- list("y.names" = y.names, "par.fixc" = par.fixc, 
                         "env.data" = env.data, "Invertebrates" = Invertebrates, 
                         "par.invtraits.orig" = par.invtraits.orig, "par.unc" = par.unc,
                         "list.metadata.catch" = list.metadata.catch)
    file.name <- paste0(dir.catch, catch, "_list.inputs.variables.parameters.rds")
    saveRDS(list.results, file = file.name)
    return(list.results)
}


run.streambugs.catch <- function(y.names.par.catch, tout, res.add = F, name.run, dir.output, run.C = T, write.plot.results = F, ...){
    
    cat("\nRunning streambugs for:\n", 
        y.names.par.catch$list.metadata.catch$Catchment, 
        "with",  length(y.names.par.catch$Invertebrates) ,"taxa in",
        length(y.names.par.catch$y.names$reaches), "reaches\n\n")
    # res.add <- T
    # y.names.par.catch <- list.variables.par.catch[[1]]
    y.names <- y.names.par.catch$y.names
    par.fixc <- y.names.par.catch$par.fixc
    reaches.orig <- y.names.par.catch$y.names$reaches
    list.warnings <- list()
    vect.na.reaches <- c()
    
    # simulate reaches separately: ####
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~-

    start <- proc.time()

    for ( i in 1:length(y.names$reaches)){
        # i <- 1
      reach <- y.names$reaches[i]
      ind.reach <- y.names$y.reaches==reach

      cat("Running Streambugs for reach:", reach,"\n")
      
      res.i <- run.streambugs(y.names = y.names$y.names[ind.reach],
                              times   = tout,
                              par     = par.fixc,
                              inp     = NA, # ecr modified this for this application
                              C       = run.C,
                              # method = "euler",
                              return.res.add = res.add, # takes time, set to false if we only want simulations without additional output
                              # file.def = paste("output/Toss/streambugs_modelstruct_",name.run,".dat",sep=""),
                              # file.res = paste("output/Toss/streambugs_results_",name.run,".dat",sep=""),
                              # file.add = paste(dir.output, "/add_output_", name.run, ".dat", sep=""),
                              verbose = F #,
                              # atol = 1e-08,
                              # rtol = 1e-08
                              # ...
                              )$res

      # check if NA or not in results
      if(any(is.na(res.i))){ 
        
        # if there are any NA in results, write warning in metadata.file to be printed
        vect.na.reaches <- append(vect.na.reaches, reach)
        na.taxa <- gsub(paste0(reach, "_random_"), "", names(which(colSums(is.na(res.i))>0)))
        na.warning <- paste0("WARNING: NA in results for reach ", reach)
        # list.warnings[[reach]] <- list(na.warning, "Taxa" = na.taxa)
        list.warnings[[reach]] <- na.warning
        cat("\n", na.warning, "\n")
          
        # remove reach from parameters (update y.names below)
        ind.reach.par <- which(grepl(reach, names(par.fixc)))
        par.fixc.update <- par.fixc[-ind.reach.par]
        y.names.par.catch$par.fixc <- par.fixc.update # overwrite parameters in global y.names.par variable
        
        # remove reach from environmental data
        ind.reach.env.data <- which(y.names.par.catch$env.data$ReachID == reach)
        env.data.update <- y.names.par.catch$env.data[-ind.reach.env.data,]
        y.names.par.catch$env.data <- env.data.update # overwrite env.data in global y.names.par variable
        
      } else {
        # if no NA, bind results together
        if(!exists("res")){ # test if object "res" already exists
          res <- res.i
        } else{
          # can't use "left_join" because res is class "streambugs"
          # res <- left_join(res, res.i, by = "time")
          res <- cbind(res,res.i) # original code, advantage use base package, but can't join results with NAs
        }  
      }
    }
    
    # remove duplicated columns "time"
    ind.time <- which(colnames(res)=="time")
    ind.time <- ind.time[-1]
    res      <- res[,-ind.time]

    duration <- proc.time()-start
    cat(duration,"\n")
    
    # if reaches with NA, update y.names without reaches causing problems/NAs
    if(length(vect.na.reaches) > 0) {
      # recover elements of y.names for this catchment
      POM <- y.names$taxa[1:2]
      Algae <- y.names$taxa[3]
      Habitats <- y.names$habitats
      inv.catch <- y.names.par.catch$Invertebrates
      reaches.update <- reaches.orig[-which(reaches.orig %in% vect.na.reaches)]
      
      # construct updated y.names
      y.names.update <- construct.statevariables(POM=POM, Algae=Algae, Invertebrates=inv.catch,
                                                 Reaches=reaches.update, Habitats=Habitats)
      y.names.update <- decode.statevarnames(y.names.update)
      
      # overwrite original y.names
      y.names.par.catch$y.names <- y.names.update
    }
    
    # # simulate all reaches together: ####
    # # ~~~~~~~~~~~~~~~~~~~~~~~~~~~-
    # res <- run.streambugs(y.names = y.names, #y.names$y.names,
    #                       times   = tout,
    #                       par     = par.fixc,
    #                       inp     = NA, # ecr modified this for this application
    #                       C       = run.C,
    #                       # method = "euler",
    #                       return.res.add = res.add, # takes time, set to false if we only want simulations without additional output
    #                       # file.def = paste("output/Toss/streambugs_modelstruct_",name.run,".dat",sep=""),
    #                       # file.res = paste("output/Toss/streambugs_results_",name.run,".dat",sep=""),
    #                       # file.add = paste(dir.output, "/add_output_", name.run, ".dat", sep=""),
    #                       verbose = TRUE)#,
    #                       # ...)
    # 
    # print("Streambugs simulation finished.")

    if(res.add){ print("Calculating additional results.")}
    if(write.plot.results){
        
        # # plot simulation results as time series
        # file.name <- paste(dir.output, name.run,"_IndivBiomass_StreambugsPlot.pdf",sep="")
        # print("Plots in file:")
        # print(file.name)
        # pdf(file.name,width=18,height=16)
        # plot(res,par.fixc,inp)
        # dev.off()
        
        # write simulation results as time series in csv file
        streambugs.results <- as.matrix(res)
        filename <- paste(dir.output, name.run, "_Results.csv", sep="")
        write.table(streambugs.results, filename, sep = ",", row.names = F)
    }
    
    y.names.par.catch[["res"]] <- res
    y.names.par.catch$list.metadata.catch[["simulation_time"]] <- duration
    y.names.par.catch$list.metadata.catch[["warnings"]] <- list.warnings
    
    return(y.names.par.catch)
    
}



