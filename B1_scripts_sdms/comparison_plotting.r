# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Init ----
Sys.setenv(LANG="EN")
set.seed(13)   # Always set seed to a lucky number
getwd()        # Show working directory. It needs to be the location of 'main.r'
rm(list=ls())  # Free work space
graphics.off() # Clean graphics display


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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

if (!require("jsonlite")){install.packages("jsonlite"); library("jsonlite")}



# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Custom R scripts ----
source("ml_models.r")
source("plotting.r")
source("performances_assessment.r")
source("global_variables.r")


# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Creating the plots for the different experiments ----


# gauss_temp DONE
list.exp.gauss      <- list("baseline"          = "baseline_13_03_2024_15h05",
                            "low noise"         = "gaussian_temp_1_14_03_2024_08h26",
                            "mid noise"         = "gaussian_temp_3_14_03_2024_08h31",
                            "high noise"        = "gaussian_temp_5_14_03_2024_08h36")

create.comparison.plots("gauss_temp", list.exp.gauss, 
                        file.prev.taxa="8Catch_1416Sites_curve.curr.temp-10_interc.orgmic.sapro4_10Yea_3651Steps_PrevalenceTaxonomy_ThreshPresAbs.csv",
                        taxon="Gammaridae",
                        env.fact="tempmaxC")



# miss_gammaridae
list.exp.miss.gamma <- list("baseline"          = "baseline_21_07_2023_12h00",
                            "low noise"         = "gammaridae_missdetection_10_22_07_2023_04h08",
                            "mid noise"         = "gammaridae_missdetection_25_22_07_2023_10h15",
                            "high noise"        = "gammaridae_missdetection_50_22_07_2023_16h37")

create.comparison.plots("miss_gammaridae", list.exp.miss.gamma)



# miss_all_taxa
list.exp.miss.all   <- list("baseline"          = "baseline_21_07_2023_12h00",
                            "low noise"         = "all_taxa_missdetection_10_22_07_2023_23h01",
                            "mid noise"         = "all_taxa_missdetection_25_23_07_2023_05h43",
                            "high noise"        = "all_taxa_missdetection_50_23_07_2023_12h09")

create.comparison.plots("miss_all_taxa", list.exp.miss.all)



# subset
list.exp.subset     <- list("baseline"          = "baseline_21_07_2023_12h00",
                            "low noise"         = "subset_2000_24_07_2023_00h40",
                            "mid noise"         = "subset_1000_23_07_2023_21h24",
                            "high noise"        = "subset_500_23_07_2023_19h25")

create.comparison.plots("subset", list.exp.subset)



# remove_fact
list.remove.fact    <- list("baseline"          = "baseline_21_07_2023_12h00",
                            "low noise"         = "env_fact_removed_1_24_07_2023_03h30",
                            "mid noise"         = "env_fact_removed_4_24_07_2023_06h16",
                            "high noise"        = "env_fact_removed_8_24_07_2023_08h30")

create.comparison.plots("remove_fact", list.remove.fact)



# add_fact
#list.add.fact       <- list("baseline"          = "baseline_21_07_2023_12h00",
#                            "low noise"         = "add_env_fact_1_24_07_2023_10h41",
#                            "mid noise"         = "TODO",
#                            "high noise"        = "TODO")

#create.comparison.plots("add_fact", list.add.fact)



# gauss_temp_vel
list.gauss.temp.vel <- list("baseline"          = "baseline_21_07_2023_12h00",
                            "low noise"         = "gaussian_velocity_3_24_07_2023_23h09",
                            "mid noise"         = "gaussian_temp_3_21_07_2023_16h42",
                            "high noise"        = "gaussian_temp_velocity_3_25_07_2023_06h38")

create.comparison.plots("gauss_temp_vel", list.gauss.temp.vel)
  
