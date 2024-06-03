# Load libraries:
# ===============

if ( !require(streambugs) )
{
  install.packages("streambugs")
  library(streambugs)
}

if ( !require(deSolve) )
{
  install.packages("deSolve")
  library(deSolve)
}

if ( !require(dplyr) )
{
  install.packages("dplyr")
  library(dplyr)
}

# if ( !require(vioplot) )
# {
#   install.packages("vioplot")
#   library(vioplot)
# }

# load streambugs core
# --------------------

# path.core <- "library/core/"
# 
# source(paste(path.core,"streambugs.r",sep=""))
# source(paste(path.core,"streambugs_rhs.r",sep=""))
# source(paste(path.core,"streambugs_aux.r",sep=""))
# source(paste(path.core,"get_par_trait_cond.r",sep=""))
# 
# rm(path.core)

# load streambugs application files

#source("library/construct_statevariables.r")
source("library/construct_envpars.r")
source("library/sysanal.r")
source("library/construct_invpars_traits.r")
source("library/construct_invpars_traits.r")
source("library/construct_traitpars.r")
source("library/construct_traitpars_species.r")
source("library/construct_traitpars_genus.r")
source("library/construct_traitpars_family.r")
source("library/scale_traits.r")
source("library/derive_par_env_traits.r")

source("library/define_par_dist.r")
source("library/get_par_dist_def.r")
source("library/generate_pars_samp_matrix.r")
source("library/convert_CSusPOM.r")

source("library/stoichcalc_function.r")
source("library/assign_par_stoich_p.r")

source("library/plot_foodweb.r")
source("library/plot_streambugs_unc.r")

source("library/construct_posterior_parsamp.r")
source("library/construct_matrix_res_add_t.r")
source("library/get_res_add_t.r")

source("library/divide_parameters.r")
source("library/transform_priors.r")

source("library/generate_par_samp_vector.r")

source("library/calc_prior_dens.r")


# application specific files: 

# source("library/application specific NIS/boxplot_freqobs_pobs.r")
# source("library/application specific NIS/boxplot_freqobs_pobs_binary.r")
source("library/application specific NIS/check_observations.r")
source("library/application specific NIS/calc_max_correct_modres.r")
source("library/application specific NIS/compare_observed_predicted.r")
source("library/application specific NIS/compare_observed_predicted_sample.r")
source("library/application specific NIS/get_factors_inverts_wrong_predicted.r")
source("library/application specific NIS/calc_D_crit.r")
source("library/application specific NIS/check_survival.r")

source("library/application specific NIS/calc_loglikelihood.r")
source("library/application specific NIS/calc_max_correct_modres.r")
source("library/application specific NIS/calc_res_steadystate.r")
source("library/application specific NIS/calc_freq_observation.r")

source("library/application specific NIS/simulate_and_return_likeli_moresites.r")
source("library/application specific NIS/simulate_and_return_likeli_catchments.r")
source("library/application specific NIS/simulate_streambugs_return_resadd.r")
source("library/application specific NIS/simulation_parsamp.r") 
# source("library/application specific NIS/mcmc_simulation_caseE.r") 
# source("library/application specific NIS/mcmc_caseE_parallel.r")
# source("library/application specific NIS/mcmc_simulation_caseE_parallel.r")


