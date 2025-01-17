# Script to run Dengvaxia survival models, fitting to severe data --------------
#  - Each model variant (15 total) is the same stan model, with different
#    parameters turned on and off, using the flags.
# -  The final model (10) is given at the top, all other models are commented 
#    out but can be run below. 
#  - n_it defines the number of iterations to model is run for, set here to 1000
#    as a demonstration, but run for 10000 for all models presented in the 
#    paper. 
# ------------------------------------------------------------------------------


setwd("Chapter3")
rm(list = ls())
dir.create("CYD/output/severe")

# source functions 
file.sources = paste0("CYD/R/", list.files(path = "CYD/R/"))
sapply(file.sources, source)
n_it = 1000 # change to 10000 to recreate thesis results 

# M10 in thesis but started model runs from 0 
# M9 - M6 but serotype specific L ----------------------------------------------
run_CYD_model (
  L_K = 1,
  include_beta = 2,
  include_eps = 1,
  psi_J = 1,
  tau_K = 0,
  mono_lc_SN = 1,
  mono_lc_MU = 2,
  delta_KJ = 2,
  include_pK2 = 0,
  MU_test_SN = 1,
  folder = "severe/M9",
  n_it = n_it,
  adapt_delta = 0.9,
  stan_model = "CYD_model_severe.stan",
  baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
  cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
  mu =  read.csv("CYD/data/processed/mu.csv"),
  severe = T
)


# 
# # M0 - start with best fitting CYD-TDV model without severe data ---------------
# run_CYD_model (
#   mono_lc_SN = 1,
#   L_K = 1, 
#   tau_K = 0,
#   mono_lc_MU = 2,
#   include_beta = 0,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M0",
#   n_it = n_it,
#   adapt_delta = 0.8,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# # M1 - M1 but single L ---------------------------------------------------------
# run_CYD_model (
#   L_K = 0, 
#   mono_lc_SN = 1,
#   tau_K = 0,
#   mono_lc_MU = 2,
#   include_beta = 0,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M1",
#   n_it = n_it,
#   adapt_delta = 0.99,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# # M2 - M1 but age-group offset -------------------------------------------------
# run_CYD_model (
#   include_beta = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M2",
#   n_it = n_it,
#    adapt_delta = 0.99,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# # M3 - M1 but age-specific psi -------------------------------------------------
# run_CYD_model (
#   psi_J = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   include_beta = 0,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M3",
#   n_it = n_it,
#    adapt_delta = 0.99,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# # M4 - M3 but include epsilon --------------------------------------------------
# run_CYD_model (
#   include_eps = 1,
#   psi_J = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   include_beta = 0,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M4",
#   n_it = n_it,
#    adapt_delta = 0.99,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# 
# # M5 - M4 but include single beta ----------------------------------------------
# run_CYD_model (
#   include_beta = 1,
#   include_eps = 1,
#   psi_J = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M5",
#   n_it = n_it,
#    adapt_delta = 0.99,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# 
# # M6 - M5 but include outcome specific beta ------------------------------------
# run_CYD_model (
#   include_beta = 2,
#   include_eps = 1,
#   psi_J = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M6",
#   n_it = n_it,
#   adapt_delta = 0.99,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# 
# # M7 - M5 but no epsilon -------------------------------------------------------
# run_CYD_model (
#   include_eps = 0,
#   include_beta = 1,
#   psi_J = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M7",
#   n_it = n_it,
#    adapt_delta = 0.99,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# 
# # M8 - M6 but no epsilon -------------------------------------------------------
# run_CYD_model (
#   include_beta = 2,
#   include_eps = 0,
#   psi_J = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M8",
#   n_it = n_it,
#   adapt_delta = 0.9,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# 
# # M9 - M6 but serotype specific L ----------------------------------------------
# run_CYD_model (
#   L_K = 1, 
#   include_beta = 2,
#   include_eps = 1,
#   psi_J = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M9",
#   n_it = n_it,
#   adapt_delta = 0.9,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# M10 - M6 but serotype specific tau ----------------------------------------------
# run_CYD_model (
#   tau_K = 1, 
#   include_beta = 2,
#   include_eps = 1,
#   psi_J = 1,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M10",
#   n_it = n_it,
#   adapt_delta = 0.9,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# # M11 - M9 but drop epsilon  --------------------------------------------------
# run_CYD_model (
#   include_eps = 0,
#   L_K = 1, 
#   include_beta = 2,
#   psi_J = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M11",
#   n_it = n_it,
#   adapt_delta = 0.9,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
#
# # M12 - M9 but drop beta  -----------------------------------------------------
# run_CYD_model (
#   include_beta = 0,
#   L_K = 1, 
#   include_eps = 1,
#   psi_J = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M12",
#   n_it = n_it,
#   adapt_delta = 0.9,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# 
# # M13 - M9 but single beta  -----------------------------------------------------
# run_CYD_model (
#   include_beta = 1,
#   L_K = 1, 
#   include_eps = 1,
#   psi_J = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M13",
#   n_it = n_it,
#   adapt_delta = 0.9,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# 
# # M14 - M9 but drop psi_j  ----------------------------------------------
# run_CYD_model (
#   psi_J = 0,
#     L_K = 1, 
#   include_beta = 2,
#   include_eps = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   delta_KJ = 2, 
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M14",
#   n_it = n_it,
#   adapt_delta = 0.94,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
# 
# 
# # M15 - M9 but single delta ----------------------------------------------------
# run_CYD_model (
#   delta_KJ = 0, 
#   L_K = 1, 
#   include_beta = 2,
#   include_eps = 1,
#   psi_J = 1,
#   tau_K = 0,
#   mono_lc_SN = 1,
#   mono_lc_MU = 2,
#   include_pK2 = 0,
#   MU_test_SN = 1,
#   folder = "severe/M15",
#   n_it = n_it,
#   adapt_delta = 0.9,
#   stan_model = "CYD_model_severe.stan",
#   baseline_SP = read.csv("CYD/data/processed/baseline_SP.csv"),
#   cases =  readRDS("CYD/data/processed/cases_stan_format.RDS"),
#   mu =  read.csv("CYD/data/processed/mu.csv"),
#   severe = T
# )
