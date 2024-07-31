# Script to run sensitivity analysis on final survival model (M32) -------------

setwd("Chapter2")

# source functions
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)
n_it = 10000
# ------------------------------------------------------------------------------
# SA1 - M32 but HI = 1 
run_model (
  HI =1, 
  mono_lc_MU = 2,
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "final/SA1",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# SA2 - M32 but HI = 6
run_model (
  HI = 6, 
  mono_lc_MU = 2,
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "final/SA2",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# SA3 - M32 but MU can't test SN
run_model (
  MU_test_SN = 0,
  mono_lc_MU = 2,
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_SN = 1,
  tau_K = 1,
  folder = "final/SA3",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)

# SA4 - M32 but MO inf only 
run_model (
  MU_symp = 0, 
  mono_lc_MU = 2,
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "final/SA4",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)



# SA5 - M32 but HI - 24
run_model (
  HI = 24, 
  mono_lc_MU = 2,
  include_pK3 = 0,
  rho_K = 0,
  include_beta = 3,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "final/SA5",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
)
