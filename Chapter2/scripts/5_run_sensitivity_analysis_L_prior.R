# Script to run sensitivity analysis on choice of L prior  ---------------------
setwd("Chapter2")

# Run final model (32) with SD of 0.25,0.5,0.75,1,1.25,1.5,1.75,2 
# A file will be saved with each BF in it 
 
# source functions
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)
n_it = 10000


# SD 0.25 
run_model (
  include_pK3 = 0,
  L_sd = 0.25, 
  rho_K = 0,
  include_beta = 3,
  mono_lc_MU = 2,
  mono_lc_SN = 1,
  tau_K = 1,
  MU_test_SN = 1,
  folder = "final/BF_0.25",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4)),
  BF = T, 
  plot = F
)

# SD 0.5 
run_model (
  L_sd = 0.5, 
  rho_K = 0,
  include_beta = 3,
 mono_lc_MU = 2,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 0,
  MU_test_SN = 1,
  folder = "final/BF_0.5",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4)),
  BF = T, 
 plot = F
)

# SD 0.75 
run_model (
  L_sd = 0.75, 
  rho_K = 0,
  include_beta = 3,
 mono_lc_MU = 2,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 0,
  MU_test_SN = 1,
  folder = "final/BF_0.75",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4)),
  BF = T, 
 plot = F
)


# SD 1 
run_model (
  L_sd = 1, 
  rho_K = 0,
  include_beta = 3,
 mono_lc_MU = 2,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 0,
  MU_test_SN = 1,
  folder = "final/BF_1",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4)),
  BF = T, 
 plot = F
)

# SD 1.25 
run_model (
  L_sd = 1.25, 
  rho_K = 0,
  include_beta = 3,
 mono_lc_MU = 2,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 0,
  MU_test_SN = 1,
  folder = "final/BF_1.25",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4)),
  BF = T, 
 plot = F
)

# SD 1.5 
run_model (
  L_sd = 1.5, 
  rho_K = 0,
  include_beta = 3,
 mono_lc_MU = 2,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 0,
  MU_test_SN = 1,
  folder = "final/BF_1.5",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4)),
  BF = T, 
 plot = F
)

# SD 1.75 
run_model (
  L_sd = 1.75, 
  rho_K = 0,
  include_beta = 3,
 mono_lc_MU = 2,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 0,
  MU_test_SN = 1,
  folder = "final/BF_1.75",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4)),
  BF = T, 
 plot = F
)

# SD 2 
run_model (
  L_sd = 2, 
  rho_K = 0,
  include_beta = 3,
 mono_lc_MU = 2,
  mono_lc_SN = 1,
  tau_K = 1,
  include_pK3 = 0,
  MU_test_SN = 1,
  folder = "final/BF_2",
  n_it = n_it,
  adapt_delta = 0.77,
  stan_model = "TAK_model.stan",
  baseline_SP = read.csv("data/processed/seropositive_by_age_baseline.csv"),
  VCD =  read.csv("data/processed/vcd_data.csv"),
  hosp = read.csv("data/processed/hosp_data.csv") ,
  mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4)),
  BF = T, 
 plot = F
)
