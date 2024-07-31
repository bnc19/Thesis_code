# function to format stan data if simulated data -------------------------------
format_stan_data_sim = function(sim_data,
                                pop_formatted,
                                sim_param) {
  
  list2env(sim_param, env = environment())
  
  T = length(time)
  D = length(VCD_years)
  A = sum(VCD_years < 4)
  
# check that each time point has at least one case, or add one at random
# otherwise stan won't return AR as multinomial can't have no trials
add_one_case =  function(cases){ 
  
  if(any(colSums(cases) == 0)) {
    index = which(colSums(cases) == 0)
    cases[sample(1:nrow(cases),1),index] = 1 
  } 
  
  
  return(cases)
}

add_one_case_vec =  function(cases){ 
  
  if(any(cases == 0)) {
    index = which(cases == 0)
    cases[index] = 1 
  } 
  
  
  return(cases)
}

  VCD_D = sim_data$VCD_D$Y
  VCD_D = add_one_case_vec(VCD_D)
  
  HOSP_D = sim_data$HOSP_D$Y
  HOSP_D = add_one_case_vec(HOSP_D)
  
  # calculate # hosp by serostatus + trial arm + serotype, over time
  hosp_BVK4_m = array(sim_data$cases$HOSP_BVK4$Y, dim = c(B * K * V, 4))  # (D=4)
  hosp_BVK4_m = add_one_case(hosp_BVK4_m)

  # calculate # VCD by serostatus + trial arm + serotype, over time
  VCD_BVKD_m = array(sim_data$cases$VCD_BVKD$Y, dim = c(B * K * V, D))
  VCD_BVKD_m = add_one_case(VCD_BVKD_m)
  
  # calculate # hosp by serostatus + trial arm + age-group, over time
  hosp_BVJA_m = array(sim_data$cases$HOSP_BVJA$Y, dim = c(B * V * J, A))
  hosp_BVJA_m = add_one_case(hosp_BVJA_m)
  
  # calculate # VCD by serostatus + trial arm + age-group, over time
  VCD_BVJA_m = array(sim_data$cases$VCD_BVJA$Y, dim = c(B * V * J, A))
  VCD_BVJA_m = add_one_case(VCD_BVJA_m)
  
  # VCD age and serotype year and 2
  VCD_KJ2_m = array(sim_data$cases$VCD_KJ2$Y, dim = c(K * J, 2))
  VCD_KJ2_m = add_one_case(VCD_KJ2_m)
  
  # hosp age and serotype year and 2
  hosp_KJ2_m = array(sim_data$cases$HOSP_KJ2$Y, dim = c(K * J, 2))
  hosp_KJ2_m = add_one_case(hosp_KJ2_m)
  
  # VCD by serostatus, trial arm for t_5 (to censor population t_6)
  N_VCD_BV5_m  = array(sim_data$N_VCD_BV5$Y, dim = c(2, 2))
  
  # calculate # pop
  pop_BVJD = array(pop_formatted$N_pop_BVJD$N, dim = c(B, V, J, D))
  pop_BVD = array(pop_formatted$N_pop_BVD$N, dim = c(B, V, D))

  # data
  stan_data = list(
    time = time,
    T = T,
    J = J,
    K = K,
    B = B,
    V = V,
    D = D,
    A = A,
    C = C,
    HI = HI,
    R = R,
    include_pK3 = include_pK3,
    include_eps = include_eps,
    include_beta = include_beta,
    mono_lc_SN = mono_lc_SN,
    mono_lc_MU = mono_lc_MU,
    rho_K = rho_K,
    L_K = L_K,
    w_CK = w_CK,
    alpha_CK = alpha_CK,
    tau_K = tau_K,
    L_sd = L_sd,
    L_mean = L_mean,
    lower_bound_L = lower_bound_L,
    MU_test_SN = MU_test_SN,
    MU_symp = MU_symp,
    enhancement = enhancement,
    mu = mu,
    SP_J = sim_data$SP_J,
    VCD_D = VCD_D,
    HOSP_D = HOSP_D,
    pop_J = pop_formatted$N_pop_J,
    pop = pop_BVJD,
    pop_BVD = pop_BVD,
    VCD_BVKD = VCD_BVKD_m,
    VCD_BVJA = VCD_BVJA_m,
    HOSP_BVJA = hosp_BVJA_m,
    HOSP_BVK4 = hosp_BVK4_m,
    HOSP_KJ2 = hosp_KJ2_m,
    VCD_KJ2 = VCD_KJ2_m,
    N_VCD_BV5 = N_VCD_BV5_m
  )
  
  return(stan_data)
}

# function to summarise posterior estimates ------------------------------------

get_posterior = function(stan_fit,i,file_path) {
  
  # extract posterior estimates
  fit_ext = stan_fit$draws(format = "df")
  i1 = which(names(fit_ext) == "ll")
  posterior_chains = fit_ext[names(fit_ext)[1:i1]]
  
  posterior = summarise_draws(posterior_chains)
  
  posterior = posterior %>%  as.data.frame() %>% 
    select(variable, mean, q5, q95) %>%  
    mutate(across(mean:q95, ~round(.x,3)))
  
  write.csv(posterior, paste0(file_path,  "/posterior_", i, ".csv"))
  
  return(posterior)

}

# extract AR  ------------------------------------------------------------------

save_AR = function(stan_fit, i,file_path) {
  # extract posterior estimates
  fit_ext = stan_fit$draws(format = "df")
  AR = which(grepl("AR" , names(fit_ext)))
  AR_out = fit_ext[AR] %>%  as.data.frame()
  saveRDS(AR_out, paste0(file_path,  "/AR_", i, ".RDS"))
}

# function to get extract binomial confidence intervals ------------------------

bin_conf_AR = function(case_data) {
  case_data %>%
    mutate(
      mean =  binconf(Y, N, method = "exact")[, 1] * 100,
      lower = binconf(Y, N, method = "exact")[, 2] * 100,
      upper = binconf(Y, N, method = "exact")[, 3] * 100
    )  %>%
    select(-Y,-N) %>%
    mutate(type = "data")}

# function to calculate attack rates from formatted cases ----------------------

calculate_data_attack_rates = function(case_data, sim_cases) {
  
  # get true populations 
  pop = lapply(case_data, select, N)

  # select simulated cases to plot 
  sim_cases2 = sim_cases$cases
    
  # add outcome and true population to simulated cases 
  sim_cases_pop=list()
  for(i in 1:6) sim_cases_pop[[i]] = bind_cols(sim_cases2[[i]], pop[[i]])
  
  sim_attack_rates = lapply(sim_cases_pop, bin_conf_AR) %>%
    bind_rows() %>% # if missing then not dissagregated by that variables 
    mutate(age = factor(ifelse(is.na(age), "all", paste(age)), 
                        levels = c(levels(age), 'all'))) %>% 
    mutate(serostatus = factor(ifelse(is.na(serostatus), "all", paste(serostatus)), 
                               levels = c(levels(serostatus), 'all'))) %>% 
    mutate(serotype = factor(ifelse(is.na(serotype), "all", paste(serotype)), 
                             levels = c(levels(serotype), 'all'))) %>%  
    mutate(trial = factor(ifelse(is.na(trial), "all", paste(trial)), 
                          levels = c(levels(trial), 'all')))
  
}

# function to plot attack rates for simulated data -----------------------------

plot_sim_attack_rate = function(case_data, 
                                sim_cases,
                                file_path,
                                AR) {
  
# get model and simulated attack rates   
sim_attack_rates =  calculate_data_attack_rates(case_data, sim_cases)
model_attack_rates = extract_model_results(AR)

# plot age
  
  AR_BVJRD = sim_attack_rates %>%
    filter(age != "all", serotype == "all")
  
  AR_BVJRD_model =  model_attack_rates %>%
    filter(group == "AR_BVJRD") %>%
    separate(name, into = c("serostatus", "trial", "age", "outcome", "time")) %>%
    mutate(
      serostatus = factor(serostatus, labels = c("SN", "SP")),
      trial = factor(trial, labels = c("P", "V")),
      outcome = factor(outcome, labels = c("symp", "hosp")),
      age = factor( age, labels = c("4-5yrs", "6-11yrs", "12-16yrs")),
    ) %>%
    filter(time <= 4) %>% 
    mutate(time = ifelse(time == 1, 12,
                         ifelse(time == 2, 18,
                           ifelse(time == 3, 24, 36 ))))
  
  AR_BVJRD_plot = AR_BVJRD_model %>%  
    bind_rows(AR_BVJRD) %>%
    ggplot(aes(x = time, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = age,
        group = interaction(type, age)
      ),
      position = position_dodge(width = 2.5),
      size = 4
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, age),
        linetype = type,
        color = age
      ),
      position = position_dodge(width =  2.5),
      width =  0.4,
      linewidth = 1
    ) +
    facet_grid(outcome + serostatus ~ trial, scales = "free") +
    labs(x = "Month", y = "Attack rate (%)") +
    scale_x_continuous(breaks = c(12, 18, 24, 36))
  
# plot serotype
  
  AR_BVKRD = sim_attack_rates %>%
    filter(age == "all", serotype != "all") 
  
  AR_BVKRD_model =  model_attack_rates %>%
    filter(group == "AR_BVKRD") %>%
    separate(name,
             into = c("serostatus", "trial", "serotype", "outcome", "month")) %>%
    filter(outcome == 1) %>% # symp has each time point
    mutate(month = ifelse(month == 1, 12,
                          ifelse(month == 2, 18,
                           ifelse(month == 3, 24,
                            ifelse(month == 4, 36,
                             ifelse(month == 5, 48, 54)))))) %>%
    bind_rows(separate(  # add hosp which has fewer time points
      filter(model_attack_rates, group == "AR_BVKHD"),
      name, into = c("serostatus", "trial", "serotype", "time"))) %>%
    mutate(
      serostatus = factor(serostatus, labels = c("SN", "SP")),
      trial = factor(trial, labels = c("P", "V")),
      serotype = factor( serotype, labels = c("D1", "D2", "D3", "D4")),
    ) %>%
    mutate(month = ifelse(is.na(time), month,
                          ifelse(
                            time == 1, 24,
                            ifelse(time == 2, 36,
                                   ifelse(time == 3, 48, 54))
                          ))) %>%
    mutate(outcome = as.factor(ifelse(is.na(outcome), "hosp", "symp"))) %>%
    select(-time) %>% rename("time" = month) 
  
  
  AR_BVKRD_plot = AR_BVKRD_model %>%
    bind_rows(AR_BVKRD) %>%
    ggplot(aes(x = time, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = serotype,
        group = interaction(type, serotype)
      ),
      position = position_dodge(width = 3),
      size = 4
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, serotype),
        linetype = type,
        color = serotype
      ),
      position = position_dodge(width =  3),
      width =  0.4,
      linewidth = 1
    ) +
    facet_grid(outcome + serostatus ~ trial, scales = "free") +
    labs(x = "Month", y = "Attack rate (%)") +
    scale_x_continuous(breaks = c(12, 18, 24, 36, 48, 54))
  
# Plot by age and serotype
  
  AR_KJRD = sim_attack_rates %>%
    filter(age != "all", serotype != "all")
  
  AR_KJRD_model  =  model_attack_rates %>%
    filter(group == "AR_KJRD") %>%
    separate(name, into = c("serotype", "age", "outcome", "time")) %>%
    mutate(
      time = ifelse(time == 1, 12, 24),
      outcome = factor(outcome, labels = c("symp", "hosp")),
      serotype = factor(serotype, labels = c("D1", "D2", "D3", "D4")),
      age = factor(age, labels = c("4-5yrs", "6-11yrs", "12-16yrs"))
    )
    
  AR_KJRD_plot = AR_KJRD_model %>%
    bind_rows(AR_KJRD) %>%
    ggplot(aes(x = time, y = mean)) +
    geom_point(
      aes(
        shape = type,
        color = serotype,
        group = interaction(type, serotype)
      ),
      position = position_dodge(width = 2),
      size = 3
    ) +
    geom_errorbar(
      aes(
        ymin = lower ,
        ymax = upper ,
        group = interaction(type, serotype),
        linetype = type,
        color = serotype
      ),
      position = position_dodge(width =  2),
      width =  0.4,
      linewidth = 1
    ) +
    facet_grid(age ~ outcome, scales = "free") +
    labs(x = "Month", y = "Attack rate (%)") +
    scale_x_continuous(breaks = c(12, 24))
  
  # save
  
  ggsave(
    plot = AR_BVKRD_plot,
    filename =  paste0(file_path, "/AR_plot_BVKRD.jpg"),
    height = 40,
    width = 80,
    units = "cm",
    dpi = 600,
    scale = 0.6
  )
  
  ggsave(
    plot = AR_BVJRD_plot,
    filename =  paste0(file_path, "/AR_plot_BVJRD.jpg"),
    height = 40,
    width = 80,
    units = "cm",
    dpi = 600,
    scale = 0.6
  )
  
  
  ggsave(
    plot = AR_KJRD_plot,
    filename =  paste0(file_path, "/AR_plot_KJRD.jpg"),
    height = 40,
    width = 60,
    units = "cm",
    dpi = 600,
    scale = 0.6
  )
  
}

# Overall function to fit to simulated data in order to check
# the model is identifiable ----------------------------------------------------

fit_simulated_data = function(model,
                              list_data,
                              n_it = 4000,
                              n_chains = 4,
                              adapt_delta = 0.77,
                              i, 
                              init_values = function() {
                                list(
                                  hs = runif(2, 1, 10),
                                  hl = runif(1, 50, 200),
                                  ts = runif(2, -11, 11),
                                  lambda_D = replicate(6, runif(4, 0.00001, 0.01)),
                                  p = runif(3, 0.1, 0.9),
                                  pK3 = runif(4, 0.1, 0.9),
                                  gamma = runif(1, 0.1, 0.9),
                                  rho = runif(4, 0.1, 4),
                                  phi = runif(1, 0.1, 0.5),
                                  delta = runif(4, 0.1, 0.9),
                                  L = runif(4, 0, 4),
                                  w = runif(4, 1, 5),
                                  alpha = runif(4, 0, 2),
                                  beta = runif(2, 0, 2),
                                  tau = runif(4, 1, 4),
                                  lc = replicate(4, runif(3, 2, 10)),
                                  sens = runif(1, .8, 1),
                                  spec = runif(1, 0.95, 1),
                                  epsilon = runif(1, 1, 4)
                                )
                              },
                              file_path = "output/simulations/") {
    

comp_model = cmdstan_model(model, stanc_options = list("O1"))

file_path = paste0(file_path, "n_", i)
dir.create(file_path)

# fit to simulated data
stan_fit = comp_model$sample(
  data = list_data,
  chains = n_chains,
  parallel_chains = n_chains,
  iter_warmup = floor(n_it / 2),
  iter_sampling = floor(n_it / 2),
  thin = 1,
  init = init_values,
  seed = 14,
  refresh = 500,
  adapt_delta = adapt_delta
)

# save posterior chains and AR 
posterior = get_posterior(stan_fit,i,file_path)
save_AR(stan_fit, i,file_path)

}

# combine the simulated parameter with the estimated parameter 
combine_fit_sim = function(post, sim){
  
  select_param = sim[1:which(names(sim) == "omega")]
  select_param2 = select_param[c(
    "hs",
    "hl",
    "ts",
    "lambda_D",
    "p",
    "epsilon",
    "gamma",
    "rho",
    "phi",
    "delta",
    "L",
    "tau",
    "w",
    "lc",
    "alpha",
    "beta",
    "sens",
    "spec",
    "omega"
  )]
  
  params_df = select_param2 %>%
    unlist() %>%
    as.data.frame() %>%
    rownames_to_column("variables") %>%
    rename("fixed" = '.')
  
  both = post %>%
    filter(!variable %in% c("lp__", "ll",
                            "pK3[1]", "pK3[2]", "pK3[3]", "pK3[4]",
                            "kappa", 
                            "pSP[1]", "pSP[2]", "pSP[3]")) %>%
    bind_cols(params_df) %>%
    filter(
      !variable %in% c(
        "epsilon",
        "rho[2]",
        "rho[3]",
        "rho[4]",
        "L[2]",
        "L[3]",
        "L[4]",
        "w[2]",
        "w[3]",
        "w[4]",
        "lc[3,1]",
        "lc[1,2]",
        "lc[3,2]",
        "lc[1,3]",
        "lc[3,3]",
        "lc[1,4]",
        "lc[3,4]",
        "alpha[2]",
        "alpha[3]",
        "alpha[4]",
        "beta[2]"
      )
    )
  
  return(both)
}

# plot posterior individually 
plot_posterior = function(both,
                          file_path)  {
  
  p1 = both %>%
    filter(mean >= 1) %>%
    ggplot(aes(x = variable, y = mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = q5, ymax = q95)) +
    geom_point(aes(y = fixed), color = "red")
  
  
  p2 = both %>%
    filter(mean < 1 & !grepl("lambda", variables)) %>%
    ggplot(aes(x = variable, y = mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = q5, ymax = q95)) +
    geom_point(aes(y = fixed), color = "red")
  
  
  p3 = both %>%
    filter(grepl("lambda", variables)) %>%
    ggplot(aes(x = variable, y = mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = q5, ymax = q95)) +
    geom_point(aes(y = fixed), color = "red")
  
  grid = plot_grid(p1, p2, p3, ncol = 3)
  
 ggsave(
    plot = grid,
    paste0(file_path, "/posterior_plot.jpg"),
    height = 20,
    width = 60,
    unit = "cm"
  )
  
return(grid)
}