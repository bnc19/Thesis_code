

# Function to run the model equilibration
run_WHO_eq = function(
    i_country, # country index
    i_R0,  # transmission strength index 
    ns,    # number of posterior samples of VE estimates to use
    nparticles, # number of stochastic simulations per posterior sample
    posterior_file, # input file - vaccine posterior samples 
    nthreads,  # number of threads to run at once
    outf # output file 
           ) {
    gen_dengue = denguetak::model
    
    cn <= i_country
    vacc_cov =
      0.8 ## needed for config_params but not used for equilb runs
    vacc.age =
      9 ## needed for config_params but not used for equilb runs
    
    pars_nv =
      config_params(posterior_file, ns, vacc_cov, 1, vacc.age, 0, 0, 10)
    pars_v = pars_nv
    for (i in 1:ns)
      pars_v[[i]]$ZeroVE = 0
    
    
    ## equil runs
    
    for (j in 1:ns) {
      pars_nv[[j]]$R0_req = R0[i_R0, cn]
      pars_nv[[j]]$Mwt = R0[i_R0, cn]
    }
    
    mod = gen_dengue$new(
      pars = pars_nv,
      time = 0,
      n_particles = nparticles,
      n_threads = nthreads,
      seed = 1L,
      pars_multi = TRUE,
      deterministic = FALSE
    )
    
    # run to end of 2023
    mod.state = mod$run(run_steps_pre_vacc)
    
    # save final model conditions 
    save_loc = paste0(outf, "/", countries[cn])
    if (!dir.exists(save_loc))
      dir.create(save_loc)
    save_loc = paste0(save_loc, "/equil")
    if (!dir.exists(save_loc))
      dir.create(save_loc)
    saveRDS(mod.state, paste0(save_loc, "/equil_state_R0", i_R0, ".rds"))
    
  }


# Function to run the transmission model including vaccination

run_WHO_vacc =
  function(i_country, # country index
           i_R0,  # transmission strength index 
           ns, # number of posterior samples of VE estimates to use
           nparticles, # number of stochastic simulations per posterior sample
           posterior_file, # input file - vaccine posterior samples 
           nthreads,  # number of threads to run at once
           vacc.age, # age of vaccination
           vacc_cov, # vaccination coverage 
           DoVEInf, # include vaccine efficacy against infection (= 1)
           vacc_by_serostatus, # pre-vaccination screening (=1)
           mdy, # period of vaccination decay 
           outf  # output file 
           ) {
    
    
    gen_dengue = denguetak::model
    
    cn <= i_country
    
    # get model parameters including posterior samples 
    pars_v =
      config_params(posterior_file,
                    ns,
                    vacc_cov,
                    1,
                    vacc.age,
                    DoVEInf,
                    vacc_by_serostatus,
                    mdy)
    
    for (i in 1:ns)
      pars_v[[i]]$ZeroVE = 0
    
    # import initial conditions (having run model to equilibrium)
    mod.state =
      readRDS(paste0(outf, "/", countries[cn], "/equil/equil_state_R0", i_R0, ".rds"))
    
    ### Real VE runs
    
    for (j in 1:ns) {
      pars_v[[j]]$R0_req = R0[i_R0, cn]
      pars_v[[j]]$Mwt = R0[i_R0, cn]
    }
    
    # set up model 
    mod = gen_dengue$new(
      pars = pars_v,
      time = 0,
      n_particles = nparticles,
      n_threads = 40L,
      seed = 1L,
      pars_multi = TRUE,
      deterministic = FALSE
    )
    
    # set initial conditions to end of 2023 from results of model equilibration 
    mod$update_state(pars = pars_v,
                     state = mod.state,
                     time = run_steps_pre_vacc)
    mod.var.list = mod$info()[[1]]
    out_vars =
      mod.var.list$index[grepl("out", names(mod.var.list$index))]
    N_out = length(out_vars)
    
    # run to end of 2100
    mod.state.v = mod$run(run_steps_post_vacc)
    result.v = out_vars
    for (j in 1:N_out)
      result.v[[j]] = mod.state.v[out_vars[[j]], , ]
    
    save_loc = paste0(outf, "/", countries[cn])
    if (!dir.exists(save_loc))
      dir.create(save_loc)
    save_loc = paste0(save_loc, "/vacc")
    if (!dir.exists(save_loc))
      dir.create(save_loc)
    save_loc = paste0(save_loc, "/vc_", vacc_cov)
    if (!dir.exists(save_loc))
      dir.create(save_loc)
    save_loc =
      paste0(save_loc, "/", ifelse(DoVEInf == 1, "VI", "VS"), "_D", mdy)
    if (!dir.exists(save_loc))
      dir.create(save_loc)
    saveRDS(result.v,
            paste0(save_loc, "/vacc_res_", vacc.age, "_R0", i_R0, ".rds"))
    
  }

# Function to run the transmission model without vaccination

run_WHO_novacc =
  function(i_country, # country index
           i_R0,  # transmission strength index 
           ns, # number of posterior samples of VE estimates to use
           nparticles, # number of stochastic simulations per posterior sample
           posterior_file, # input file - vaccine posterior samples 
           nthreads,  # number of threads to run at once
           vacc.age, # age of vaccination
           vacc_cov, # vaccination coverage 
           vacc_by_serostatus, # pre-vaccination screening (=1)
           outf  # output file 
           ) {
    
    gen_dengue = denguetak::model
    
    cn <= i_country
    
    # get model parameters including posterior samples 
    pars_nv =
      config_params(posterior_file,
                    ns,
                    vacc_cov,
                    1,
                    vacc.age,
                    1,
                    vacc_by_serostatus,
                    10)
    
    # import initial conditions (having run model to equilibrium)
    mod.state =
      readRDS(paste0(outf, "/", countries[cn], "/equil/equil_state_R0", i_R0, ".rds"))
    
    ## Zero VE vacc runs
    for (j in 1:ns) {
      pars_nv[[j]]$R0_req = R0[i_R0, cn]
      pars_nv[[j]]$Mwt = R0[i_R0, cn]
    }
    
    # set up model 
    mod = gen_dengue$new(
      pars = pars_nv,
      time = 0,
      n_particles = nparticles,
      n_threads = nthreads,
      seed = 1L,
      pars_multi = TRUE,
      deterministic = FALSE
    )
    
    # set initial conditions to end of 2023 from results of model equilibration 
    mod$update_state(pars = pars_nv,
                     state = mod.state,
                     time = run_steps_pre_vacc)
    
    
    mod.var.list = mod$info()[[1]]
   
    out_vars =
      mod.var.list$index[grepl("out", names(mod.var.list$index))]
    
     N_out = length(out_vars)
    
    # run to end of 2100
    mod.state.nv = mod$run(run_steps_post_vacc)
    
    result.nv = out_vars
    for (j in 1:N_out)
      result.nv[[j]] = mod.state.nv[out_vars[[j]], , ]
    
    save_loc = paste0(outf, "/", countries[cn])
    if (!dir.exists(save_loc))
      dir.create(save_loc)
    save_loc = paste0(save_loc, "/no_vacc")
    if (!dir.exists(save_loc))
      dir.create(save_loc)
    save_loc = paste0(save_loc, "/vc_", vacc_cov)
    if (!dir.exists(save_loc))
      dir.create(save_loc)
    saveRDS(result.nv,
            paste0(save_loc, "/no_vacc_res_", vacc.age, "_R0", i_R0, ".rds"))
  }
