

# Function to save posterior summary of parameters -----------------------------
save_post_summary = function(stan_fit,file_path){
  posterior = summary(stan_fit, probs = c(0.025, 0.975))$summary
  write.csv(posterior, paste0(file_path, "/posterior.csv"))
}


# Function to calculate posterior mean and 95% CrI of parameter values ---------

summarise_posterior_chains = function(stan_fit)  {
  
  posterior_chains= extract_post_chains(stan_fit)
  
  summary = data.frame(
    mean  = apply(posterior_chains, 2,  mean),
    lower = apply(posterior_chains, 2, quantile, probs = 0.025),
    upper = apply(posterior_chains, 2, quantile, probs = 0.975)
  )
  return(summary)
  
}

save_posterior = function(file_path){

  library(dplyr)
  library(tidyr)
  
  
  VE_dis = readRDS(paste0(file_path, "fit_ext.RDS"))
  
  lambda = VE_dis$lambda_D %>%  as.data.frame() %>% 
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value) %>%  
    summarise(
      lower = round(quantile(value, 0.025), 3),
      mean = round(mean(value),3),
      upper = round(quantile(value, 0.975), 3)
    ) %>% 
    mutate(param = "lambda")
  
  p = VE_dis$p %>%  as.data.frame() %>% 
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value) %>%  
    summarise(
      lower = round(quantile(value, 0.025), 2),
      mean = round(mean(value),2),
      upper = round(quantile(value, 0.975), 2)
    ) %>% 
    mutate(param = "historic exposure")
  
  h  = VE_dis$hs %>%  as.data.frame() %>% 
    bind_cols(VE_dis$ts %>%  as.data.frame()) %>%
    bind_cols(VE_dis$hl %>%  as.data.frame()) %>% 
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value) %>%  
    summarise(
      lower = round(quantile(value, 0.025), 2),
      mean = round(mean(value),2),
      upper = round(quantile(value, 0.975), 2)
    ) %>% 
    mutate(param = c("ts", "hl", rep("hs", 2)))
  
  
  rho = VE_dis$rho %>%  as.data.frame() %>% 
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value) %>%  
    summarise(
      lower = round(quantile(value, 0.025), 2),
      mean = round(mean(value),2),
      upper = round(quantile(value, 0.975),2)
    ) %>% 
    mutate(param = "rho") 

  
  delta = VE_dis$delta %>%  as.data.frame() %>% 
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value) %>%  
    summarise(
      lower = round(quantile(value, 0.025), 2),
      mean = round(mean(value),2),
      upper = round(quantile(value, 0.975), 2)
    ) %>% 
    mutate(param = "delta") 

  gamma = VE_dis$gamma %>%  as.data.frame() %>% 
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value) %>%  
    summarise(
      lower = round(quantile(value, 0.025), 2),
      mean = round(mean(value),2),
      upper = round(quantile(value, 0.975), 2)
    ) %>% 
    mutate(param = "gamma") 
  
  L = VE_dis$L %>%  as.data.frame() %>% 
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value) %>%  
    summarise(
      lower = round(quantile(value, 0.025), 2),
      mean = round(mean(value),2),
      upper = round(quantile(value, 0.975), 2)
    ) %>% 
    mutate(param = "L") 
  
  w = VE_dis$w %>%  as.data.frame() %>% 
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value) %>%  
    summarise(
      lower = round(quantile(value, 0.025), 2),
      mean = round(mean(value),2),
      upper = round(quantile(value, 0.975), 2)
    ) %>% 
    mutate(param = "w")
  
  epsilon = VE_dis$epsilon %>%  as.data.frame() %>% 
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value) %>%  
    summarise(
      lower = round(quantile(value, 0.025), 2),
      mean = round(mean(value),2),
      upper = round(quantile(value, 0.975), 2)
    ) %>% 
    mutate(param = "epsilon")
  
  Z = VE_dis$Z %>%  as.data.frame() %>% 
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    group_by(name) %>% 
    mutate(value = value) %>%  
    summarise(
      lower = round(quantile(value, 0.025), 2),
      mean = round(mean(value),2),
      upper = round(quantile(value, 0.975), 2)
    ) %>% 
    mutate(param = "Z")
  
  
  out = bind_rows(lambda, p,h, L, w, rho, delta, gamma, Z, epsilon)
  
  out = out %>% 
    unite("Posterior", lower, upper, sep = "-") %>%
    unite("Posterior", mean,Posterior, sep = " (") %>% 
    mutate(Posterior = paste0(Posterior, ")"))
    
  
  write.csv(out, paste0(file_path, "post.csv") )
  
}


