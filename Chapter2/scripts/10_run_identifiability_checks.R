# Simulated case data and fit the model 10 times to check for identifiability
# issues. Plot the resulting fits 
setwd("Chapter2")
# packages 
library(purrr)
library(truncnorm)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)  
library(cowplot)
library(bayesplot)
library(cmdstanr)
library(posterior)
library(Hmisc)
library(glue)
options(mc.cores = parallel::detectCores())
age_fill = scales::brewer_pal(palette = "Blues")(5)[2:5]

setwd(paste0(getwd(), "/TAK_VE"))

# source functions
file.sources = paste0("R/", list.files(path = "R/"))
sapply(file.sources, source)

# number of repeats
N_sim = 20

# compile fixed model 
comp_sim_model = cmdstan_model("models/fixed_model.stan", 
                               stanc_options = list("O1"))

# data
pop_formatted = readRDS("data/processed/pop_data.RDS")
mu =  array(read.csv("data/processed/n0_new.csv")$mean, dim = c(2, 4))
case_data = readRDS("data/processed/case_data.RDS")

# Draw from priors ------------------------------------------------------------
draw_priors = function(){
 
   # variables 
  VCD_years = c(12, 18, 24, 36, 48, 54) / 12
  time = 1:54
  B = V = R = 2
  K = A = 4
  J = C =3
  HI = 12
  steps = c(12, 6, 6, 12, 12, 6)
  D = length(steps)
  T = length(time)
  
  # Empty arrays
  hs = ts = p = pK3 = tau = w =  L = beta = lambda_mD =  lambda_m = pm =  c()
  rho = alpha = delta = hK3 = hK1 = hK2 = q = shape1 = shape2 = c()
  lc = array(NA, dim = c(C, K))
  lambda_D = array(NA, dim = c(K, D))
  
  # flags 
  include_pK3 = 0
  include_eps = 0
  include_beta = 3
  mono_lc_SN = 1
  mono_lc_MU = 2
  rho_K = 0
  L_K = 0
  w_CK = 0
  alpha_CK = 0
  tau_K = 1
  MU_test_SN = 1
  MU_symp = 1
  enhancement = 1
  L_sd = 1
  L_mean =0
  lower_bound_L = 0
  
   # titre priors 
  hs[1] = rtruncnorm(1, a = 0, mean = 2, sd = .5)
  hs[2] = rtruncnorm(1, a = 0, mean = 4, sd = .5)
  hl    = rtruncnorm(1, a = 12, mean = 84, sd = 20)
  ts[1] = rtruncnorm(1, a = -12, mean = -2, sd = .5)
  ts[2] = rtruncnorm(1, a = -12, mean = 0, sd = .5)
  
  # initial conditions 
  
  # constrain pK3 by lambda
  for (k in 1:K) for (d in 1:D)  lambda_D[k, d] = rtruncnorm(1, a =0, mean= 0.0005, sd= 0.005)
  for (d in 1:D) lambda_mD[d] = mean(lambda_D[, d])
  
  # mean lambda across serotype and time
  lambda_m =  weighted.mean(lambda_mD, steps)
  pm = (1 - exp(-14 * 12 * lambda_m))  # mean p for a single serotype
  if (pm < 1e-3) pm = 1e-3
  if (pm > 0.99) pm = 0.99
  
  # simplified rearrangement of beta
  vari = 0.02 * pm  * (1 - pm) 
  shape1 = (((1 - pm) / vari) - (1 / pm)) * (pm ^ 2)
  shape2 = shape1 * (1 / pm - 1)
  
  # constrain p1 and p2 for simulation so values are in line with pk3
  p[1] = rbeta(1, shape1, shape2)
  p[2] = rbeta(1, shape1, shape2)
  p[3] = rbeta(1, shape1, shape2)
  
  while (p[3] < p[2]) { # make sure p increases with age
    p[3] = rbeta(1, shape1, shape2)
  }
  
  while (p[2] < p[1]) { # make sure p increases with age
    p[1] = rbeta(1, shape1, shape2)
  }
  
  # draw prob parameters 
  gamma = rtruncnorm(1, a = 0, b = 1, mean = 0.4, sd = 0.15)
  for (k in 1:K) rho[k] = rtruncnorm(1, a = 0, mean = 2.5, sd = 1)
  for (k in 1:K) delta[k] = rtruncnorm(1, a = 0, b = 1, mean = 0.35, sd = 0.1)
  phi = rtruncnorm(1, a = 0, b = 1, mean = 0.2, sd = 0.1)
  epsilon = rtruncnorm(1, a = 1, mean = 2, sd = 2)
  
  # draw RR parameters 
  for (k in 1:K) {
    L[k] = rtruncnorm(1, a = lower_bound_L, mean = 0.2, sd = 1.5)
    tau[k] = rtruncnorm(1, a = 0, mean = 2, sd = 1)
    w[k] = rtruncnorm(1, a = 0, mean = 1, sd = 1)
    lc[1, k] = rtruncnorm(1, a = 0, mean = 5, sd = 1.5)
    lc[2, k] = rtruncnorm(1, a = 0, mean = 7, sd = 1.5)
    lc[3, k] = rtruncnorm(1, a = 0, mean = 7, sd = 1.5)
    alpha[k] = rtruncnorm(1, a = 0, mean = 0.25, sd = 0.5)
    omega = rtruncnorm(1, a = 0, mean = 1, sd = 0.5)
    kappa = rtruncnorm(1, a = 0, mean = 1, sd = 0.5)
  }
  
  for (j in 1:J - 1) beta[j]  = rnorm(1, 1, 2)
  
  # draw sens and spec 
  sens = rtruncnorm(1, a = 0, b = 1, mean = 0.9, sd = 0.05)
  spec = rtruncnorm(1, a = 0, b = 1, mean = 0.995, sd = 0.01)
  
  # calculate # pop by trial arm, serostatus, age group, over time
  N_pop_BVJD = pop_formatted$N_pop_BVJD %>%  
    arrange(year, age, trial , serostatus) 
  
  pop_BVJD = array(N_pop_BVJD$N, dim = c(B, V, J, D))
  
  return(list_data =  list(
    hs = hs,
    hl = hl,
    ts = ts,
    p = p, 
    lambda_D = lambda_D,
    gamma = gamma,
    rho = rho,
    delta = delta,
    phi = phi,
    L = L,
    tau = tau ,
    w = w,
    lc = lc,
    alpha = alpha,
    beta = beta,
    sens = sens,
    spec = spec,
    epsilon = epsilon,
    omega = omega,
    kappa = kappa,
    C = C,
    B = B,
    K = K,
    D = D,
    V = V,
    J = J,
    T = T,
    R = R,
    HI = HI,
    include_eps = include_eps,
    include_beta = include_beta,
    mono_lc_SN = mono_lc_SN,
    mono_lc_MU = mono_lc_MU,
    rho_K = rho_K,
    L_K = L_K,
    w_CK = w_CK,
    alpha_CK = alpha_CK,
    tau_K = tau_K,
    MU_test_SN = MU_test_SN,
    MU_symp = MU_symp,
    enhancement = enhancement,
    lower_bound_L = lower_bound_L,
    VCD_years = VCD_years,
    time = time,
    mu = mu,
    pop_J = pop_formatted$N_pop_J,
    pop = pop_BVJD,
    include_pK3 = include_pK3,
    L_sd = L_sd,
    L_mean = L_mean
  )
  )
}

set.seed(14)
param_sim = replicate(
  N_sim,
  draw_priors(),
  simplify = F
)

sim_cases = lapply(param_sim, simulate_data, comp_sim_model)

# plot simulated cases to check distributions are sort of reasonable 
sim_cases_to_plot = bind_rows(map(map(sim_cases,1), bind_rows), .id = "ni")

sim_cases_to_plot2 = sim_cases_to_plot %>% 
  ungroup() %>%  # convert to character to remove NAs
  # if missing then not disaggregated by that variable 
  mutate_if(is.factor, as.character) %>%
  mutate(trial = ifelse(is.na(trial), "both", trial),
         serostatus = ifelse(is.na(serostatus), "both", serostatus),
         serotype = ifelse(is.na(serotype), "all", serotype),
         age = ifelse(is.na(age), "all", age)) 


sim_cases_to_plot_f = factor_VCD(sim_cases_to_plot2)

sim =  bind_rows(sim_cases_to_plot_f) %>%   # combine simulated 
  mutate(ni = factor(ni, labels = c("true", paste("sim", 1:N_sim)), 
                     levels = c("true", 1:N_sim))) # plot true at the top 

symp_plot_sim_true = sim %>% 
  filter(outcome == "symptomatic") %>% 
  group_by(serostatus, trial, age, time, ni) %>%  # dont plot by serotype
  summarise(Y = sum(Y)) %>% 
  filter( !(age == "all" & time < 48)) %>%  # no age specific data after 48 months
  ggplot(aes(x = time, y = Y)) +
  geom_bar(aes(fill = age),
           position = "stack", stat = "identity") +
  facet_grid(serostatus + trial ~ ni, scales = "free_y") +
  theme(legend.position = "top") +
  scale_fill_manual(values = age_fill) +
  theme(strip.text.x = element_text(face = "bold"))  + 
  scale_x_continuous(limits = c(0,54), breaks = c(0, 24, 48)) +
  labs(y = "Simulated symptomatic cases", x = "Month")

hosp_plot_sim_true = sim %>% 
  filter(outcome != "symptomatic") %>% 
  group_by(serostatus, trial, age, time, ni) %>%  # dont plot by serotype
  summarise(Y = sum(Y)) %>% 
  filter( !(age == "all" & time < 48)) %>%  # no age specific data after 48 months
  ggplot(aes(x = time, y = Y)) +
  geom_bar(aes(fill = age),
           position = "stack", stat = "identity")+
  facet_grid(serostatus + trial ~ ni, scales = "free_y") +
  theme(legend.position = "top") +
  scale_fill_manual(values = age_fill) +
  theme(strip.text.x = element_text(face = "bold")) + 
  scale_x_continuous(breaks = c(0, 24, 48)) +
  labs(y = "Simulated hospitalised cases", x = "Month") 

# Save simulated case data 

ggsave(
  plot = symp_plot_sim_true, 
  "output/figures/symp_plot_sim.jpg",
  height = 30,
  width = 60,
  unit = "cm",
  scale = 0.6,
  dpi = 600
)

ggsave(
  plot = hosp_plot_sim_true, 
  "output/figures/hosp_plot_sim.jpg",
  height = 30,
  width = 60,
  unit = "cm",
  scale = 0.6,
  dpi = 600
)


# fit to the simulated data ----------------------------------------------------
# format simulated data to fit to
list_data = list()

for(i in 1:N_sim){
  list_data[[i]] = format_stan_data_sim(
    sim_param = param_sim[[i]], 
    sim_data  =  sim_cases[[i]],
    pop_formatted = pop_formatted
  )
}

for(i in 1:N_sim){
  fit_simulated_data(
    list_data = list_data[[i]],
    adapt_delta = 0.8,
    i = i, 
    model = "models/TAK_model.stan",
    n_it = 2000,
    file_path = "output/updated_simulations/")
}


for(i in 1:N_sim){  
  plot_sim_attack_rate(
    case_data = case_data,
    sim_cases = sim_cases[[i]],
    AR = readRDS(paste0("output/updated_simulations/n_", i, "/AR_", i, ".RDS")),
    file_path = paste0("output/updated_simulations/n_", i)
  )
} 


################################################################################
# CAN GO FROM HERE IF JUST PLOTTING IDENTIFIABILITY OF PARMETERS  
# Plus get param sim above
################################################################################

# plot overall figure of parameter reconstruction 
# collect all posterior 
index_files = which(grepl("n_", list.files(path = "output/updated_simulations/")))
post_sources = paste0("output/updated_simulations/",
                      list.files(path = "output/updated_simulations/")[index_files],
                      "/posterior.csv")

n1 = gsub("/posterior.csv", "", gsub("output/updated_simulations/n_", "", post_sources))
names(post_sources) = n1
post_sources2 = post_sources[as.character(sort(as.numeric(n1)))]


post = list()
for(i in 1:length(index_files)) post[[i]] = read.csv(paste0("output/updated_simulations/n_",i,"/posterior_", i, ".csv"))

combined = mapply(combine_fit_sim, post = post, sim = param_sim[c(1:length(index_files))], SIMPLIFY = F)

for(i in 1:length(index_files)) combined[[i]]$iteration = i

labels = c(
  "hs 1",
  "hs 2",
  "hl",
  "ts 1",
  "ts 2",
  "lambda[1][1]",
  "lambda[2][1]",
  "lambda[3][1]",
  "lambda[4][1]",
  "lambda[1][2]",
  "lambda[2][2]",
  "lambda[3][2]",
  "lambda[4][2]",
  "lambda[1][3]",
  "lambda[2][3]",
  "lambda[3][3]",
  "lambda[4][3]",
  "lambda[1][4]",
  "lambda[2][4]",
  "lambda[3][4]",
  "lambda[4][4]",
  "lambda[1][5]",
  "lambda[2][5]",
  "lambda[3][5]",
  "lambda[4][5]",
  "lambda[1][6]",
  "lambda[2][6]",
  "lambda[3][6]",
  "lambda[4][6]",
  "p[1]",
  "p[2]",
  "p[3]",
  "gamma",
  "rho",
  "phi",
  "delta[1]",
  "delta[2]",
  "delta[3]",
  "delta[4]",
  "L",
  "tau[1]",
  "tau[2]",
  "tau[3]",
  "tau[4]",
  "w",
  "n50[1]",
  "n50[2][1]",
  "n50[2][2]",
  "n50[2][3]",
  "n50[2][4]",
  'alpha',
  "beta",
  "sens",
  "spec",
  'omega'
)


levels = c(
  "hs[1]",
  "hs[2]",
  "hl",
  "ts[1]",
  "ts[2]",
  "lambda_D[1,1]",
  "lambda_D[2,1]",
  "lambda_D[3,1]",
  "lambda_D[4,1]",
  "lambda_D[1,2]",
  "lambda_D[2,2]",
  "lambda_D[3,2]",
  "lambda_D[4,2]",
  "lambda_D[1,3]",
  "lambda_D[2,3]",
  "lambda_D[3,3]",
  "lambda_D[4,3]",
  "lambda_D[1,4]",
  "lambda_D[2,4]",
  "lambda_D[3,4]",
  "lambda_D[4,4]",
  "lambda_D[1,5]",
  "lambda_D[2,5]",
  "lambda_D[3,5]",
  "lambda_D[4,5]",
  "lambda_D[1,6]",
  "lambda_D[2,6]",
  "lambda_D[3,6]",
  "lambda_D[4,6]",
  "p[1]",
  "p[2]",
  "p[3]",
  "gamma",
  "rho[1]",
  "phi",
  "delta[1]",
  "delta[2]",
  "delta[3]",
  "delta[4]",
  "L[1]",
  "tau[1]" ,
  "tau[2]" ,
  "tau[3]" ,
  "tau[4]" ,
  "w[1]",
  "lc[1,1]" ,
  "lc[2,1]" ,
  "lc[2,2]" ,
  "lc[2,3]" ,
  "lc[2,4]" ,
  "alpha[1]" ,
  "beta[1]" ,
  "sens",
  "spec",
  'omega'
)  
     
combined2 = combined %>%  
  bind_rows() %>% 
  mutate(variable = factor(variable, 
                           levels = levels, 
                           labels = labels))

p1 = combined2 %>% 
  filter(!grepl("lambda", variable)) %>%   # drop lambda 
  filter(!grepl("ts", variable)) %>% 
  filter(!grepl("hs", variable)) %>% 
  filter(variable != "hl") %>% 
  ggplot(aes(x = factor(iteration), y = fixed)) +
  geom_point(aes(color="true parameter"), size = 2) +
  geom_point(aes(y = mean, color = "posterior parameter"), size = 2) +
  geom_errorbar(aes(ymin =q5, ymax = q95), linewidth = 2, alpha = 0.3, width = 0) +
  facet_wrap(.~ variable, scales = "free", ncol = 4,
             labeller = label_parsed) +
  ylab(" ") + xlab("Simulation") +
  scale_color_manual(" ", values=c("grey", "red")) +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size=18),
        axis.text.x=element_blank())

p2 =combined2 %>% 
  filter(grepl("lambda", variable)) %>%   #  lambda 
  ggplot(aes(x = factor(iteration), y = fixed)) +
  geom_point(aes(color="true parameter"), size = 2) +
  geom_point(aes(y = mean, color = "posterior parameter"), size = 2) +
  geom_errorbar(aes(ymin =q5, ymax = q95), linewidth = 2, alpha = 0.3, width = 0) +
  facet_wrap(.~ variable, scales = "free", ncol = 4,
             labeller = label_parsed) +
  ylab(" ") + xlab("Simulation") +
  scale_color_manual(" ", values=c("grey", "red")) +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size=18),
        axis.text.x=element_blank())

ggsave(
  plot = p1, 
  "output/figures/simulation_plot.jpg",
  height = 27,
  width = 30,
  unit = "cm",
  dpi = 600
)
ggsave(
  plot = p2, 
  "output/figures/lambda_simulation_plot.jpg",
  height = 20,
  width = 30,
  unit = "cm",
  dpi = 600
)
