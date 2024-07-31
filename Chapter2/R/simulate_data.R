# function to generate parameter inputs from priors ----------------------------

simulate_parameters = function(pop_formatted,
                               mu,
                               VCD_years = c(12, 18, 24, 36, 48, 54) / 12,
                               time = 1:54,
                               B = 2,
                               K = 4,
                               V = 2,
                               R = 2,
                               J = 3,
                               C = 3,
                               A = 4,
                               HI = 12,
                               include_pK3 = 0,
                               include_eps = 0,
                               include_beta = 0,
                               mono_lc_SN = 0,
                               mono_lc_MU = 0,
                               rho_K = 0,
                               L_K = 0,
                               w_CK = 0,
                               alpha_CK = 0,
                               tau_K = 0,
                               L_sd = 1,
                               L_mean =0,
                               lower_bound_L =0,
                               MU_test_SN = 0,
                               MU_symp = 1,
                               enhancement = 1, 
                               lambda_prior = 4,
                               lc50_SN_mean = 4.5,
                               lc50_SP_mean = 6.5) {
# set up
steps = c(12, 6, 6, 12, 12, 6)
D = length(steps)
T = length(time)
  

# Empty arrays
hs = ts = p = pK3 = tau = w =  L = beta = lambda_mD =  lambda_m = pm =  c()
rho = alpha = delta = hK3 = hK1 = hK2 = q = shape1 = shape2 = c()
lc = array(NA, dim = c(C, K))
lambda_D = array(NA, dim = c(K, D))
  
# Draw from priors ------------------------------------------------------------
  
# titre priors 
hs[1] = rtruncnorm(1, a = 0, mean = 1.65, sd = 0.5)
hs[2] = rtruncnorm(1, a = 0, mean = 4.20, sd = 0.5)
hl    = rtruncnorm(1, a = 12, mean = 84, sd = 12)
ts[1] = rtruncnorm(1, a = -12, mean = -2.21, sd = 0.5)
ts[2] = rtruncnorm(1, a = -12, mean = 0.15, sd = 0.5)

# initial conditions 

# constrain pK3 by lambda
for (k in 1:K) for (d in 1:D)  lambda_D[k, d] = rlnorm(1, -7, lambda_prior)

for (d in 1:D) lambda_mD[d] = mean(lambda_D[, d])

# mean lambda across serotype and time
lambda_m =  weighted.mean(lambda_mD, steps)
pm = (1 - exp(-14 * 12 * lambda_m))  # mean p for a single serotype
if (pm < 1e-3) pm = 1e-3
if (pm > 0.99) pm = 0.99

# simplfied rearrangement of beta
shape1 = 49 * pm
shape2 = shape1 * (1 / pm - 1)

# constrain p1 and p2 for simulation so values are in line with pk3
p[1] = rbeta(1, shape1, shape2)
p[2] = rbeta(1, shape1, shape2)
p[3] = rbeta(1, shape1, shape2)

while (p[2] < p[3]) { # make sure p increases with age
  p[2] = rbeta(1, shape1, shape2)
}

while (p[2] < p[1]) { # make sure p increases with age
  p[1] = rbeta(1, shape1, shape2)
}

# draw prob parameters 
gamma = rtruncnorm(1, a = 0, b = 1, mean = 0.85, sd = 0.2)
for (k in 1:K) rho[k] = rtruncnorm(1, a = 0, mean = 1.97, sd = 0.40)
for (k in 1:K) delta[k] = rtruncnorm(1, a = 0, b = 1, mean = 0.25, sd = 0.1)
phi = rtruncnorm(1, a = 0, b = 1, mean = 0.25, sd = 0.05)
epsilon = rtruncnorm(1, a = 1, mean = 1, sd = 1)

# draw RR parameters 
for (k in 1:K) {
  L[k] = rtruncnorm(1, a = lower_bound_L, mean = L_mean,sd = L_sd)
  tau[k] = rtruncnorm(1, a = 0, mean = 1, sd = 1)
  w[k] = rtruncnorm(1, a = 0, mean = 1, sd = 2)
  lc[1, k] = rtruncnorm(1, a = 0, mean = lc50_SN_mean, sd = 1)
  lc[2, k] = rtruncnorm(1, a = 0, mean = lc50_SP_mean, sd = 1)
  lc[3, k] = rtruncnorm(1, a = 0, mean = lc50_SP_mean, sd = 1)
  alpha[k] = rtruncnorm(1, a = 0, mean = 0, sd = 2)
}

for (j in 1:J - 1) beta[j]  = rnorm(1, 0, 2)

# draw sens and spec 
sens = rtruncnorm(1, a = 0, b = 1, mean = 0.9, sd = 0.05)
spec = rtruncnorm(1, a = 0, b = 1, mean = 0.995, sd = 0.01)

# calculate # pop by trial arm, serostatus, age group, over time
N_pop_BVJD = pop_formatted$N_pop_BVJD %>%  
  arrange(year, age, trial , serostatus) 

pop_BVJD = array(N_pop_BVJD$N, dim = c(B, V, J, D))

return(
  list(
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
    L_mean= L_mean
  )
)
}

# function to extract simulated data from stan fit -----------------------------

extract_sim_cases = function(sim_fit) {

# get case data
fit_ext = sim_fit$draws(format = "df")
cases_ext = fit_ext[which(grepl("C_BVKJRD", names(fit_ext)))]
  
# format
cases = cases_ext %>%  
  as.data.frame() %>%
  mutate(ni = row_number()) %>%
  pivot_longer(cols = -ni)  %>%
  separate(name, into = c("group", "name"), sep = "\\[") %>%
  separate(name, into = c("name", NA), sep = "\\]") %>%
  separate(name, into = c("serostatus", "trial", 
                          "serotype","age", 
                          "outcome", "time")) %>% 
  mutate(
    serostatus = factor(serostatus, labels = c("SN", "SP")),
    trial = factor(trial, labels = c("P", "V")),
    serotype = factor(serotype, labels = c("D1", "D2", "D3", "D4")),
    age = factor(age, labels = c("4-5yrs", "6-11yrs", "12-16yrs")),
    outcome = factor(outcome, labels = c('symp', "hosp")),
    time = ifelse(time == 1, 12,
                  ifelse(time == 2, 18,
                         ifelse(
                           time == 3, 24,
                           ifelse(time == 4, 36,
                                  ifelse(time == 5, 48, 54)))))) %>%
  select(-c(ni, group)) %>%
  mutate(value = round(value, 0)) # integers only 
  
# seropositive
SP_J =  round(as.numeric(fit_ext[which(grepl("N_SP", names(fit_ext)))]), 0)
  
  return(list(cases, SP_J))
}

# Function to format cases after extracting them from the stan fixed model -----

format_sim_cases = function(sim_fit) {

cases = extract_sim_cases(sim_fit)[[1]]
SP = extract_sim_cases(sim_fit)[[2]]
  
# match to published trial data

VCD_BVKD = cases %>%  
  filter(outcome == "symp") %>%
  group_by(serostatus, trial, serotype, time, outcome) %>%
  summarise(Y = sum(value)) %>% 
  arrange(time, serostatus, trial , serotype)
  
VCD_BVJA = cases %>%  
  filter(outcome == "symp", time < 48) %>%
  group_by(serostatus, trial, age, time, outcome) %>%
  summarise(Y = sum(value)) %>%  
  arrange(time, serostatus, trial , age)

VCD_KJ2 = cases %>%  
  filter(outcome == "symp", time < 36) %>%
  mutate(time = ifelse(time == 12, 12, 24)) %>%
  group_by(serotype, age, time, outcome) %>% 
  summarise(Y = sum(value)) %>%  
  arrange(time, serotype, age)
  
HOSP_BVK4 = cases %>%  
  filter(outcome == "hosp") %>%
  mutate(time = ifelse(time == 12 | time == 18 | time == 24, 24, time)) %>%
  group_by(serostatus, trial, serotype, time, outcome) %>%
  summarise(Y = sum(value)) %>%  
  arrange(time, serostatus, trial , serotype)
  
HOSP_BVJA  = cases %>%  
  filter(outcome == "hosp", time < 48) %>%
  group_by(serostatus, trial, age, time, outcome) %>%
  summarise(Y = sum(value)) %>%  
  arrange(time, serostatus, trial , age)
  
HOSP_KJ2 = cases %>%  
  filter(outcome == "hosp", time < 36) %>%
  mutate(time = ifelse(time == 12, 12, 24)) %>%
  group_by(serotype, age, time, outcome) %>%
  summarise(Y = sum(value)) %>%  
  arrange(time, serotype, age)

VCD_D = cases %>% 
  filter(outcome == "symp") %>%
  group_by(time) %>%
  summarise(Y = sum(value)) %>% 
  arrange(time)

HOSP_D = cases %>% 
  filter(outcome == "hosp") %>%
  group_by(time) %>%
  summarise(Y = sum(value))%>% 
  arrange(time)
  
# cases in t_5 (not age-specific) for censoring
N_VCD_BV5  = cases %>%  
  filter(outcome == "symp", time == 48) %>%
  group_by(serostatus, trial) %>%
  summarise(Y = sum(value)) %>% 
  arrange(trial, serostatus)
  
# check distribution
plot_cases = cases %>%
  ggplot(aes(x = time, y = value)) +
  geom_bar(aes(fill = serotype),
           position = "stack", stat = "identity") +
  facet_grid(outcome + age ~ serostatus + trial, scales = "free_y") +
  theme(legend.position = "none")
 
out = list(
  cases = list(
  VCD_BVKD = VCD_BVKD,
  VCD_BVJA = VCD_BVJA,
  VCD_KJ2 = VCD_KJ2,
  HOSP_BVK4 = HOSP_BVK4,
  HOSP_BVJA = HOSP_BVJA,
  HOSP_KJ2 = HOSP_KJ2),
  VCD_D = VCD_D,
  HOSP_D = HOSP_D,
  N_VCD_BV5 = N_VCD_BV5,
  SP_J = SP,
  plot_cases = plot_cases
)

# if prior draws return NA for cases 
isNA = sapply(out, function(x) all(is.na(unlist(x))))

if (any(isNA)) {
  return(NA)
} else {
  return(out)
}

}

# overall function to simulate data --------------------------------------------
simulate_data = function (
    sim_parameters, sim_model){
  
# run fixed model 
sim_fit = sim_model$sample(
  data = sim_parameters,
  fixed_param = T,
  seed = 14,
  iter_sampling = 1,
  chains = 1
)
  
# format data to match published case data
sim = format_sim_cases(sim_fit)

}
