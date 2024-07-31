# Run prior predictive checks
setwd("Chapter2")

# Simulate case data from the model priors 

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
options(mc.cores = parallel::detectCores())
age_fill = scales::brewer_pal(palette = "Blues")(5)[2:5]

# source functions
file.sources = paste0("TAK_VE/R/", list.files(path = "TAK_VE/R/"))
sapply(file.sources, source)

# data 
case_data = readRDS("TAK_VE/data/processed/case_data.RDS")

# parameters 
N_sim = 10 
set.seed(14)


# simulate parameters from exact priors  
PPC_param_sim = replicate(
  N_sim,
  simulate_parameters(
    pop_formatted = readRDS("TAK_VE/data/processed/pop_data.RDS"),
    mu =  array(read.csv("TAK_VE/data/processed/n0_new.csv")$mean, dim = c(2, 4)),
    include_pK3 = 0,
    rho_K = 0,
    include_beta = 3,
    mono_lc_MU = 1,
    mono_lc_SN = 1,
    tau_K = 1,
    MU_test_SN = 1,
    lambda_prior = 2,
    lc50_SN_mean = 4.5,
    lc50_SP_mean = 6.5
  ),
  simplify = F
)

# simulate case data from parameters  
comp_sim_model = cmdstan_model("TAK_VE/models/fixed_model.stan", 
                               stanc_options = list("O1"))

PPC_sim_cases = lapply(PPC_param_sim, simulate_data, comp_sim_model)

# return selection of plots 
sim_cases_to_plot = bind_rows(map(map(PPC_sim_cases,1), bind_rows), .id = "ni")

sim_cases_to_plot2 = sim_cases_to_plot %>% 
  ungroup() %>%  # convert to character to remove NAs
  # if missing then not disaggregated by that variable 
  mutate_if(is.factor, as.character) %>%
  mutate(trial = ifelse(is.na(trial), "both", trial),
         serostatus = ifelse(is.na(serostatus), "both", serostatus),
         serotype = ifelse(is.na(serotype), "all", serotype),
         age = ifelse(is.na(age), "all", age)) 
  

sim_cases_to_plot_f = factor_VCD(sim_cases_to_plot2)

sim_and_true = case_data %>%
  bind_rows() %>%
  rename(time = month) %>%  # named it wrong when simulating data
  mutate(ni = "true") %>% # true case data
  bind_rows(sim_cases_to_plot_f) %>%   # add simulated 
  mutate(ni = factor(ni, labels = c("true", paste("sim", 1:10)), 
                     levels = c("true", 1:10))) # plot true at the top 
  
symp_plot_sim_true = sim_and_true %>% 
  filter(outcome == "symptomatic",
         serostatus != "both",
         trial != "both") %>% 
  group_by(serostatus, trial, age, time, ni) %>%  # dont plot by serotype
  summarise(Y = sum(Y)) %>% 
  filter(age != "all") %>% 
  ggplot(aes(x = time, y = Y)) +
  geom_bar(aes(fill = age),
           position = "stack", stat = "identity") +
  facet_grid(serostatus + trial ~ ni, scales = "free_y") +
  theme(legend.position = "top") +
  scale_fill_manual(values = age_fill) +
  theme(strip.text.x = element_text(face = "bold"),
        legend.title = element_blank())+
  scale_x_continuous(breaks = c(12,24,36)) +
  xlab("Month") + ylab("Cases")

hosp_plot_sim_true = sim_and_true %>% 
  filter(outcome != "symptomatic",
         serostatus != "both",
         trial != "both") %>% 
  group_by(serostatus, trial, age, time, ni) %>%  # dont plot by serotype
  summarise(Y = sum(Y)) %>% 
  filter(age != "all") %>% 
  ggplot(aes(x = time, y = Y)) +
  geom_bar(aes(fill = age),
           position = "stack", stat = "identity")+
  facet_grid(serostatus + trial ~ ni, scales = "free_y") +
  theme(legend.position = "top") +
  scale_fill_manual(values = age_fill) +
  theme(strip.text.x = element_text(face = "bold"),
        legend.title = element_blank())+
  scale_x_continuous(breaks = c(12,24,36))+
  xlab("Month") + ylab("Cases")



ggsave(plot = symp_plot_sim_true, "TAK_VE/output/figures/prior_checks_symp.jpg", 
       height = 12,  width = 18, unit = "cm" )

ggsave(plot = hosp_plot_sim_true, "TAK_VE/output/figures/prior_checks_hosp.jpg",
       height = 12, width = 18, unit = "cm" )
