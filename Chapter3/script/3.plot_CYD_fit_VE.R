# Script to plot the main CYD figures including model fits and the vaccine  
# efficacy estimates (symp/hosp and symp/sev)

setwd("Chapter3")
rm(list = ls())
dir.create("CYD/output/figures")

library(tidyverse)
library(Hmisc)
library(wesanderson)
library(readxl)
library(cowplot)

# source files 
file.sources = paste0("CYD/R/", list.files(path = "CYD/R/"))
sapply(file.sources, source)
path = "CYD/output/M12/"
path_s = "CYD/output/severe/M9/" # (M10 in plots but started fitting from M0)

# Data 
VE = readRDS(paste0(path, "VE.RDS"))
AR = readRDS(paste0(path, "AR.RDS"))
VE_sev = readRDS(paste0(path_s, "VE.RDS"))
AR_sev = readRDS(paste0(path_s, "AR.RDS"))
VCD =  readRDS("CYD/data/processed/cases_stan_format.RDS")

# colours 
age_fill = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
serotype_fill = c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59", "#111111")
trial_fill = c("#C51B8A", "#99CC99")
cols = c("#67A9CF", "#C51B8A", "#99CC99", "#FFCC99")


theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)))

# plot attack rates 
AR_data = lapply(VCD, calc_CYD_attack_rates)
AR_model = extract_CYD_model_results(AR)  
AR_sev_model = extract_CYD_model_results(AR_sev)  

# plot attack rates ------------------------------------------------------------
# plot symp attack rate by serotype and trial arm
Sy_AR_plot_VK = AR_sev_model %>%
  filter(group == "V_AR_VK") %>%
  separate(name, into = c("arm", "serotype")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serotype = factor(serotype, labels = c("DENV1","DENV2", "DENV3", "DENV4"))) %>% 
  bind_rows(AR_data$Sy_VK) %>%
  ggplot(aes(x = arm, y = mean)) +
  geom_point(aes(shape = type, color = serotype,group = interaction(type, serotype)),
             position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, serotype),
                    linetype = type, color = serotype),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = " ", y = "Symptomatic \nattack rate (%)") +
  scale_color_manual(values = serotype_fill) +
  guides(shape = "none", linetype = "none")

# plot symp attack rate by age, serostatus and trial arm -----------------------
Sy_AR_plot_BVJ = AR_sev_model %>%
  filter(group == "V_AR_BVJ") %>%
  separate(name, into = c("serostatus", "arm", "age")) %>% 
  mutate(serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         arm = factor(arm, labels = c("placebo", "vaccine")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
  bind_rows(AR_data$Sy_BVJ) %>%
  unite(c(arm, serostatus), col = "x", sep = "\n") %>%
  ggplot(aes(x = x, y = mean)) +
  geom_point(aes(shape = type, color = age, group = interaction(type, age)),
             position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, age),
                    linetype = type, color = age),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = " ", y = "Symptomatic \nattack rate (%)") +
  scale_color_manual(values = age_fill) 


# plot hosp attack rate by age, serostatus, serotype and trial arm -------------
H_AR_plot_BVKJ = AR_sev_model %>%
  filter(group == "H_AR_BVKJ") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "age")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serotype = factor(serotype, labels = c("DENV1","DENV2", "DENV3", "DENV4")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>% 
  bind_rows(AR_data$Ho_BVKJ) %>%
  ggplot(aes(x = serostatus, y = mean)) +
  geom_point(aes(shape = type, color = serotype, group = interaction(type, serotype)),
             position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, serotype),
                    linetype = type, color = serotype),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = " ", y = "Hospitalisation \nattack rate (%)") +
  scale_color_manual(values = serotype_fill) +
  facet_grid(arm~ age) +
  guides(shape = "none", linetype = "none")


# plot hosp attack rate by age, serostatus, trial arm and time -----------------
H_AR_plot_BVJD =  AR_sev_model %>%
  filter(group == "H_AR_BVJD") %>%
  separate(name, into = c("serostatus", "arm", "age", "time")) %>% 
  mutate(serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         arm = factor(arm, labels = c("placebo", "vaccine")),
         age = factor(age, labels = c("2-8yrs", "9-16yrs")),
         time = factor(time, labels = c("1-13", "14-24", "25-36", "37-60"))) %>%
  bind_rows(AR_data$Ho_BVJD ) %>%
  ggplot(aes(x = time, y = mean)) +
  geom_point(aes(shape = type, color = age, group = interaction(type, age)),
             position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, age),
                    linetype = type, color = age),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  facet_grid(serostatus ~ arm) +
  labs(x = "Month", y = "Hospitalisation \nattack rate (%)") +
  scale_color_manual(values = age_fill) 

# plot severe attack rates -----------------------------------------------------

# get severe attack rate by age if seropositive and vaccine 
Se_AR_SP_V = AR_sev_model %>%
  filter(group == "S_AR_spvJ") %>% 
  separate(name, into = c("age")) %>% 
  mutate(age = factor(age, labels = c("2-8yrs", "9-16yrs")),
         arm = factor("vaccine"),
         serotype = factor("all"),
         serostatus = factor("seropositive")) %>%
  bind_rows(AR_data$Se_BVKJ) %>%  
  filter(arm == "vaccine", serostatus == "seropositive") 

# plot severe attack rate by serotype, trial arm, serotype and age
Se_AR_BVKJ = AR_sev_model %>%
  filter(group == "S_AR_BVKJ") %>%
  separate(name, into = c("serostatus", "arm", "serotype", "age")) %>% 
  mutate(
    arm = factor(arm, labels=c("placebo", "vaccine")),
    serotype = factor(serotype, labels =c("DENV1","DENV2","DENV3", "DENV4")),
    serostatus = factor(serostatus, labels = c("seronegative", "seropositive")), 
    age = factor(age, labels = c("2-8yrs", "9-16yrs"))) %>%  
  bind_rows(AR_data$Se_BVKJ) %>%  
  filter(!(arm == "vaccine" & serostatus == "seropositive")) %>%  
  bind_rows(Se_AR_SP_V) %>%  
  ggplot(aes(x = serostatus, y = mean)) +
  geom_point(aes(shape = type, color = serotype, group = interaction(type, serotype)),
             position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, serotype),
                    linetype = type, color = serotype),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = " ", y = "Severe \nattack rate (%)") +
  scale_color_manual(values = serotype_fill) +
  facet_grid(age ~ arm)  +
  guides(shape = "none", linetype = "none")

# plot severe attack rate by age, arm, time in SN   ----------------------------
Se_AR_VJD = AR_sev_model %>%
  filter(group == "S_AR_snVJD") %>%
  separate(name, into = c("arm", "age", "time")) %>%
  mutate(
    time = factor(time, labels = c("1-13", "14-24", "25-36", "37-60")),
    arm = factor(arm, labels = c("placebo", "vaccine")),
    age = factor(age, labels = c("2-8yrs", "9-16yrs"))
  ) %>%  
  bind_rows(AR_data$Se_VJD) %>% 
  ggplot(aes(x = time, y = mean)) +
  geom_point(aes(shape = type, color = age, group = interaction(type, age)),
             position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, group = interaction(type, age),
                    linetype = type, color = age),
                position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = "Month", y = "Severe \nattack rate (%)") +
  facet_wrap(~ arm) +
  scale_color_manual(values = age_fill) 

# save AR plot -----------------------------------------------------------------
S_AR_grid = cowplot::plot_grid(
  Sy_AR_plot_BVJ,
  Sy_AR_plot_VK, ncol =1,
  rel_heights = c(1,1.5),
  labels = c("a", "b"))


H_AR_grid = cowplot::plot_grid(H_AR_plot_BVJD,
                               H_AR_plot_BVKJ,
                               ncol = 1,
                               labels = c("a", "b"))

SE_AR_grid = cowplot::plot_grid(Se_AR_VJD,
                                Se_AR_BVKJ,
                                ncol = 1,
                                rel_heights = c(1,1.5),
                                labels = c("a", "b"))

ggsave(
  plot = S_AR_grid,
  filename =  "CYD/output/figures/S_AR_grid.png",
  height = 12,
  width = 18,
  units = "cm",
  dpi = 600
)

ggsave(
  plot = H_AR_grid,
  filename =  "CYD/output/figures/H_AR_grid.png",
  height = 16,
  width = 18,
  units = "cm",
  dpi = 600
)


ggsave(
  plot = SE_AR_grid,
  filename =  "CYD/output/figures/SE_AR_grid.png",
  height = 12,
  width = 18,
  units = "cm",
  dpi = 600
)

# plot symp hosp VE ----------------------------------------------------------------------
VE_model = extract_CYD_model_results(VE)

VE_f =  VE_model %>%
  filter(group == "VE_BKRT") %>% 
  separate(name, into = c("serostatus", "serotype", "outcome", "month")) %>%
  mutate(
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised")), 
    serotype = factor(serotype, 
                      labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    month = as.numeric(month),
    serostatus = factor(serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic"))) 
VE_plot = VE_f %>% 
ggplot(aes(x = month , y = mean)) +
  geom_line(aes(color = outcome)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = outcome), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 60,12)) +
  facet_grid(serostatus ~ serotype, scales = "free") +
  theme(legend.position = c(0.9,0.77),
        legend.title = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
      scale_color_manual(values =  cols) +
  scale_fill_manual(values =  cols) 

ggsave(
  plot = VE_plot,
  filename =  "CYD/output/figures/VE_plot_C.png",
  height = 12,
  width = 18,
  units = "cm",
  dpi = 600
)


# plot VE  ---------------------------------------------------------------------
VE_model_sev = extract_CYD_model_results(VE_sev)

VE_severe = VE_model_sev %>% 
  filter(group == "VE") %>% 
  separate(name, into = c("serostatus", "serotype", "age", "outcome", "month")) %>%
  mutate(
    age = factor(age, labels = c("2-8yrs", "9-16yrs")), 
    outcome = factor(outcome, labels = c("symptomatic", "severe")), 
    serotype = factor(serotype, 
                      labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    month = as.numeric(month),
    serostatus = factor(serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic")))

VE_severe_plot = VE_severe %>%  
  ggplot(aes(x = month , y = mean)) +
  geom_line(aes(color = outcome)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = outcome), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 60,12)) +
  facet_grid(serostatus ~ age ~ serotype, scales = "free") +
  theme(legend.position = c(0.88,0.9)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
  scale_color_manual(values =  cols) +
  scale_fill_manual(values =  cols) 

ggsave(
  plot = VE_severe_plot,
  filename =  "CYD/output/figures/severe_VE_plot_C.png",
  height = 17,
  width = 18,
  units = "cm",
  dpi = 600
)



# plot VE compared to SA on CYD cut-off period --------------------------------
cols = c("#67A9CF", "#C51B8A", "#99CC99", "#FFCC99")

# Data 
SA = readRDS("CYD/output/SA1/VE.RDS")
SA_model = extract_CYD_model_results(SA)

SA_f = SA_model %>%
  filter(group == "VE_BKRT") %>% 
  separate(name, into = c("serostatus", "serotype", "outcome", "month")) %>%
  mutate(
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised")), 
    serotype = factor(serotype, 
                      labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    month = as.numeric(month),
    serostatus = factor(serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  mutate(cutoff = "48 months")


VE_f = VE_f %>% 
  mutate(cutoff = "60 months")


out = VE_f %>% 
  bind_rows(SA_f) %>% 
  ggplot(aes(x = month, y  = mean)) + 
  geom_line(aes(color = cutoff)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = cutoff), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 60,12)) +
  facet_grid(outcome + serostatus ~ serotype, scales = "free") +
  theme(legend.position = c(0.9,0.89),
        legend.title = element_blank()) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
  scale_color_manual(values =  cols) +
  scale_fill_manual(values =  cols) 

ggsave(
  plot = out,
  filename =  "CYD/output/figures/SA_cutoff.png",
  height = 16,
  width = 18,
  units = "cm",
  dpi = 600
)


# plot VE compared to SA on symptomatic probabilities ---------------------------

