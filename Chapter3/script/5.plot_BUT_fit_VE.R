# Script to plot the main Butantan-DV figures including model fit and the 
# vaccine efficacy estimates

setwd("Chapter3")
dir.create("BUT/output/figures")
rm(list=ls())

# load functions 
library(tidyverse)
library(Hmisc)
library(wesanderson)
library(readxl)
library(cowplot)

# colours 
age_fill = scales::brewer_pal(palette = "Blues")(4)[c(2:4)]
serotype_fill = c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59")
trial_fill = c("#C51B8A", "#99CC99")

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)))


# source files 
file.sources = paste0("BUT/R/", list.files(path = "BUT/R/"))
sapply(file.sources, source)
path = "BUT/output/M7/"

# Data 
VE = readRDS(paste0(path, "VE.RDS"))
AR = readRDS(paste0(path, "AR.RDS"))
cases = readRDS("BUT/data/processed/cases_stan_format.RDS")

# plot attack rates 
# add aggregated populations to data and calculate attack rates

AR_age_data = calc_BUT_attack_rates(cases$Sy_BVJ)
AR_serotype_data = calc_BUT_attack_rates(cases$Sy_BVK)
AR_model = extract_BUT_model_results(AR)

# plot serotype serostatus attack rate 
AR_plot_BVK = AR_model %>%
  filter(group == "AR_BVK") %>%
  separate(name, into = c("serostatus", "arm", "serotype")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         serotype = factor(serotype, labels = c(paste0("DENV", 1:2)))) %>% 
  bind_rows(AR_serotype_data) %>%
  ggplot(aes(x = arm, y = mean)) +
  geom_point(aes(shape = type,color = serotype,
                 group = interaction(type, serotype)),
    position = position_dodge(width = 0.5),size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    group = interaction(type, serotype),
      linetype = type, color = serotype),
    position = position_dodge(width =  0.5),width =  0.4,linewidth = 1) +
  labs(x = "Trial arm", y = "Symptomatic \nattack rate (%)") +
  facet_wrap(~ serostatus) + 
  theme(legend.position =c(0.86,0.82)) +
  guides(color = guide_legend(ncol = 2),
         shape = guide_legend(ncol = 2)) +
  scale_color_manual(values = serotype_fill)
  

# plot symp attack rate by age and trial arm
AR_plot_BVJ = AR_model %>%
  filter(group == "AR_BVJ") %>%
  separate(name, into = c("serostatus", "arm", "age")) %>% 
  mutate(arm = factor(arm, labels = c("placebo", "vaccine")),
         serostatus = factor(serostatus, labels = c("seronegative", "seropositive")),
         age = factor(age, labels = c("2-6yrs", "7-17yrs", "18-59yrs"))) %>% 
  bind_rows(AR_age_data) %>%
  ggplot(aes(x = arm, y = mean)) +
  geom_point(aes(shape = type, color = age,
                 group = interaction(type, age)),
    position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper, 
                    group = interaction(type, age),
      linetype = type, color = age),
    position = position_dodge(width =  0.5), width =  0.4, linewidth = 1) +
  labs(x = "Trial arm", y = "Symptomatic \nattack rate (%)") + 
  guides(shape = "none",linetype = "none") +
  theme(legend.position =c(0.92,0.75)) +
  facet_wrap(~ serostatus)  + 
  scale_color_manual(values = age_fill)


g2 = plot_grid(AR_plot_BVK,  AR_plot_BVJ,
               labels = c("a", "b"), ncol = 1,
               rel_heights  = c(1.1,1))

ggsave(
  plot = g2,
  filename =  "BUT/output/figures/main_fit_B.png",
  height = 12,
  width = 18,
  units = "cm",
  dpi = 600
)



# plot VE 

VE_model = extract_BUT_model_results(VE)

VE_plot =  VE_model %>%
  separate(name, into = c("Serostatus", "Serotype", "Age", "Month")) %>%
  filter(Age == 1) %>% # all the same 
  mutate(
    Serotype = factor(Serotype, 
                      labels = c("DENV1", "DENV2")),
    Month = as.numeric(Month),
    Serostatus = factor(Serostatus, 
                        labels = c("seronegative", "monotypic", "multitypic"))) %>% 
  ggplot(aes(x = Month , y = mean)) +
  geom_line(aes(color = Serotype)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Serotype), alpha = 0.4) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  scale_x_continuous(breaks = seq(0, 24,6)) +
  facet_wrap(~Serostatus) +
  theme(legend.position = c(0.07,0.17)) +
  scale_y_continuous(limits = c(0,100)) +
  scale_color_manual(values = serotype_fill) +
  scale_fill_manual(values = serotype_fill)


ggsave(
  plot = VE_plot,
  filename =  "BUT/output/figures/ve_B.png",
  height = 6,
  width = 18,
  units = "cm",
  dpi = 600
)
