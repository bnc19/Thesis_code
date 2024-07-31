# Script to plot the main figure including model fit and the vaccine efficacy 
# estimate (Fig. 2) 

setwd("Chapter2")

# load functions 
library(tidyverse)
library(Hmisc)

# colours 
age_fill = scales::brewer_pal(palette = "Blues")(4)[2:4]
age_fill2 = scales::brewer_pal(palette = "Blues")(4)[c(2,4)]
serotype_fill = c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59")
trial_fill = c("#C51B8A", "#99CC99")

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.key.size = unit(0.5, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = unit(c(0,0,0,0), "cm")
    ) 
  
) 

# source files 
file.sources = paste0("TAK_VE/R/", list.files(path = "TAK_VE/R/"))
sapply(file.sources, source)
path = "TAK_VE/output/final/M32/"

# Data 
VE = readRDS(paste0(path, "VE.RDS"))
VE_model = extract_model_results(VE) 

# plot VE 
VE_all =  VE_model %>%
  filter(group == "VE") %>%
  separate(name, into = c("serostatus", "serotype", "age","outcome", "month")) %>%
  mutate(
    age = ifelse(age ==1, 1, 2), # age groups 2 and 3 have same VE 
    age = factor(age, labels = c("4-5yrs","6-16yrs")),
    outcome = factor(outcome, labels = c("symptomatic", "hospitalised")),
    serotype = factor(serotype, labels = c("DENV1", "DENV2", "DENV3", "DENV4")),
    month = as.numeric(month),
    serostatus = factor(serostatus, labels = c("seronegative", "monotypic", "multitypic"))) 


VE_symp = VE_all %>% 
  filter(outcome == "symptomatic") %>% 
  ggplot(aes(x = month , y = mean)) +
  geom_line(aes(color = age), linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age), alpha = 0.5) +
  labs(x = "Month", y = "Vaccine efficacy \nagainst symptomatic disease (%)") +
  scale_x_continuous(breaks = seq(0, 54,12)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
  facet_grid(serostatus ~ serotype, scales = "free_y" ) + 
  theme(legend.position = c(0.06,0.744)) + 
  scale_fill_manual(values = age_fill2) +
  scale_colour_manual(values = age_fill2) 

VE_hosp = VE_all %>% 
  filter(outcome != "symptomatic") %>% 
  ggplot(aes(x = month , y = mean)) +
  geom_line(aes(color = age), linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = age), alpha = 0.5) +
  labs(x = "Month", y = "Vaccine efficacy \nagainst hospitalisation (%)") +
  scale_x_continuous(breaks = seq(0, 54,12)) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", linewidth=1) +
  facet_grid(serostatus ~ serotype, scales = "free_y" ) + 
  theme(legend.position = c(0.07,0.76)) + 
  scale_fill_manual(values = age_fill2) +
  scale_colour_manual(values = age_fill2) 


ggsave(
  plot = VE_symp,
  filename =  "TAK_VE/output/figures/main_VE_symp.png",
  height = 10,
  width = 18,
  units = "cm",
  dpi = 600
)


ggsave(
  plot = VE_hosp,
  filename =  "TAK_VE/output/figures/main_VE_hosp.png",
  height = 10,
  width = 18,
  units = "cm",
  dpi = 600
)

