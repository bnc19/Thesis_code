# Plot model calibration 
setwd("Chapter2")

# set up 
library(tidyverse)
library(Hmisc)
library(wesanderson)

age_fill = scales::brewer_pal(palette = "Blues")(4)[2:4]
serotype_fill =  c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59")
trial_fill = c("#C51B8A", "#99CC99")


file.sources = paste0("TAK_VE/R/", list.files(path = "TAK_VE/R/"))
sapply(file.sources, source)
path = "TAK_VE/output/final/M32/"


theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = unit(c(0,0,0,0), "cm")
    ) 
  
) 

# symptomatic data 
AR = readRDS(paste0(path, "AR.RDS"))
cases =  readRDS("TAK_VE/data/processed/case_data.RDS")


# Calculate model attack rates 
AR_model = extract_model_results(AR)  

# Calculate data attack rates 
calc_AR = function(data){
  data %>%  
  mutate(mean =  binconf(Y,N, method = "exact")[,1] * 100,  
         lower = binconf(Y,N, method = "exact")[,2] * 100, 
         upper = binconf(Y,N, method = "exact")[,3] * 100) %>%  
  mutate(type = "data")
}

cases_AR = lapply(cases, calc_AR)

# Plot attack rates by age
AR_plot_BVJRD =  AR_model %>%
  filter(group == "AR_BVJRD") %>%
  separate(name, into = c("serostatus", "trial", "age", "outcome", "time")) %>% 
  mutate(serostatus = ifelse(serostatus== 1, "seronegative", "seropositive"),
         trial = ifelse(trial == 1, "placebo", "vaccine"),
         outcome = factor(outcome, labels = c("symptomatic", "hospitalised")),
         age = factor(age,
                      labels = c("4-5yrs", "6-11yrs", "12-16yrs"),
                      levels = c(1,2,3)),) %>% 
  mutate(month = ifelse(time == 1, 12, ifelse(time == 2, 18, 
                                              ifelse(time == 3,24,
                                                     ifelse(time == 4, 36,
                                                            ifelse( time == 5, 48, 54)))))) %>%
  filter(month <= 36) %>% 
  bind_rows(cases_AR$VCD_BVJA, cases_AR$HOSP_BVJA) %>%
  ggplot(aes(x = month, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = age,
      group = interaction(type, age)
    ),
    position = position_dodge(width = 3),
    size = 2
  ) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, age),
      linetype = type,
      color = age
    ),
    position = position_dodge(width =  3),
    width =  0.4,
    linewidth = 1
  ) +
  facet_grid(outcome+serostatus ~ trial, scales = "free") +
  labs(x = "Month", y = "Attack rate (%)") +
  theme(legend.position = "top") +
  scale_x_continuous(breaks = c(12,18,24,36)) +
  scale_colour_manual(values = age_fill) 

# Plot serotype attack rates 

  AR_plot_BVKRD =  AR_model %>%
  filter(group == "AR_BVKRD") %>%
  separate(name, into = c("serostatus", "trial", "serotype", "outcome", "month")) %>% 
  filter(outcome == 1) %>% # symphas each time point 
  mutate(month = ifelse(month == 1, 12, 
                        ifelse(month == 2, 18, 
                               ifelse(month == 3,24,
                                      ifelse(month == 4, 36, 
                                             ifelse( month == 5, 48, 54)))))) %>%
  bind_rows(separate(filter(AR_model, group == "AR_BVKHD"), # add hosp which has fewer time points 
                     name, into = c("serostatus", "trial", "serotype", "time"))) %>%  
  mutate(serostatus = ifelse(serostatus== 1, "seronegative", "seropositive"),
         trial = ifelse(trial == 1, "placebo", "vaccine"),
         outcome = ifelse(outcome == 1, "symptomatic", "hospitalised"),
         serotype = factor(serotype, labels= c("DENV1", "DENV2", "DENV3", "DENV4"))) %>% 
  mutate(month = ifelse(is.na(time), month,
                        ifelse(time == 1, 24,
                        ifelse(time == 2, 36,
                               ifelse(time == 3, 48, 54))))) %>%
  mutate(outcome = ifelse(is.na(outcome), "hospitalised", outcome)) %>%  
  bind_rows(cases_AR$N_hosp_BVK4, cases_AR$VCD_BVKD) %>%
  mutate(outcome = factor(outcome, levels = c("symptomatic", "hospitalised"))) %>% 
  ggplot(aes(x = month, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = serotype,
      group = interaction(type, serotype)
    ),
    position = position_dodge(width =6),
    size = 2
  ) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, serotype),
      linetype = type,
      color = serotype
    ),
    position = position_dodge(width =  6),
    width =  0.4,
    linewidth = 1
  ) +
  facet_grid(outcome+serostatus ~ trial, scales = "free" ) +
  labs(x = "Month", y = "Attack rate (%)") +
  theme(legend.position = "top")+
  scale_x_continuous(breaks = c(12,18,24,36, 48,54)) +
  scale_colour_manual(values = serotype_fill) 

# Plot by age and serotype ---- 

AR_plot_KJRD =  AR_model %>%
  filter(group == "AR_KJRD") %>%
  separate(name, into = c("serotype", "age", "outcome", "month")) %>% 
  mutate(month = ifelse(month == 1, 12, 24),
         outcome = factor(outcome, labels = c("symptomatic", "hospitalised")),
         serotype = factor(serotype, labels= c("DENV1", "DENV2", "DENV3", "DENV4")),
         age = factor(age, labels = c("4-5yrs", "6-11yrs", "12-16yrs"))) %>% 
  bind_rows(cases_AR$VCD_KJ2, cases_AR$HOSP_KJ2) %>%
  ggplot(aes(x = month, y = mean)) +
  geom_point(
    aes(
      shape = type,
      color = serotype,
      group = interaction(type, serotype)
    ),
    position = position_dodge(width = 4),
    size = 2
  ) +
  geom_errorbar(
    aes(
      ymin = lower ,
      ymax = upper ,
      group = interaction(type, serotype),
      linetype = type,
      color = serotype
    ),
    position = position_dodge(width =  4),
    width =  0.4,
    linewidth = 1
  ) +
  facet_grid(outcome~age, scales="free") +
  labs(x = "Month", y = "Attack rate (%)") +
  theme(legend.position = "top")+
  scale_x_continuous(breaks = c(12,24)) +
  scale_colour_manual(values = serotype_fill) 


# output 
ggsave(
  plot = AR_plot_BVJRD,
  filename =  "TAK_VE/output/figures/AR_plot_BVJRD.png",
  height = 12,
  width = 18,
  units = "cm",
  dpi = 600,
)

ggsave(
  plot = AR_plot_BVKRD,
  filename =  "TAK_VE/output/figures/AR_plot_BVKRD.png",
  height = 12,
  width = 18,
  units = "cm",
  dpi = 600,
)


ggsave(
  plot = AR_plot_KJRD,
  filename =  "TAK_VE/output/figures/AR_plot_KJRD.png",
  height = 8,
  width = 18,
  units = "cm",
  dpi = 600
)


