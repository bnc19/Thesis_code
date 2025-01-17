# Script to plot VE of Dengvaxia and Qdenga side by side

setwd("Chapter3")
rm(list = ls())
library(tidyverse)
source("compare_vaccines/R/plotting_functions.R")

# colours 
cols = c("#67A9CF", "#C51B8A", "#99CC99", "#FFCC99")

theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.spacing.y = unit(0, "pt"),
      legend.margin = margin(0, 0, 0, 0)))


# files 
TAK_file = "/Users/bethancracknelldaniels/Desktop/R projects/Thesis_code/Chapter2/output/M32/VE.RDS"
CYD_file = "CYD/output/M12/VE.RDS"

set.seed(123)

CYD_VE = readRDS(CYD_file)
TAK_VE = readRDS(TAK_file)

# VE 
C_VE_model = extract_model_results(CYD_VE)  
T_VE_model = extract_model_results(TAK_VE)  
C_VE_model$Vaccine = "Dengvaxia"
T_VE_model$Vaccine = "Qdenga"

# combine both vaccines
comb_VE = C_VE_model %>%
  bind_rows(T_VE_model) %>% 
  filter(group == "VE_BKRT") %>%  
  separate(name, into = c("serostatus", "serotype","outcome","time")) %>%  
  mutate(
    month = as.numeric(time),
    serostatus = factor(serostatus, levels = 1:3, labels = c("seronegative", "monotypic", "multitypic")),
    serotype = factor(serotype, levels = 1:4, labels = c(paste0("DENV", 1:4))),
    outcome = factor(outcome, levels = 1:2,  labels = c("symptomatic", "hospitalised"))
    )

plot_ve_over_time = comb_VE %>%   
  filter(month < 55) %>% 
  ggplot(aes(x = month, y = mean)) +
  geom_line(aes(color = Vaccine)) +
  geom_ribbon(aes( ymin = lower, ymax = upper, fill = Vaccine), alpha = 0.5) +
  labs(x = "Month", y = "Vaccine efficacy (%)") +
  geom_hline(yintercept = 0, linetype = "dashed",color = "black", linewidth=1) +
  facet_grid(serostatus + outcome ~ serotype, scales = "free") + 
  theme(legend.position = "top")+
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols)


ggsave(
  plot = plot_ve_over_time,
  filename = "compare_vaccines/output/VE_overtime.jpg",
  height = 18,
  width = 18,
  units = "cm",
  dpi = 300
)

