# Plot main sensitivity analysis on HI, mu test SP, mu symp

setwd("Chapter2")

library(tidyverse)
library(wesanderson)
fill = c("#000000", "#67A9CF", "#C51B8A", "#99CC99",  "#FFCC99", "#D55E00")


# collect all SA
index_files = which(grepl("SA", list.files(path = "output/final/")))
SA_sources = paste0("output/final/", list.files(path = "output/final/")[index_files], "/VE.RDS")
VE = lapply(SA_sources, readRDS)
VE[["baseline"]] = readRDS("output/final/M32/VE.RDS")

VE = lapply(VE, sample_n, 100)

get_SA = function(SA, parameter = "VE_BKRT"){
  
  index = which(grepl(parameter , names(SA) ))
  model_output = SA[index]
  
  out = model_output %>% 
    as.data.frame() %>%  
    mutate(ni = row_number()) %>%  
    pivot_longer(cols = - ni) %>%  
    mutate(value = value * 100) %>% 
    separate(name, into = c(NA, NA, 
                            "serostatus", "serotype",
                            "outcome", "time")) %>%  
    group_by(serotype, outcome, serostatus, ni) %>%  
    summarise(
      value = mean(value, na.rm = T)
      ) %>%  
    group_by(serotype, outcome, serostatus) %>% 
    summarise(
      lower = quantile(value, 0.025, na.rm = T),
      mean = mean(value, na.rm = T),
      upper = quantile(value, 0.975, na.rm = T) # uncertainty around time and mean 
    )
  
  return(out)
}


SA_formatted = lapply(VE,get_SA)

SA_formatted[[1]]$SA = "1"
SA_formatted[[2]]$SA = "2"
SA_formatted[[3]]$SA = "3"
SA_formatted[[4]]$SA = "4"
SA_formatted[[5]]$SA = "5"
SA_formatted[[6]]$SA = "baseline"


SA_plots = SA_formatted %>%  
  bind_rows() %>%  
  filter(serostatus != 3) %>%  
  mutate(serotype = factor(serotype, levels = 1:4,
                           labels = c("DENV1", "DENV2", "DENV3", "DENV4"))) %>% 
  mutate(serostatus = factor(serostatus, 
                             labels = c("seronegative", "monotypic"))) %>%  
  mutate(outcome = factor(outcome, levels = 1:2,
                          labels = c("symptomatic",
                                     "hospitalised"))) %>% 
  mutate(SA =factor(SA, levels = c("baseline", "1", "2", "5", "3", "4"),
                    labels = c("Main model", 
                               "1 month heterotypic\nimmunity",
                               "6 month heterotypic\nimmunity",
                               "24 month heterotypic\nimmunity",
                               "Perfect classification of\nmultitypics at baseline",
                               "No symptomatic\npost-secondary infections"))) %>% 
  ggplot(aes(y = mean, x = serotype)) +
  geom_point(aes(color = SA), position = position_dodge(0.5)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper, color = SA),
                position = position_dodge(0.5), width = 0.5, size =1) +
  facet_grid(serostatus ~ outcome, scales = "free") +
  ylab("Vaccine efficacy (%)") + xlab("Sensitivity analysis") +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() +
  theme(
    text = element_text(size = 18),
    legend.position ="top",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    legend.spacing.y = unit(0, "pt"),
    legend.margin = margin(0, 0, 0, 0)
  )    +
  scale_color_manual(values = fill) 

ggsave(
  SA_plots,
  filename = "output/figures/sens_analysis.jpg",
  width = 25,
  height = 20,
  unit = "cm"
)
