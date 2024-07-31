################################################################################
# Plots for manuscript 
################################################################################
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

setwd("Chapter4")
theme_set(
  theme_bw() +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      plot.margin = unit(c(0.5, 0, 0, 0), "cm"),
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      legend.margin = margin(0,0,0,0),
      legend.spacing.x = unit(0, "mm"),
      legend.spacing.y = unit(0, "mm")
    )
)


source("R/process_simulations.R")
pathout = "/outputs/processed/"

countries = "BRA"
vacc.age = 6
vacc.cov = 0.8
R0 = paste0("R0", 1:9)
serop9 = paste0("sp9 ", 1:9*10, "%")
R0_sel = R0[c(2,4,6,8)]
col =  c("#67A9CF", "#C51B8A", "#99CC99", "#FFCC99",  "#9467bd")

col_sero = c("#BDC9E1", "#D55E00", "#CC79A7", "#016C59")

# import data 
all_1000_3_max = readRDS(paste0(pathout, "serotype_dynamics_nsim_1000_3_max.RDS"))
all_1000_2_min = readRDS(paste0(pathout, "serotype_dynamics_nsim_1000_2_min.RDS"))

# tidy data --------------------------------------------------------------------

# maximum NV DENV3 scenarios 
tidy_rds_3 = all_1000_3_max %>% 
  as_tibble() %>% 
  pivot_longer(names(all_1000_3_max)[1]:names(all_1000_3_max)[length(all_1000_3_max)]) %>% 
  separate(name, into = c("country", "R0"), sep = "_") %>% 
  mutate(serop9 =  as.numeric(substr(R0, 3, nchar(R0)))) %>% 
  unnest(cols = c(value)) %>% 
  mutate(sim = factor(sim, levels= c("nv", "ni", "ns"),
                      labels = c("NV","VI", "VS")),
         incidence = "max") %>% 
  mutate(incidence_f = factor(incidence, 
                              levels = c("min", "max"),
                              labels = c("Lowest NV DENV3 attack rate", 
                                         "Highest NV DENV3 attack rate")))



# minimum NV DENV2 scenarios 
tidy_rds_2_min = all_1000_2_min %>% 
  as_tibble() %>% 
  pivot_longer(names(all_1000_2_min)[1]:names(all_1000_2_min)[length(all_1000_2_min)]) %>% 
  separate(name, into = c("country", "R0"), sep = "_") %>% 
  mutate(serop9 =  as.numeric(substr(R0, 3, nchar(R0)))) %>% 
  unnest(cols = c(value)) %>% 
  mutate(sim = factor(sim, levels= c("nv", "ni", "ns"),
                      labels = c("NV","VI", "VS")),
         incidence = "min")

# Serotype dynamic plots -------------------------------------------------------

# 1) get two example simulations from the highest DENV2/3 simulations and plot dynamics 

set.seed(10) 
samp = sample(1:1000, 1)

# D2
serotype_dyn_BRA_example_plot_min2 = tidy_rds_2_min %>% # select BRA and single max simulation 
  filter(R0 %in% R0_sel, sim != "VS", country == "BRA", 
         iteration ==  samp,  incidence == "min") %>%  
  mutate(serop9 = paste0("sp9 ", serop9 * 10, "%")) %>% # make SP 9 discrete 
  ggplot(aes(x = time, y = cases/1000)) + # thousands of cases 
  geom_line(aes(color = serotype, linetype = sim), linewidth = 1) +
  facet_wrap( ~ serop9, ncol=4) +
  ylab("Cases (thousand)") +
  xlab("Year") +
  ggtitle("Lowest NV DENV2 attack rate") +
  scale_color_manual(values = col_sero) +
  scale_x_continuous(breaks=c(10,20))+
  theme(legend.position = "top")

# D3 
serotype_dyn_BRA_example_plot_max3 = tidy_rds_3 %>% # select BRA and single max simulation 
  filter(R0 %in% R0_sel, sim != "VS", country == "BRA", 
         iteration ==  samp,  incidence == "max") %>%  
  mutate(serop9 = paste0("sp9 ", serop9*10, "%")) %>% # make SP 9 discrete 
  ggplot(aes(x = time, y = cases/1000)) + # thousands of cases 
  geom_line(aes(color = serotype, linetype = sim), linewidth = 1) +
  facet_wrap( ~ serop9, ncol =4) +
  ylab("Cases (thousand)") +
  xlab("Year") +
  ggtitle("Highest NV DENV3 attack rate") +
  scale_color_manual(values = col_sero) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks=c(10,20))


# 2) Get the % each serotype dominates by transmission setting for DENV2/3 min and max

dom_serotype_D2_min_max = tidy_rds_2_min %>%  
  group_by(country, R0, time, iteration, sim,incidence) %>% 
  mutate(max = max(cases)) %>% # at each year, which serotype has largest incidence 
  filter(cases == max) %>% 
  group_by(country, sim, serotype, serop9, incidence, R0) %>% # group by everything but time and iteration 
  summarise(dominant = n() / (20000) * 100)   # how many times does each serotype dominate at each time 

dom_serotype_D3_min_max = tidy_rds_3 %>%  
  group_by(country, R0, time, iteration, sim,incidence) %>% 
  mutate(max = max(cases)) %>% # at each year, which serotype has largest incidence 
  filter(cases == max) %>% 
  group_by(country, sim, serotype, serop9, incidence, R0) %>% # group by everything but time and iteration 
  summarise(dominant = n() / (20000) * 100)   # how many times does each serotype dominate at each time 


dom_serotype_plot_D2_min  = dom_serotype_D2_min_max %>%    
  filter(R0 %in% R0_sel, sim != "VS", country == "BRA", incidence == "min") %>% 
  mutate(serop9 = paste0("sp9 ", serop9*10, "%")) %>% # make SP 9 discrete 
  ggplot(aes(x = sim, y = dominant)) +
  geom_col(aes(fill = serotype)) + #, position = position_dodge2(0.2)
  ylab("Dominant \nserotype (%)") +
  xlab(" ")  +
  theme(legend.position = "none") +
  scale_fill_manual(values = col_sero) +
  facet_wrap(~serop9, ncol = 4)

dom_serotype_plot_D3_max  = dom_serotype_D3_min_max %>%    
  filter(R0 %in% R0_sel, sim != "VS", country == "BRA", incidence == "max") %>% 
  mutate(serop9 = paste0("sp9 ", serop9*10, "%")) %>% # make SP 9 discrete 
  ggplot(aes(x = sim, y = dominant)) +
  geom_col(aes(fill = serotype)) + #, position = position_dodge2(0.2)
  ylab("Dominant \nserotype (%)") +
  xlab(" ")  +
  theme(legend.position = "none")+
  scale_fill_manual(values = col_sero) +
  facet_wrap(~serop9, ncol = 4)


# 3) calculate how many times DENV2 serotype becomes dominant   
DENV2_outbreak = tidy_rds_2_min %>%  
  filter(incidence == "min", sim != "VS") %>% 
  group_by(country, time, iteration, sim, serop9) %>% 
  mutate(max = max(cases)) %>% # at each year, which serotype has largest incidence 
  filter(cases == max & serotype == "DENV2") %>% 
  group_by(country, iteration, sim, serop9, incidence) %>% # for each country, iteration, sim and transmission setting
  arrange(time) %>%  # arrange by time 
  mutate(outbreaks = sum(diff(time) > 1) + 1)  # how many times does serotype change over time, + 1 is number of different dominant serotypes 

# calculate how many times DENV3 is dominant 
DENV3_outbreak = tidy_rds_3 %>%  
  filter(incidence == "max", sim != "VS") %>% 
  group_by(country, time, iteration, sim, serop9) %>% 
  mutate(max = max(cases)) %>% # at each year, which serotype has largest incidence 
  filter(cases == max & serotype == "DENV3") %>% 
  group_by(country, iteration, sim, serop9) %>% # for each country, iteration, sim and transmission setting
  arrange(time) %>%  # arrange by time 
  mutate(outbreaks = sum(diff(time) > 1) +1) # how many times is DENV3 dominant, not following the year it was already dominant 


# plot dominant periods 
DENV2_outbreak_plot = DENV2_outbreak %>% 
  mutate(outbreaks = factor(outbreaks)) %>% 
  group_by(serop9,R0, country, sim, incidence, outbreaks) %>%  # over 1000 simulations, how many times does D3 dominate 
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>% 
  filter(country == "BRA", R0 %in% R0_sel) %>% 
  mutate(serop9 = paste0("sp9 ", serop9*10, "%")) %>% # make SP 9 discrete 
  ggplot(aes(x = sim, y = proportion * 100)) +
  geom_col(aes(fill = outbreaks)) + #, position = position_dodge2(0.2)
  ylab("Dominant \nDENV2 periods (%)") +
  xlab("") +
  scale_fill_manual(values=col)  +
  facet_wrap(~serop9, ncol = 4)  

DENV3_outbreak_plot = DENV3_outbreak %>% 
  group_by(serop9, country, sim) %>%  # over 1000 simulations, how many times does D3 dominate 
  mutate(outbreaks = factor(outbreaks)) %>% 
  group_by(serop9,R0, country, sim, incidence, outbreaks) %>%  # over 1000 simulations, how many times does D3 dominate 
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>% 
  filter(country == "BRA", R0 %in% R0_sel) %>% 
  mutate(serop9 = paste0("sp9 ", serop9*10, "%")) %>% # make SP 9 discrete 
  ggplot(aes(x = sim, y = proportion*100)) +
  geom_col(aes(fill = outbreaks)) + #, position = position_dodge2(0.2)
  xlab(" ") +
  ylab("Dominant \nDENV3 periods (%)") +
  scale_fill_manual(values=col)  +
  facet_wrap(~serop9, ncol = 4) 


# combine plots ----------------------------------------------------------------


DENV2_dynamics = plot_grid(serotype_dyn_BRA_example_plot_min2,
                           plot_grid(
                             dom_serotype_plot_D2_min,
                             DENV2_outbreak_plot,
                             ncol = 2, 
                             rel_widths = c(1,1),
                             labels = c("b", "c"),
                             axis = "tblr",  
                             align = "hv"
                           ), 
                           rel_heights = c(1.5,1),
                           labels = c("a", " "),
                           ncol = 1)



DENV3_dynamics = plot_grid(serotype_dyn_BRA_example_plot_max3,
                           plot_grid(dom_serotype_plot_D3_max,
                                     DENV3_outbreak_plot,
                                     ncol = 2, 
                                     rel_widths = c(1,1),
                                     labels = c("b", "c"),
                                     axis = "tblr",  
                                     align = "hv"
                           ), 
                           rel_heights = c(1.5,1),
                           labels = c("a", " "),
                           ncol = 1)





ggsave(
  filename =  paste0(save_file, "/DENV2_serotype_dynamics.png"),
  DENV2_dynamics,
  device = "png",
  width = 18,
  height = 12,
  units = "cm",
  dpi = 600
) 


ggsave(
  filename =  paste0(save_file, "/DENV3_serotype_dynamics.png"),
  DENV3_dynamics,
  device = "png",
  width = 18,
  height = 12,
  units = "cm",
  dpi = 600
) 


