# Script to plot 2 year VE for all vaccines by serostatus and 
# serotype (but not age)

setwd("Chapter3")
rm(list = ls())

dir.create("compare_vaccines/figures")
source("compare_vaccines/R/plotting_functions.R")
library(tidyverse)
n = 500

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


TAK_file = "/Users/bethancracknelldaniels/Desktop/R projects/Thesis_code/Chapter2/output/M32/VE.RDS"
CYD_file = "CYD/output/M12/VE.RDS"
BUT_file = "BUT/output/M7/VE.RDS"

# Read in best fitting models and sample n iterations 

BUT_VE = as.data.frame(sapply(readRDS(BUT_file), sample, n))

# Select VE without age 
CYD_VE = select(readRDS(CYD_file),contains("VE_BKRT"))
TAK_VE = select(readRDS(TAK_file),contains("VE_BKRT"))

CYD_VE2 = as.data.frame(sapply(CYD_VE, sample, n))
TAK_VE2 = as.data.frame(sapply(TAK_VE, sample, n))

# summarise mean and 95% CrI VE by serostatus, serotype at 24 months
tidy_BUT_VE = format_2_year_VE(fit=BUT_VE, vaccine = "BUT")
tidy_CYD_VE = format_2_year_VE(fit=CYD_VE2, vaccine = "CYD")
tidy_TAK_VE = format_2_year_VE(fit=TAK_VE2, vaccine = "TAK")


tidy_BUT_VE$Vaccine = "Butantan-DV"
tidy_CYD_VE$Vaccine = "Dengvaxia"
tidy_TAK_VE$Vaccine = "Qdenga"

year2_comp = tidy_TAK_VE %>% 
  bind_rows(tidy_BUT_VE, tidy_CYD_VE) %>%  
  ggplot(aes(x = Vaccine, y= mean, group = Serotype)) +
  geom_point(aes(color = Serotype),
             position = position_dodge(0.4), size = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = Serotype),
                position = position_dodge(0.4), width = 0.5) +
  facet_wrap(~Serostatus, ncol =3) + 
  geom_hline(yintercept = 0, linetype = 2) +
  theme(legend.position = c(0.92,0.3),
        axis.text.x = element_text(angle = 25, hjust = 1)) +
  ylab("Vaccine efficacy (%)") +
  ggtitle("Cumulative vaccine efficacy across months 1 to 24") +
  scale_color_manual(values = serotype_fill)

ggsave(
  plot = year2_comp,
  filename = "compare_vaccines/output/year_2_VE_comp.jpg",
  height = 8,
  width = 18,
  units = "cm",
  dpi = 600
)

# summarise mean and 95% CrI VE by outcome, seerostatus, serotype upto 4.5 yrs
tidy_CYD_VE_4.5 = format_4.5_year_VE(CYD_VE2)
tidy_TAK_VE_4.5 = format_4.5_year_VE(TAK_VE2)

tidy_CYD_VE_4.5$Vaccine = "Dengvaxia"
tidy_TAK_VE_4.5$Vaccine = "Qdenga"

year_4.5_comp = tidy_TAK_VE_4.5 %>% 
  bind_rows(tidy_CYD_VE_4.5) %>%  
  ggplot(aes(x = Vaccine, y= mean, group = Serotype)) +
  geom_point(aes(color = Serotype),
             position = position_dodge(0.4)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = Serotype),
                position = position_dodge(0.4), width = 0.5) +
  facet_grid(Outcome~Serostatus) + 
  geom_hline(yintercept = 0, linetype = 2) +
  theme( legend.position ="top"  )  + 
  ylab("Vaccine efficacy (%)") +
  ggtitle("Cumulative vaccine efficacy across months 1 to 54") +
  scale_color_manual(values = serotype_fill)

ggsave(
  plot = year_4.5_comp,
  filename = "compare_vaccines/output/year_4.5_comp.jpg",
  height = 10,
  width = 18,
  units = "cm",
  dpi = 600
)

