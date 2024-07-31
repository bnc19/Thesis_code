# Plot the sensitivity analysis around choice of L prior, including Bayes Factor
# support for enhancement 

setwd("Chapter2")

library(tidyverse)

# collect all BF
index_files = which(grepl("BF", list.files(path = "output/final/")))
bf_sources = paste0("output/final/", list.files(path = "output/final/")[index_files], "/BF.csv")
BF = bind_rows(lapply(bf_sources, read.csv))

ll_source = paste0("output/final/", 
                   list.files(path = "output/final/")[index_files], "/posterior.csv")

post = bind_rows(lapply(ll_source, read.csv))

ll = post %>%  
  filter(variable == "ll") %>% 
  select(mean, q5, q95)

L = post %>%  
  filter(variable == "L[1]") %>% 
  select(mean, q5, q95) %>% 
  rename(L_post_mean = mean,
         L_q5 = q5, 
         L_q95 = q95)

BF = cbind(BF, ll,L)

p1 = BF %>% 
  ggplot(aes(x = L_sd, y = BF )) +
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  ylab("Bayes factor") + xlab("") +
  scale_x_continuous(lim = c(0.25,2), breaks =seq(0.25,2, by = 0.25))

p2 = BF %>% 
  pivot_longer(cols = c(post_dens, prior_dens)) %>% 
  mutate(name = ifelse(name == "post_dens", "posterior density", "prior density")) %>% 
  ggplot(aes(x = L_sd, y =  value )) +
  geom_point(aes(color = name)) + 
  geom_line(aes(color = name)) +
  ylab("") + xlab("") +
  theme_bw()+
  theme(legend.position = c(0.85,0.8),
        legend.key.size = unit(0.2, "cm"),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
  scale_x_continuous(lim = c(0.25,2), breaks =seq(0.25,2, by = 0.25)) 
  
p3 = BF %>% 
  ggplot(aes(x =L_sd, y =  mean )) +
  geom_point() + 
  geom_line() + 
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0)) +
# geom_errorbar(aes(ymin=q5, ymax= q95), width =0.1) +  # CI don't look good 
  ylab("Log-likelihood") + xlab("")+
  scale_x_continuous(lim = c(0.25,2), breaks =seq(0.25,2, by = 0.25))

p4 = BF %>%
  ggplot(aes(x = L_sd, y =  L_post_mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = L_q5, ymax = L_q95), width = 0.05) +
  theme_bw() +
  theme(axis.text.x  = element_text(size = 14),
        axis.title.x = element_text(size = 14)) +
  ylab("L posterior \n(enhancement parameter) ") + xlab("L prior standard deviation") +
  scale_x_continuous(lim = c(0.2, 2.03), breaks = seq(0.25, 2, by = 0.25))


out = cowplot::plot_grid(p1, NULL, p2, NULL,
                         p3, NULL, p4, ncol=1, axis="tblr", align = "hv",
                         rel_heights = c(1,-0.22, 1, -0.2, 1, -0.2, 1),
                         labels = c("a", "" , "b", ' ',"c",'' ,"d"))

ggsave(
  plot = out,
  filename =  "output/figures/L_prior_sens_analysis.png",
  height = 28,
  width = 25,
  units = "cm",
  dpi = 600,
  scale = 0.7
)

