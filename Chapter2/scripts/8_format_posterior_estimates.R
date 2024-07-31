# Format posterior distribution for table 
setwd("Chapter2")
library(tidyverse)

# source data 

posterior = read.csv("output/final/M32/posterior.csv")


tidy_post = posterior %>% 
  select(variable, mean, q5, q95) %>%  
  filter(!variable %in% c('lp__','pK3[1]','pK3[2]','pK3[3]','pK3[4]',
                          "epsilon", 'rho[2]','rho[3]','rho[4]',
                          'L[2]','L[3]','L[4]', 'w[2]','w[3]','w[4]',
                          "lc[1,2]", "lc[3,2]", "lc[1,3]", "lc[3,3]", 
                          "lc[1,4]", "lc[3,4]", 'w[2]','w[3]','w[4]',
                          'alpha[2]','alpha[3]','alpha[4]', "beta[2]",
                          "omega", "kappa",
                          "pSP[1]", "pSP[2]", "pSP[3]", "ll")) 

# only present min and max FOI and present in scientific notation 
lambda_post = tidy_post %>%
  filter(grepl("lambda", variable)) %>%
  filter(mean == max(mean) |
           mean == min(mean)) %>%
  mutate_if(is.numeric, formatC, format = "e", digits = 2)

out = tidy_post %>%  
  filter(!grepl("lambda", variable)) %>%  
  mutate_if(is.numeric, round, 2) %>%  
  mutate_if(is.numeric, formatC, 2, format="f") %>%  
  mutate_if(is.numeric, as.factor) %>% 
  bind_rows(lambda_post) %>% 
  unite("X", mean:q5, sep = " (") %>%  
  unite("X", X:q95, sep = " to ") %>% 
  mutate(X = paste0(X, ")")) %>% 
  separate(variable, into = c("variable", "dep", "dep2")) 


write.csv(out, "output/final/M32/posterior_formatted.csv")  


