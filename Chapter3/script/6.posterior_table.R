# Script to format the posterior parameter estimates for table

setwd("Chapter3")
library(tidyverse)

# source data 

posterior_CYD = read.csv("CYD/output/M12/posterior.csv")

tidy_post_CYD = posterior_CYD %>% 
  select(variable, mean, q5, q95) %>%  
  filter(!variable %in% c('lp__','pK2[1]','pK2[2]','pK2[3]','pK2[4]',
                          "delta[3]", "delta[4]", 
                          "epsilon", 'rho[2]','rho[3]','rho[4]',
                          'w[2]','w[3]','w[4]', "tau[2]", "tau[3]", "tau[4]",
                          "lc[1,2]", "lc[1,3]", "lc[1,4]",
                          "lc[3,1]", "lc[3,2]", "lc[3,3]","lc[3,4]",
                          "kappa",  "beta[1]", "beta[2]","beta[3]", "beta[4]",
                          "scale_FOI",
                          'alpha[2]','alpha[3]','alpha[4]',
                          "pSP[1]", "pSP[2]", "pSP[3]", "ll")) 

# only present min and max FOI and present in scientific notation 
lambda_post_CYD = tidy_post_CYD %>%
  filter(grepl("lambda", variable)) %>%
  mutate_if(is.numeric, formatC, format = "e", digits = 2)


out_CYD = tidy_post_CYD %>%  
  filter(!grepl("lambda", variable)) %>%  
  mutate_if(is.numeric, round, 2) %>%  
  mutate_if(is.numeric, as.character) %>% 
  bind_rows(lambda_post_CYD) %>% 
  unite("X", mean:q5, sep =" (") %>%  
  unite("X", X:q95, sep = " to ") %>% 
  mutate(X = paste0(X, ")")) %>% 
  separate(variable, into = c("variable", "dep", "dep2")) 

write.csv(out_CYD, "CYD/output/M12/posterior_formatted.csv")  


# Format Dengvaxia severe fit for table in thesis 

posterior_CYD_sev = read.csv("CYD/output/severe/M9/posterior.csv")

tidy_post_CYD_sev = posterior_CYD_sev %>% 
  select(variable, mean, q5, q95) %>%  
  filter(!variable %in% c('lp__','pK2[1]','pK2[2]','pK2[3]','pK2[4]',
                          "delta[3]", "delta[4]", 
                          'rho[2]','rho[3]','rho[4]',
                          'w[2]','w[3]','w[4]', "tau[2]", "tau[3]", "tau[4]",
                          "lc[1,2]", "lc[1,3]", "lc[1,4]",
                          "lc[3,1]", "lc[3,2]", "lc[3,3]","lc[3,4]",
                          "kappa", "beta[3]", "beta[4]",
                          "scale_FOI",
                          'alpha[2]','alpha[3]','alpha[4]',
                          "pSP[1]", "pSP[2]", "pSP[3]", "ll")) 


# only present min and max FOI and present in scientific notation 
lambda_post_CYD_sev = tidy_post_CYD_sev %>%
  filter(grepl("lambda", variable)) %>%
  mutate_if(is.numeric, formatC, format = "e", digits = 2)


out_CYD_sev = tidy_post_CYD_sev %>%  
  filter(!grepl("lambda", variable)) %>%  
  mutate_if(is.numeric, round, 2) %>%  
  mutate_if(is.numeric, as.character) %>% 
  bind_rows(lambda_post_CYD_sev) %>% 
  unite("X", mean:q5, sep =" (") %>%  
  unite("X", X:q95, sep = " to ") %>% 
  mutate(X = paste0(X, ")")) %>% 
  separate(variable, into = c("variable", "dep", "dep2")) 

write.csv(out_CYD_sev, "CYD/output/severe/M9/posterior_formatted.csv")  


# Format Butantan-DV posterior distribution for table in thesis  

# source data 

posterior_BUT = read.csv("BUT/output/M7/posterior.csv")

tidy_post_BUT = posterior_BUT %>% 
  select(variable, mean, q5, q95) %>%  
  filter(!variable %in% c('lp__', "epsilon", 'rho[2]',  'L[1]',
                          'L[2]', 'w[2]','w[3]', "lc[1,2]",  
                          "lc[3,1]", "lc[3,2]", "beta[1]","beta[2]", "kappa",
                          "pSP[1]", "pSP[2]", "pSP[3]", "ll")) 

# FOI  in scientific notation 
lambda_post_BUT = tidy_post_BUT %>%
  filter(grepl("lambda", variable)) %>%
  mutate_if(is.numeric, formatC, format = "e", digits = 2)

out_BUT = tidy_post_BUT %>%  
  filter(!grepl("lambda", variable)) %>%  
  mutate_if(is.numeric, round, 2) %>%  
  mutate_if(is.numeric, as.character) %>% 
  bind_rows(lambda_post_BUT) %>% 
  unite("X", mean:q5, sep =" (") %>%  
  unite("X", X:q95, sep = " to ") %>% 
  mutate(X = paste0(X, ")")) %>% 
  separate(variable, into = c("variable", "dep", "dep2")) 


write.csv(out_BUT, "BUT/output/M7/posterior_formatted.csv")  
